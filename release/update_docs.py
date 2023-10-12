#!/usr/bin/env python3

import json
import subprocess
import sys
from os.path import join, splitext, basename
from contextlib import contextmanager
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from shutil import copytree



def main(sympy_doc_git, doc_html_zip, version, dev_version, push=None):
    """Run this as ./update_docs.py SYMPY_DOC_GIT DOC_HTML_ZIP VERSION [--push]

    !!!!!!!!!!!!!!!!!
    NOTE: This is intended to be run as part of the release script.
    NOTE: This script will automatically push to the sympy_doc repo.
    !!!!!!!!!!!!!!!!!

    Args
    ====

    SYMPY_DOC_GIT: Path to the sympy_doc repo.
    DOC_HTML_ZIP: Path to the zip of the built html docs.
    VERSION: Version string of the release (e.g. "1.6")
    DEV_VERSION: Version string of the development version (e.g. "1.7.dev")
    --push (optional): Push the results (Warning this pushes direct to github)

    This script automates the "release docs" step described in the README of the
    sympy/sympy_doc repo:

    https://github.com/sympy/sympy_doc#release-docs
    """
    if push is None:
        push = False
    elif push == "--push":
        push = True
    else:
        raise ValueError("Invalid arguments")

    update_docs(sympy_doc_git, doc_html_zip, version, dev_version, push)


def update_docs(sympy_doc_git, doc_html_zip, version, dev_version, push):

    # We started with a clean tree so restore it on error
    with git_rollback_on_error(sympy_doc_git, branch='gh-pages') as run:

        # Delete docs for the last version
        run('git', 'rm', '-rf', 'latest')

        # Extract new docs in replacement
        extract_docs(sympy_doc_git, doc_html_zip)

        # Commit new docs
        run('git', 'add', 'latest')
        run('git', 'commit', '-m', 'Add sympy %s docs' % version)

        # Update versions.json
        with open(join(sympy_doc_git, 'versions.json'), 'w') as f:
            json.dump({'latest': version, 'dev': dev_version}, f)
        run('git', 'diff')
        run('git', 'add', 'versions.json')
        run('git', 'commit', '-m', 'Update versions.json')

        if push:
            run('git', 'push')
        else:
            print('Results are committed but not pushed')


@contextmanager
def git_rollback_on_error(gitroot_path, branch='master'):

    def run(*cmdline, **kwargs):
        """Run subprocess with cwd in sympy_doc"""
        print()
        print('Running: $ ' + ' '.join(cmdline))
        print()
        return subprocess.run(cmdline, cwd=gitroot_path, check=True, **kwargs)

    unclean_msg = "The git repo should be completely clean before running this"

    try:
        run('git', 'diff', '--exit-code') # Error if tree is unclean
    except subprocess.CalledProcessError:
        raise ValueError(unclean_msg)
    if run('git', 'clean', '-n', stdout=subprocess.PIPE).stdout:
        raise ValueError(unclean_msg)

    run('git', 'checkout', branch)
    run('git', 'pull')

    bsha_start = run('git', 'rev-parse', 'HEAD', stdout=subprocess.PIPE).stdout
    sha_start = bsha_start.strip().decode('ascii')

    try:
        yield run
    except Exception as e:
        run('git', 'reset', '--hard', sha_start)
        raise e from None

def extract_docs(sympy_doc_git, doc_html_zip):

    subdirname = splitext(basename(doc_html_zip))[0]

    with TemporaryDirectory() as tempdir:
        print()
        print('Extracting docs to ' + tempdir)
        print()
        ZipFile(doc_html_zip).extractall(tempdir)

        print()
        print('Copying to sympy_doc/latest')
        print()
        srcpath = join(tempdir, subdirname)
        dstpath = join(sympy_doc_git, 'latest')
        copytree(srcpath, dstpath)

if __name__ == "__main__":
    main(*sys.argv[1:])
