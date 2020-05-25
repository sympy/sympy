#!/usr/bin/env python3

import subprocess
import sys
from os.path import join, splitext, basename
from contextlib import contextmanager
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from shutil import copytree


def main(sympy_doc_git, doc_html_zip, version, push=None):
    """Run this as ./update_docs.py SYMPY_DOC_GIT DOC_HTML_ZIP VERSION [--push]

    !!!!!!!!!!!!!!!!!
    NOTE: This is intended to be run as part of the release script.
    NOTE: This script will automatically push to the sympy_doc repo.
    !!!!!!!!!!!!!!!!!

    Args
    ====

    SYMPY_DOC_GIT: Path to the sympy_doc repo.
    DOC_HTML_ZIP: Path to the zip of the built html docs.
    VERSION: Version string (e.g. "1.6")
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

    update_docs(sympy_doc_git, doc_html_zip, version, push)



def update_docs(sympy_doc_git, doc_html_zip, version, push):

    def run(*cmdline):
        """Run subprocess with cwd in sympy_doc"""
        print()
        print('Running: $ ' + ' '.join(cmdline))
        print()
        return subprocess.run(cmdline, cwd=sympy_doc_git, check=True)

    @contextmanager
    def hard_reset_on_error():
        """Revert any uncomitted changes on error"""
        try:
            yield
        except Exception as e:
            run('git', 'reset', '--hard')
            raise e from None

    run('git', 'diff', '--exit-code') # Error if tree is unclean

    # We started with a clean tree so restore it on error
    with hard_reset_on_error():

        run('git', 'checkout', 'gh-pages')
        run('git', 'pull')

        update_releases_txt(sympy_doc_git, version)

        run('git', 'diff') # Show change to releases.txt
        run('git', 'add', 'releases.txt')
        run('git', 'commit', '-m', 'Add sympy %s to releases.txt' % version)

        # Delete docs for the last version
        run('git', 'rm', '-rf', 'latest')

        # Extract new docs in replacement
        extract_docs(sympy_doc_git, doc_html_zip, version)

        run('git', 'add', 'latest')
        run('git', 'add', version)
        run('git', 'commit', '-m', 'Add sympy %s docs' % version)

        # Update indexes
        run('./generate_indexes.py')
        run('git', 'diff')
        run('git', 'commit', '-a', '-m', 'Update indexes')

        run('git', 'push')


def update_releases_txt(sympy_doc_git, version):
    """Add line to the releases.txt file"""

    print()
    print("Updating releases.txt")
    print()

    releases_txt_path = join(sympy_doc_git, 'releases.txt')

    with open(releases_txt_path) as fin:
        lines = fin.readlines()

    lines += ["{0}:SymPy {0}\n".format(version)]

    with open(releases_txt_path, "w") as fout:
        fout.writelines(lines)


def extract_docs(sympy_doc_git, doc_html_zip, version):

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

        print()
        print('Copying to sympy_doc/%s' % version)
        print()
        srcpath = join(tempdir, subdirname)
        dstpath = join(sympy_doc_git, version)
        copytree(srcpath, dstpath)


if __name__ == "__main__":
    main(*sys.argv[1:])
