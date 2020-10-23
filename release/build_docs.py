#!/usr/bin/env python3


import os
from os.path import dirname, join
from os import chdir

from helpers import run


ROOTDIR = dirname(dirname(__file__))


def main(version):
    check_version(version)
    build_html(ROOTDIR)
    build_latex(ROOTDIR)


def build_html(rootdir):
    docsdir = join(rootdir, 'doc')
    run('make', 'clean', cwd=docsdir)
    run('make', 'html', cwd=docsdir)


def build_latex(rootdir):
    docsdir = join(rootdir, 'doc')
    latexdir = join(docsdir, '_build', 'latex')
    os.environ['LATEXMKOPTS'] = '-xelatex -silent'
    run('make', 'clean', cwd=docsdir)
    run('make', 'latex', cwd=docsdir)
    run('make', 'clean', cwd=latexdir)
    run('make', 'all', cwd=latexdir)


def check_version(version):
    from sympy.release import __version__ as checked_out_version
    if version != checked_out_version:
        msg = "version %s does not match checkout %s"
        raise AssertionError(msg % (version, checked_out_version))


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
