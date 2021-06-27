#!/usr/bin/env python3


import os
from os.path import dirname, join, basename, normpath
from os import chdir
import shutil

from helpers import run


ROOTDIR = dirname(dirname(__file__))
DOCSDIR = join(ROOTDIR, 'doc')


def main(version, outputdir):
    os.makedirs(outputdir, exist_ok=True)
    build_html(DOCSDIR, outputdir, version)
    build_latex(DOCSDIR, outputdir, version)


def build_html(docsdir, outputdir, version):
    run('make', 'clean', cwd=docsdir)
    run('make', 'html', cwd=docsdir)

    builddir = join(docsdir, '_build')
    docsname = 'sympy-docs-html-%s' % (version,)
    zipname = docsname + '.zip'
    cwd = os.getcwd()
    try:
        chdir(builddir)
        shutil.move('html', docsname)
        run('zip', '-9lr', zipname, docsname)
    finally:
        chdir(cwd)
    shutil.move(join(builddir, zipname), join(outputdir, zipname))


def build_latex(docsdir, outputdir, version):
    run('make', 'clean', cwd=docsdir)
    run('make', 'latex', cwd=docsdir)

    latexdir = join(docsdir, '_build', 'latex')
    env = os.environ.copy()
    env['LATEXMKOPTS'] = '-xelatex -silent'
    run('make', 'clean', cwd=latexdir, env=env)
    run('make', 'all', cwd=latexdir, env=env)

    srcfilename = 'sympy-%s.pdf' % (version,)
    dstfilename = 'sympy-docs-pdf-%s.pdf' % (version,)
    src = join('doc', '_build', 'latex', srcfilename)
    dst = join(outputdir, dstfilename)
    shutil.copyfile(src, dst)


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
