#!/usr/bin/env python3

import os
from subprocess import check_call

def main(version, outdir):
    check_call(['bin/mailmap_update.py'])
    check_call(['bin/authors_update.py'])
    check_call(['mkdir', '-p', outdir])
    build_release_files('bdist_wheel', 'sympy-%s-py3-none-any.whl', outdir, version)
    build_release_files('sdist', 'sympy-%s.tar.gz', outdir, version)
    check_call(['release/compare_tar_against_git.py', os.path.join(outdir, 'sympy-%s.tar.gz' % (version,)), '.'])
    check_call(['release/build_docs.py', version, outdir])


def build_release_files(cmd, fname, outdir, version):
    fname = fname % (version,)
    check_call(['python', 'setup.py', cmd])
    src = os.path.join('dist', fname)
    dst = os.path.join(outdir, fname)
    check_call(['mv', src, dst])


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
