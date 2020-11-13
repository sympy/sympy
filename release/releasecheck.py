#!/usr/bin/env python3

from os.path import join, basename, normpath
from subprocess import check_call

def main(version, outdir):
    check_version(version, outdir)
    check_call(['bin/mailmap_update.py'])
    check_call(['bin/authors_update.py'])
    check_call(['mkdir', '-p', outdir])
    build_release_files('bdist_wheel', 'sympy-%s-py3-none-any.whl', outdir, version)
    build_release_files('sdist', 'sympy-%s.tar.gz', outdir, version)
    check_call(['release/compare_tar_against_git.py', join(outdir, 'sympy-%s.tar.gz' % (version,)), '.'])
    check_call(['release/test_install.py', version, outdir])
    check_call(['release/build_docs.py', version, outdir])


def build_release_files(cmd, fname, outdir, version):
    fname = fname % (version,)
    check_call(['python', 'setup.py', '-q', cmd])
    src = join('dist', fname)
    dst = join(outdir, fname)
    check_call(['mv', src, dst])


def check_version(version, outdir):
    from sympy.release import __version__ as checked_out_version
    if version != checked_out_version:
        msg = "version %s does not match checkout %s"
        raise AssertionError(msg % (version, checked_out_version))
    if basename(normpath(outdir)) != 'release-%s' % (version,):
        msg = "version %s does not match output directory %s"
        raise AssertionError(msg % (version, outdir))


if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
