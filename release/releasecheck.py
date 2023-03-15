#!/usr/bin/env python3

from os.path import join, basename, normpath
from subprocess import check_call

def main(version, prevversion, outdir):
    check_version(version, outdir)
    run_stage(['bin/mailmap_check.py', '--update-authors'])
    run_stage(['mkdir', '-p', outdir])
    build_release_files('--wheel', 'sympy-%s-py3-none-any.whl', outdir, version)
    build_release_files('--sdist', 'sympy-%s.tar.gz', outdir, version)
    run_stage(['release/compare_tar_against_git.py', join(outdir, 'sympy-%s.tar.gz' % (version,)), '.'])
    run_stage(['release/build_docs.py', version, outdir])
    run_stage(['release/sha256.py', version, outdir])
    run_stage(['release/authors.py', version, prevversion, outdir])


def green(text):
    return "\033[32m%s\033[0m" % text

def red(text):
    return "\033[31m%s\033[0m" % text

def print_header(color, *msgs):
    newlines = '\n'
    vline = '-' * 80
    print(color(newlines + vline))
    for msg in msgs:
        print(color(msg))
    print(color(vline + newlines))

def run_stage(cmd):
    cmdline = '    $ %s' % (' '.join(cmd),)

    print_header(green, 'running:', cmdline)
    try:
        check_call(cmd)
    except Exception as e:
        print_header(red, 'failed:', cmdline)
        raise e from None
    else:
        print_header(green, 'completed:', cmdline)


def build_release_files(cmd, fname, outdir, version):
    fname = fname % (version,)
    run_stage(['python', '-m', 'build', cmd])
    src = join('dist', fname)
    dst = join(outdir, fname)
    run_stage(['mv', src, dst])


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
