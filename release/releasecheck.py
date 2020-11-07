#!/usr/bin/env python3

from subprocess import check_call

def main(version, outdir):
    check_call(['bin/mailmap_update.py'])
    check_call(['bin/authors_update.py'])
    check_call(['release/build_docs.py', version, outdir])

if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
