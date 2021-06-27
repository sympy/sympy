#!/usr/bin/env python3

import os
from pathlib import Path
from subprocess import check_output

def main(version, outdir):
    outdir = Path(outdir)
    build_files = [
        outdir / f'sympy-{version}.tar.gz',
        outdir / f'sympy-{version}-py3-none-any.whl',
        outdir / f'sympy-docs-html-{version}.zip',
        outdir / f'sympy-docs-pdf-{version}.pdf',
    ]
    out = check_output(['shasum', '-a', '256'] + build_files)
    out = out.decode('ascii')
    # Remove the release/ part for printing. Useful for copy-pasting into the
    # release notes.
    out = [i.split() for i in out.strip().split('\n')]
    out = '\n'.join(["%s\t%s" % (i, os.path.split(j)[1]) for i, j in out])

    # Output to file and to screen
    with open(outdir / 'sha256.txt', 'w') as shafile:
        shafile.write(out)
    print(out)


if __name__ == "__main__":
    import sys
    sys.exit(main(*sys.argv[1:]))
