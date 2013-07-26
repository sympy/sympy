"""
This script adapts import statements in sympy/mpmath/tests to work properly.

We want to have this fully automatic, so that we don't have to do it by hand
each time we copy files from mpmath to sympy.

Usage:

$ python bin/adapt_paths.py > /tmp/x
$ patch -p0 < /tmp/x

You can use the "--dry-run" parameter to 'patch' to see if all is ok and you
should also inspect the /tmp/x that all the changes generated are actually
correct.
"""

from __future__ import print_function

from glob import glob
import re
import difflib


def get_files_mpmath():
    return glob("sympy/mpmath/tests/test_*.py")


def fix_file(filename):
    with open(filename) as f:
        orig = f.read()

    # This converts stuff like "mpmath.dps -> sympy.mpmath.dps", but will leave
    # "sympy.mpmath.dps" untouched.
    s = re.sub("(?<!sympy\.)mpmath", "sympy.mpmath", orig)

    # print differences in an unified diff format
    d = difflib.unified_diff(orig.split("\n"), s.split("\n"),
        fromfile=filename, tofile=filename + ".new", lineterm="")
    import sys
    for l in d:
        print(l)


for x in get_files_mpmath():
    fix_file(x)
