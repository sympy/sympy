#!/usr/bin/env python

"""
Program to execute tests using the py.test like interface.

The advantage over py.test is that it only depends on sympy and should just
work in any circumstances. See "sympy.test?" for documentation.
"""

from __future__ import print_function

import sys
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import re

from get_sympy import path_hack
path_hack()

epilog = """
"options" are any of the options above.
"tests" are 0 or more glob strings of tests to run.
If no test arguments are given, all tests will be run.
"""

parser = ArgumentParser(
    epilog=epilog,
    formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-v", "--verbose", action="store_true", dest="verbose", default=False)
parser.add_argument(
    "--pdb", action="store_true", dest="pdb", default=False,
    help="Run post mortem pdb on each failure")
parser.add_argument(
    "--no-colors", action="store_false", dest="colors", default=True,
    help="Do not report colored [OK] and [FAIL]")
parser.add_argument(
    "--force-colors", action="store_true", dest="force_colors", default=False,
    help="Always use colors, even if the output is not to a terminal.")
parser.add_argument(
    "-k", dest="kw", metavar="KEYWORDS", action="store", nargs='*',
    help="Only run tests matching the given keyword expressions")
parser.add_argument(
    "--tb", dest="tb", metavar="TBSTYLE", default="short",
    help="Traceback verboseness (short/no)")
parser.add_argument(
    "--random", action="store_false", dest="sort", default=True,
    help="Run tests in random order instead of sorting them.")
parser.add_argument(
    "--seed", dest="seed", type=int, metavar="SEED",
    help="Use this seed for randomized tests.")
parser.add_argument(
    '-t', '--types', dest='types', action='store', default=None,
    choices=['gmpy', 'gmpy1', 'python'],
    help='Setup ground types.')
parser.add_argument(
    '-C', '--no-cache', dest='cache', action='store_false', default=True,
    help='Disable caching mechanism.')
parser.add_argument(
    "--timeout", action="store", dest="timeout", default=False, type=int,
    help="Set a timeout for the all functions, in seconds. "
        "By default there is no timeout.")
parser.add_argument(
    "--slow", action="store_true", dest="slow", default=False,
    help="Run only the slow functions.")
parser.add_argument(
    "--no-subprocess", action="store_false", dest="subprocess", default=True,
    help="Don't run the tests in a separate subprocess. "
        "This may prevent hash randomization from being enabled.")
parser.add_argument(
    "-E", "--enhance-asserts", action="store_true", dest="enhance_asserts",
    default=False,
    help="Rewrite assert statements to give more useful error messages.")
parser.add_argument(
    '--split', action="store", type=str, default=None,
    help="Only run part of the tests. Should be of the form a/b, (e.g., 1/2)")
parser.add_argument(
    '--rerun', action="store", dest="rerun", default=0, type=int,
    help="Number of times to rerun the specified tests.")

options, args = parser.parse_known_args()

# Check this again here to give a better error message
if options.split:
    sp = re.compile(r'([0-9]+)/([1-9][0-9]*)')
    if not sp.match(options.split):
        parser.error("option --split: must be of the form a/b where a and b "
            "are integers, not %r" % options.split)

if not options.cache:
    os.environ['SYMPY_USE_CACHE'] = 'no'
if options.types:
    os.environ['SYMPY_GROUND_TYPES'] = options.types

import sympy

ok = sympy.test(*args, verbose=options.verbose, kw=options.kw,
    tb=options.tb, pdb=options.pdb, colors=options.colors,
    force_colors=options.force_colors, sort=options.sort,
    seed=options.seed, slow=options.slow, timeout=options.timeout,
    subprocess=options.subprocess, enhance_asserts=options.enhance_asserts,
    split=options.split, rerun=options.rerun)

if ok:
    sys.exit(0)
else:
    sys.exit(1)
