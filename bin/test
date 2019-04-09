#!/usr/bin/env python

"""
Program to execute tests using the py.test like interface.

The advantage over py.test is that it only depends on sympy and should just
work in any circumstances. See "sympy.test?" for documentation.
"""

from __future__ import print_function

import sys
import os
from optparse import OptionParser
import re

from get_sympy import path_hack
path_hack()

# callback to support variable length argument in optparse
# docs.python.org/2/library/optparse.html#callback-example-6-variable-arguments
def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []

    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)

    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)


parser = OptionParser()
parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
        default=False)
parser.add_option("--pdb", action="store_true", dest="pdb",
        default=False, help="Run post mortem pdb on each failure")
parser.add_option("--no-colors", action="store_false", dest="colors",
        default=True, help="Do not report colored [OK] and [FAIL]")
parser.add_option("--force-colors", action="store_true", dest="force_colors",
        default=False, help="Always use colors, even if the output is not to a terminal.")
parser.add_option("-k", dest="kw",
        help="only run tests matching the given keyword expressions",
        metavar="KEYWORDS", action="callback", callback=vararg_callback)
parser.add_option("--tb", dest="tb",
        help="traceback verboseness (short/no) [default: %default]",
        metavar="TBSTYLE", default="short")
parser.add_option("--random", action="store_false", dest="sort", default=True,
        help="Run tests in random order instead of sorting them")
parser.add_option("--seed", dest="seed", type="int",
        help="use this seed for randomized tests",
        metavar="SEED")
parser.add_option('-t', '--types', dest='types', action='store',
        default=None, choices=['gmpy', 'gmpy1', 'python'],
        help='setup ground types: gmpy | gmpy1 | python')
parser.add_option('-C', '--no-cache', dest='cache', action='store_false',
        default=True, help='disable caching mechanism')
parser.add_option("--timeout", action="store", dest="timeout",
        default=False, help="Set a timeout for the all functions, in seconds. By default there is no timeout.", type='int')
parser.add_option("--slow", action="store_true", dest="slow",
        default=False, help="Run only the slow functions.")
parser.add_option("--no-subprocess", action="store_false", dest="subprocess",
                  default=True, help="Don't run the tests in a separate "
                  "subprocess.  This may prevent hash randomization from being enabled.")
parser.add_option("-E", "--enhance-asserts", action="store_true", dest="enhance_asserts",
                  default=False, help="Rewrite assert statements to give more useful error messages.")
parser.add_option('--split', action="store", type='str', default=None,
    help="Only run part of the tests. Should be of the form a/b, e.g., 1/2")
parser.add_option('--rerun', action="store", dest="rerun",
                  default=0, help="Number of times to rerun the specified tests",
                  type='int')
parser.set_usage("test [options ...] [tests ...]")
parser.epilog = """\
"options" are any of the options above. "tests" are 0 or more glob strings of \
tests to run. If no test arguments are given, all tests will be run.\
"""

options, args = parser.parse_args()

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
