#! /usr/bin/env python

"""
Program to test that all methods/functions have at least one example doctest.

Usage:

bin/coverage_doctest.py sympy/core

or

bin/coverage_doctest.py sympy/core/basic.py

This script is based on the sage-coverage script from Sage written by William
Stein.
"""

import os
import re
import sys
from optparse import OptionParser

def parse_file(file, verbose=False):
    skipped = []
    missing_docstring = []
    missing_doctest = []
    has_doctest = []
    indirect_doctest = []
    while True:
        i = file.find("def ")
        if i == -1:
            break

        e = re.compile('\)\s*:')
        m = e.search(file[i:])
        if m is None:
            break
        j = m.end() + i

        # "j" now points to the end of the function definition.

        function_name = (' '.join(file[i:j].lstrip('def').split()))[:-1]
        bare_function_name = function_name[:function_name.find("(")]

        skip_this = False
        for skip in ['__dealloc__', '__new__', '_']:
            if function_name.startswith(skip + '('):
                skip_this = True
                break
        if function_name.startswith("_"):
            # For the time being, let's skip all "private" functions, that
            # beging with "_". Later, when our doctests are in a good shape, we
            # may doctest those too.
            skip_this = True
        if skip_this:
            if verbose:
                skipped.append(function_name)
            file = file[j:]
            continue

        k = file[j:].find('\n')
        if k == -1:
            break
        k += j
        kk = file[k+1:].find('\n')
        if kk == -1:
            break
        kk += k+1

        q0 = file[k:kk].find('"""')
        if q0 == -1:
            missing_docstring.append(function_name)
        else:
            q0 += k
            q1 = file[q0+3:].find('"""')
            if q1 == -1:
                print "ERROR: Error parsing %s" % function_name
            else:
                q1 += q0 + 3
                # the docstring is now between q0:q1
                d = file[q0:q1].find('>>>')
                if d == -1:
                    missing_doctest.append(function_name)
                else:
                    has_doctest.append(function_name)
                    if not (bare_function_name[0:2] == '__' and
                        bare_function_name[-2:] == '__'):
                        d = file[q0:q1].find(bare_function_name)
                        e = file[q0:q1].find('indirect doctest')
                        if d == -1 and e == -1:
                            indirect_doctest.append(function_name)

        file = file[j+3:]
    return skipped, missing_docstring, missing_doctest, has_doctest, \
            indirect_doctest


def coverage(filename, file, verbose=False):
    skipped, missing_docstring, missing_doctest, has_doctest, \
        indirect_doctest = parse_file(file, verbose)
    num_functions = len(missing_docstring + missing_doctest + has_doctest)
    if num_functions == 0:
        print "No functions in %s" % filename
        return
    print '-'*70
    print filename
    score = 100 * float(len(has_doctest)) / num_functions
    score = int(score)

    if missing_docstring:
        print "\nMissing documentation:\n\t * %s\n" % \
                ('\n\t * '.join(missing_docstring))
    if missing_doctest:
        print "\nMissing doctests:\n\t * %s\n" % \
                ('\n\t * '.join(missing_doctest))

    if indirect_doctest:
        print "\nIndirect doctest (function name doesn't occur in doctests):\n"\
                "\t * %s\n"%('\n\t * '.join(indirect_doctest))
        print 'Use "# indirect doctest" in the docstring to surpress this ' \
                'warning'

    print "SCORE %s: %s%% (%s of %s)" % (filename, score,
            len(has_doctest), num_functions)

    print '-'*70



def go(file, verbose=False, exact=True):
    if os.path.isdir(file):
        for F in os.listdir(file):
            go('%s/%s'%(file,F), verbose, exact=False)
        return
    if not (file.endswith('.py') or file.endswith('.pyx')) or \
        not exact and ('test_' in file or 'bench_' in file):
            return
    if not os.path.exists(file):
        print "File %s does not exist."%file
        sys.exit(1)
    f = open(file).read()
    coverage(file, f, verbose)

if __name__ == "__main__":
    bintest_dir = os.path.abspath(os.path.dirname(__file__))   # bin/cover...
    sympy_top  = os.path.split(bintest_dir)[0]      # ../
    sympy_dir  = os.path.join(sympy_top, 'sympy')  # ../sympy/
    if os.path.isdir(sympy_dir):
        sys.path.insert(0, sympy_top)

    parser = OptionParser()
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
            default=False)

    options, args = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
    else:
        for file in args:
            go(file, options.verbose)
