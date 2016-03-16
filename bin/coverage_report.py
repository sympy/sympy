#!/usr/bin/env python
"""
Script to generate test coverage reports.

Usage:

$ bin/coverage_report.py

This will create a directory covhtml with the coverage reports. To
restrict the analysis to a directory, you just need to pass its name as
argument. For example:


$ bin/coverage_report.py sympy/logic

runs only the tests in sympy/logic/ and reports only on the modules in
sympy/logic/. To also run slow tests use --slow option. You can also get a
report on the parts of the whole sympy code covered by the tests in
sympy/logic/ by following up the previous command with


$ bin/coverage_report.py -c

"""
from __future__ import print_function

import os
import re
import sys
from optparse import OptionParser

minver = '3.4'
try:
    import coverage
    if coverage.__version__ < minver:
        raise ImportError
except ImportError:
    print(
        "You need to install module coverage (version %s or newer required).\n"
        "See http://nedbatchelder.com/code/coverage/ or \n"
        "https://launchpad.net/ubuntu/+source/python-coverage/" % minver)
    sys.exit(-1)

REPORT_DIR = "covhtml"
REFRESH = False

omit_dir_patterns = ['.*tests', 'benchmark', 'examples',
                     'pyglet', 'test_external']
omit_dir_re = re.compile(r'|'.join(omit_dir_patterns))
source_re = re.compile(r'.*\.py$')


def generate_covered_files(top_dir):
    for dirpath, dirnames, filenames in os.walk(top_dir):
        omit_dirs = [dirn for dirn in dirnames if omit_dir_re.match(dirn)]
        for x in omit_dirs:
            dirnames.remove(x)
        for filename in filenames:
            if source_re.match(filename):
                yield os.path.join(dirpath, filename)


def make_report(source_dir, report_dir, use_cache=False, slow=False):
    # code adapted from /bin/test
    bin_dir = os.path.abspath(os.path.dirname(__file__))  # bin/
    sympy_top = os.path.split(bin_dir)[0]  # ../
    sympy_dir = os.path.join(sympy_top, 'sympy')  # ../sympy/
    if os.path.isdir(sympy_dir):
        sys.path.insert(0, sympy_top)
    os.chdir(sympy_top)

    import sympy
    cov = coverage.coverage()
    cov.exclude("raise NotImplementedError")
    cov.exclude("def canonize")  # this should be "@decorated"
    cov.exclude("def __mathml__")
    if use_cache:
        cov.load()
    else:
        cov.erase()
        cov.start()
        sympy.test(source_dir, subprocess=False)
        if slow:
            sympy.test(source_dir, subprocess=False, slow=slow)
        #sympy.doctest()  # coverage doesn't play well with doctests
        cov.stop()
        cov.save()

    covered_files = list(generate_covered_files(source_dir))

    if report_dir in os.listdir(os.curdir):
        for f in os.listdir(report_dir):
            if f.split('.')[-1] in ['html', 'css', 'js']:
                os.remove(os.path.join(report_dir, f))

    cov.html_report(morfs=covered_files, directory=report_dir)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c', '--use-cache', action='store_true', default=False,
                      help='Use cached data.')
    parser.add_option('-d', '--report-dir', default='covhtml',
                      help='Directory to put the generated report in.')
    parser.add_option("--slow", action="store_true", dest="slow",
                      default=False, help="Run slow functions also.")

    options, args = parser.parse_args()

    if args:
        source_dir = args[0]
    else:
        source_dir = 'sympy/'

    make_report(source_dir, **options.__dict__)

    print("The generated coverage report is in covhtml directory.")
    print("Open %s in your web browser to view the report" % os.sep.join(
        'sympy covhtml index.html'.split()))
