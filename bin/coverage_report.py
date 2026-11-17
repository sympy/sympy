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
from argparse import ArgumentParser

minver = '3.4'
try:
    import coverage
    if coverage.__version__ < minver:
        raise ImportError
except ImportError:
    print(
        "You need to install module coverage (version %s or newer required).\n"
        "See https://coverage.readthedocs.io/en/latest/ or \n"
        "https://launchpad.net/ubuntu/+source/python-coverage/" % minver)
    sys.exit(-1)


omit_dir_patterns = ['benchmark', 'examples',
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


def make_report(
    test_args, source_dir='sympy/', report_dir='covhtml', use_cache=False,
    slow=False
    ):
    # code adapted from /bin/test
    from get_sympy import path_hack
    sympy_top = path_hack()
    os.chdir(sympy_top)

    cov = coverage.coverage()
    cov.exclude("raise NotImplementedError")
    cov.exclude("def canonize")  # this should be "@decorated"
    if use_cache:
        cov.load()
    else:
        cov.erase()
        cov.start()
        import sympy
        sympy.test(*test_args, subprocess=False, slow=slow)
        #sympy.doctest()  # coverage doesn't play well with doctests
        cov.stop()
        try:
            cov.save()
        except PermissionError:
            import warnings
            warnings.warn(
                "PermissionError has been raised while saving the " \
                "coverage result.",
                RuntimeWarning
            )

    covered_files = list(generate_covered_files(source_dir))
    cov.html_report(morfs=covered_files, directory=report_dir)

parser = ArgumentParser()
parser.add_argument(
    '-c', '--use-cache', action='store_true', default=False,
    help='Use cached data.')
parser.add_argument(
    '-d', '--report-dir', default='covhtml',
    help='Directory to put the generated report in.')
parser.add_argument(
    "--slow", action="store_true", dest="slow", default=False,
    help="Run slow functions also.")
options, args = parser.parse_known_args()

if __name__ == '__main__':
    report_dir = options.report_dir
    use_cache = options.use_cache
    slow = options.slow
    make_report(
        args, report_dir=report_dir, use_cache=use_cache, slow=slow)

    print("The generated coverage report is in covhtml directory.")
    print(
        "Open %s in your web browser to view the report" %
        os.sep.join([report_dir, 'index.html'])
    )
