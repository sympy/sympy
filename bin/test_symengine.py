#!/usr/bin/env python
"""
Run tests involving symengine

These are separate from the other optional dependency tests because they need
to be run with the USE_SYMENGINE=1 environment variable set. This script does
not set the environment variable by default so that the same tests can be run
with and without symengine.

Run this as:

    $ USE_SYMENGINE=1 bin/test_symengine.py
"""

# Add the local sympy to sys.path (needed for CI)
from get_sympy import path_hack
path_hack()


class TestsFailedError(Exception):
    pass


test_list = [
    'sympy/physics/mechanics',
    'sympy/liealgebras',
]


print('Testing optional dependencies')

#
# XXX: The doctests are not tested here but there are many failures when
# running them with symengine.
#
import sympy
if not sympy.test(*test_list, verbose=True):
    raise TestsFailedError('Tests failed')
