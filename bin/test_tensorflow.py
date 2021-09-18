#!/usr/bin/env python
"""
Run tests involving tensorflow

These are separate from the other optional dependency tests because tensorflow
pins the numpy version.
"""

# Add the local sympy to sys.path (needed for CI)
from get_sympy import path_hack
path_hack()


class TestsFailedError(Exception):
    pass


test_list = doctest_list = [
    'sympy/printing/tensorflow.py',
    'sympy/printing/tests/test_tensorflow.py',
    'sympy/stats/sampling',
    'sympy/utilities/lambdify.py',
    'sympy/utilities/tests/test_lambdify.py',
]


print('Testing optional dependencies')


import sympy
if not (sympy.test(*test_list, verbose=True) and sympy.doctest(*doctest_list)):
    raise TestsFailedError('Tests failed')
