#!/usr/bin/env python
"""
Run tests for specific packages that use optional dependencies.

The optional dependencies need to be installed before running this.
"""


# Add the local sympy to sys.path (needed for CI)
from get_sympy import path_hack
path_hack()


class TestsFailedError(Exception):
    pass

print('Testing optional dependencies')

import sympy
test_list = [
    # numpy
    '*numpy*',
    'sympy/core/',
    'sympy/matrices/',
    'sympy/physics/quantum/',
    'sympy/utilities/tests/test_lambdify.py',

    # scipy
    '*scipy*',

    # llvmlite
    '*llvm*',

    # theano
    '*theano*',

    # gmpy
    'polys',

    # autowrap
    '*autowrap*',

    # ipython
    '*ipython*',

    # antlr, lfortran, clang
    'sympy/parsing/',

    # matchpy
    '*rubi*',

    # codegen
    'sympy/codegen/',
    'sympy/utilities/tests/test_codegen',
    'sympy/utilities/_compilation/tests/test_compilation',

    # cloudpickle
    'pickling',

    # pycosat
    'sympy/logic',
    'sympy/assumptions',

    #stats
    'sympy/stats',

]

blacklist = [
    'sympy/physics/quantum/tests/test_circuitplot.py',
]

doctest_list = [
    # numpy
    'sympy/matrices/',
    'sympy/utilities/lambdify.py',

    # scipy
    '*scipy*',

    # llvmlite
    '*llvm*',

    # theano
    '*theano*',

    # gmpy
    'polys',

    # autowrap
    '*autowrap*',

    # ipython
    '*ipython*',

    # antlr, lfortran, clang
    'sympy/parsing/',

    # matchpy
    '*rubi*',

    # codegen
    'sympy/codegen/',

    # pycosat
    'sympy/logic',
    'sympy/assumptions',

    #stats
    'sympy/stats',

]

if not (sympy.test(*test_list, blacklist=blacklist) and sympy.doctest(*doctest_list)):
    raise TestsFailedError('Tests failed')


print('Testing MATPLOTLIB')
# Set matplotlib so that it works correctly in headless Travis. We have to do
# this here because it doesn't work after the sympy plotting module is
# imported.
import matplotlib
matplotlib.use("Agg")
import sympy
# Unfortunately, we have to use subprocess=False so that the above will be
# applied, so no hash randomization here.
if not (sympy.test('sympy/plotting', 'sympy/physics/quantum/tests/test_circuitplot.py',
    subprocess=False) and sympy.doctest('sympy/plotting', subprocess=False)):
    raise TestsFailedError('Tests failed')


print('Testing SYMENGINE')
import sympy
if not sympy.test('sympy/physics/mechanics'):
    raise TestsFailedError('Tests failed')
if not sympy.test('sympy/liealgebras'):
    raise TestsFailedError('Tests failed')
