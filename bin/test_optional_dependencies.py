#!/usr/bin/env python
"""
Run tests for specific packages that use optional dependencies.

The optional dependencies need to be installed before running this.
"""


import pytest

# Add the local sympy to sys.path (needed for CI)
from get_sympy import path_hack
path_hack()


class TestsFailedError(Exception):
    pass


test_list = [
    # numpy
    '*numpy*',
    'sympy/core/',
    'sympy/matrices/',
    'sympy/physics/quantum/',
    'sympy/utilities/tests/test_lambdify.py',
    'sympy/physics/control/',

    # scipy
    '*scipy*',

    # matplotlib
    'sympy/plotting/',

    # llvmlite
    '*llvm*',

    # aesara
    '*aesara*',

    # jax
    '*jax*',

    # gmpy
    'sympy/polys',

    # gmpy, numpy, scipy, autowrap, matplotlib
    'sympy/external',

    # autowrap
    '*autowrap*',

    # ipython
    '*ipython*',

    # antlr, lfortran, clang
    'sympy/parsing/',

    # codegen
    'sympy/codegen/',
    'sympy/utilities/tests/test_codegen',
    'sympy/utilities/_compilation/tests/test_compilation',
    'sympy/external/tests/test_codegen.py',

    # cloudpickle
    'pickling',

    # pycosat
    'sympy/logic',
    'sympy/assumptions',

    # stats
    'sympy/stats',

    # lxml
    "sympy/utilities/tests/test_mathml.py",
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

    # matplotlib
    'sympy/plotting/',

    # llvmlite
    '*llvm*',

    # aesara
    '*aesara*',

    # gmpy
    'sympy/polys',

    # autowrap
    '*autowrap*',

    # ipython
    '*ipython*',

    # antlr, lfortran, clang
    'sympy/parsing/',

    # codegen
    'sympy/codegen/',

    # pycosat
    'sympy/logic',
    'sympy/assumptions',

    #stats
    'sympy/stats',

    # lxml
    "sympy/utilities/mathml/",
]


print('Testing optional dependencies')


from sympy import test, doctest


tests_passed = test(*test_list, blacklist=blacklist, force_colors=True)
if tests_passed is True:
    tests_passed = pytest.ExitCode.OK
doctests_passed = doctest(*doctest_list, force_colors=True)


if (tests_passed != pytest.ExitCode.OK) and not doctests_passed:
    raise TestsFailedError('Tests and doctests failed')
elif tests_passed != pytest.ExitCode.OK:
    raise TestsFailedError('Doctests passed but tests failed')
elif not doctests_passed:
    raise TestsFailedError('Tests passed but doctests failed')
