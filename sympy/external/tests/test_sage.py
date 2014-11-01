# This testfile tests SymPy <-> Sage compatibility
#
# Execute this test inside Sage, e.g. with:
# sage -python bin/test sympy/external/tests/test_sage.py
#
# This file can be tested by Sage itself by:
# sage -t sympy/external/tests/test_sage.py
# and if all tests pass, it should be copied (verbatim) to Sage, so that it is
# automatically doctested by Sage.  Note that this second method imports the
# version of SymPy in Sage, whereas the -python method imports the local version
# of SymPy (both use the local version of the tests, however).
#
# Don't test any SymPy features here. Just pure interaction with Sage.
# Always write regular SymPy tests for anything, that can be tested in pure
# Python (without Sage). Here we test everything, that a user may need when
# using SymPy with Sage.

import os
import re
import sys

from sympy.external import import_module

sage = import_module('sage.all', __import__kwargs={'fromlist': ['all']})
if not sage:
    #bin/test will not execute any tests now
    disabled = True

if sys.version_info[0] == 3:
    # Sage does not support Python 3 currently
    disabled = True

import sympy

from sympy.utilities.pytest import XFAIL


def check_expression(expr, var_symbols):
    """Does eval(expr) both in Sage and SymPy and does other checks."""

    # evaluate the expression in the context of Sage:
    if var_symbols:
        sage.var(var_symbols)
    a = globals().copy()
    # safety checks...
    assert not "sin" in a
    a.update(sage.__dict__)
    assert "sin" in a
    e_sage = eval(expr, a)
    assert not isinstance(e_sage, sympy.Basic)

    # evaluate the expression in the context of SymPy:
    if var_symbols:
        sympy.var(var_symbols)
    b = globals().copy()
    assert not "sin" in b
    b.update(sympy.__dict__)
    assert "sin" in b
    b.update(sympy.__dict__)
    e_sympy = eval(expr, b)
    assert isinstance(e_sympy, sympy.Basic)

    # Do the actual checks:
    assert sympy.S(e_sage) == e_sympy
    assert e_sage == sage.SR(e_sympy)


def test_basics():
    check_expression("x", "x")
    check_expression("x**2", "x")
    check_expression("x**2+y**3", "x y")
    check_expression("1/(x+y)**2-x**3/4", "x y")


def test_complex():
    check_expression("I", "")
    check_expression("23+I*4", "x")


@XFAIL
def test_complex_fail():
    # Sage doesn't properly implement _sympy_ on I
    check_expression("I*y", "y")
    check_expression("x+I*y", "x y")


def test_integer():
    check_expression("4*x", "x")
    check_expression("-4*x", "x")


def test_real():
    check_expression("1.123*x", "x")
    check_expression("-18.22*x", "x")


def test_E():
    assert sympy.sympify(sage.e) == sympy.E
    assert sage.e == sage.SR(sympy.E)


def test_pi():
    assert sympy.sympify(sage.pi) == sympy.pi
    assert sage.pi == sage.SR(sympy.pi)


def test_euler_gamma():
    assert sympy.sympify(sage.euler_gamma) == sympy.EulerGamma
    assert sage.euler_gamma == sage.SR(sympy.EulerGamma)


def test_oo():
    assert sympy.sympify(sage.oo) == sympy.oo
    assert sage.oo == sage.SR(sympy.oo)
    assert sympy.sympify(-sage.oo) == -sympy.oo
    assert -sage.oo == sage.SR(-sympy.oo)


def test_NaN():
    assert sympy.sympify(sage.NaN) == sympy.nan
    assert sage.NaN == sage.SR(sympy.nan)


def test_Catalan():
    assert sympy.sympify(sage.catalan) == sympy.Catalan
    assert sage.catalan == sage.SR(sympy.Catalan)


def test_GoldenRation():
    assert sympy.sympify(sage.golden_ratio) == sympy.GoldenRatio
    assert sage.golden_ratio == sage.SR(sympy.GoldenRatio)


def test_functions():
    check_expression("sin(x)", "x")
    check_expression("cos(x)", "x")
    check_expression("tan(x)", "x")
    check_expression("cot(x)", "x")
    check_expression("asin(x)", "x")
    check_expression("acos(x)", "x")
    check_expression("atan(x)", "x")
    check_expression("atan2(y, x)", "x, y")
    check_expression("acot(x)", "x")
    check_expression("sinh(x)", "x")
    check_expression("cosh(x)", "x")
    check_expression("tanh(x)", "x")
    check_expression("coth(x)", "x")
    check_expression("asinh(x)", "x")
    check_expression("acosh(x)", "x")
    check_expression("atanh(x)", "x")
    check_expression("acoth(x)", "x")
    check_expression("exp(x)", "x")
    check_expression("log(x)", "x")
    check_expression("re(x)", "x")
    check_expression("im(x)", "x")


def test_issue_4023():
    sage.var("a x")
    log = sage.log
    i = sympy.integrate(log(x)/a, (x, a, a + 1))
    i2 = sympy.simplify(i)
    s = sage.SR(i2)
    assert s == (a*log(1 + a) - a*log(a) + log(1 + a) - 1)/a

# This string contains Sage doctests, that execute all the functions above.
# When you add a new function, please add it here as well.
"""

TESTS::

    sage: test_basics()
    sage: test_basics()
    sage: test_complex()
    sage: test_integer()
    sage: test_real()
    sage: test_E()
    sage: test_pi()
    sage: test_euler_gamma()
    sage: test_oo()
    sage: test_NaN()
    sage: test_Catalan()
    sage: test_GoldenRation()
    sage: test_functions()
    sage: test_issue_4023()

"""
