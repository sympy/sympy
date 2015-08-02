from sympy import Symbol, exp, log
from sympy.calculus.singularities import singularities

from sympy.utilities.pytest import XFAIL


def test_singularities():
    x = Symbol('x', real=True)

    assert singularities(x**2, x) == ()
    assert singularities(x/(x**2 + 3*x + 2), x) == (-2, -1)


@XFAIL
def test_singularities_non_rational():
    x = Symbol('x', real=True)

    assert singularities(exp(1/x), x) == (0)
    assert singularities(log((x - 2)**2), x) == (2)


@XFAIL
def test_is_increasing():
    pass


@XFAIL
def test_is_strictly_increasing():
    pass


@XFAIL
def test_is_decreasing():
    pass


@XFAIL
def test_is_strictly_decreasing():
    pass


@XFAIL
def is_monotonic():
    pass
