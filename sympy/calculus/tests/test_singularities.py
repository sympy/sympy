from sympy import Symbol, exp, log
from sympy.calculus.singularities import singularities

from sympy.utilities.pytest import XFAIL


def test_singularities():
    x = Symbol('x', real=True)

    assert singularities(x**2, x) == ()
    assert singularities(x/(x**2 + 3*x + 2), x) == (-2, -1)
    assert singularities((x**2 - 1)/(x**3 - 1), x) == (1,)


@XFAIL
def test_singularities_non_rational():
    x = Symbol('x', real=True)

    assert singularities(exp(1/x), x) == (0)
    assert singularities(log((x - 2)**2), x) == (2)
