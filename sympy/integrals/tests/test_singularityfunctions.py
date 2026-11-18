from __future__ import annotations
from sympy.integrals.singularityfunctions import (singularityintegrate,
    _poly_singularityintegrate)
from sympy.core.function import Function, diff
from sympy.core.numbers import Rational
from sympy.core.symbol import symbols
from sympy.functions.special.singularity_functions import SingularityFunction
from sympy.integrals.integrals import integrate

x, a, n, y = symbols('x a n y')
f = Function('f')


def test_singularityintegrate():
    assert singularityintegrate(x, x) is None
    assert singularityintegrate(x + SingularityFunction(x, 9, 1), x) is None

    assert 4*singularityintegrate(SingularityFunction(x, a, 3), x) == 4*SingularityFunction(x, a, 4)/4
    assert singularityintegrate(5*SingularityFunction(x, 5, -2), x) == 5*SingularityFunction(x, 5, -1)
    assert singularityintegrate(6*SingularityFunction(x, 5, -1), x) == 6*SingularityFunction(x, 5, 0)
    assert singularityintegrate(x*SingularityFunction(x, 0, -1), x) == 0
    assert singularityintegrate((x - 5)*SingularityFunction(x, 5, -1), x) == 0
    assert singularityintegrate(SingularityFunction(x, 0, -1) * f(x), x) == f(0) * SingularityFunction(x, 0, 0)
    assert singularityintegrate(SingularityFunction(x, 1, -1) * f(x), x) == f(1) * SingularityFunction(x, 1, 0)
    assert singularityintegrate(y*SingularityFunction(x, 0, -1)**2, x) == \
        y*SingularityFunction(0, 0, -1)*SingularityFunction(x, 0, 0)


def test_singularityintegrate_polynomial_times_singularity_function():
    # issue 29707
    assert integrate((x - 4)*SingularityFunction(x, 4, 0)/4, x) == \
        SingularityFunction(x, 4, 2)/8
    assert singularityintegrate(x**2*SingularityFunction(x, 1, 1), x) == \
        SingularityFunction(x, 1, 2)/2 + 2*SingularityFunction(x, 1, 3)/3 + \
        SingularityFunction(x, 1, 4)/4
    F = singularityintegrate((3*x + 2)*SingularityFunction(x, -2, 2), x)
    assert F == 3*SingularityFunction(x, -2, 4)/4 - \
        4*SingularityFunction(x, -2, 3)/3
    for p in (-3, -1, Rational(1, 2), 5):
        assert diff(F, x).subs(x, p) == \
            ((3*x + 2)*SingularityFunction(x, -2, 2)).subs(x, p)
    # left to the Heaviside rewrite: an exponent that may be negative, a
    # second SingularityFunction or a factor that is not a polynomial
    assert _poly_singularityintegrate(x*SingularityFunction(x, 1, n), x) is None
    assert _poly_singularityintegrate(
        SingularityFunction(x, 1, 0)*SingularityFunction(x, 2, 0), x) is None
    assert _poly_singularityintegrate(f(x)*SingularityFunction(x, 1, 0), x) is None
