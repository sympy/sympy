"""Tests for tools for manipulation of rational expressions. """

from sympy.polys.rationaltools import together
from sympy.polys.rationaltools import thiele_interpolate as thiele

from sympy.core.mul import Mul
from sympy.core.numbers import Rational
from sympy.core.relational import Eq
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.trigonometric import sin
from sympy.integrals.integrals import Integral
from sympy.abc import x, y, z

A, B = symbols('A,B', commutative=False)


def test_together():
    assert together(0) == 0
    assert together(1) == 1

    assert together(x*y*z) == x*y*z
    assert together(x + y) == x + y

    assert together(1/x) == 1/x

    assert together(1/x + 1) == (x + 1)/x
    assert together(1/x + 3) == (3*x + 1)/x
    assert together(1/x + x) == (x**2 + 1)/x

    assert together(1/x + S.Half) == (x + 2)/(2*x)
    assert together(S.Half + x/2) == Mul(S.Half, x + 1, evaluate=False)

    assert together(1/x + 2/y) == (2*x + y)/(y*x)
    assert together(1/(1 + 1/x)) == x/(1 + x)
    assert together(x/(1 + 1/x)) == x**2/(1 + x)

    assert together(1/x + 1/y + 1/z) == (x*y + x*z + y*z)/(x*y*z)
    assert together(1/(1 + x + 1/y + 1/z)) == y*z/(y + z + y*z + x*y*z)

    assert together(1/(x*y) + 1/(x*y)**2) == y**(-2)*x**(-2)*(1 + x*y)
    assert together(1/(x*y) + 1/(x*y)**4) == y**(-4)*x**(-4)*(1 + x**3*y**3)
    assert together(1/(x**7*y) + 1/(x*y)**4) == y**(-4)*x**(-7)*(x**3 + y**3)

    assert together(5/(2 + 6/(3 + 7/(4 + 8/(5 + 9/x))))) == \
        Rational(5, 2)*((171 + 119*x)/(279 + 203*x))

    assert together(1 + 1/(x + 1)**2) == (1 + (x + 1)**2)/(x + 1)**2
    assert together(1 + 1/(x*(1 + x))) == (1 + x*(1 + x))/(x*(1 + x))
    assert together(
        1/(x*(x + 1)) + 1/(x*(x + 2))) == (3 + 2*x)/(x*(1 + x)*(2 + x))
    assert together(1 + 1/(2*x + 2)**2) == (4*(x + 1)**2 + 1)/(4*(x + 1)**2)

    assert together(sin(1/x + 1/y)) == sin(1/x + 1/y)
    assert together(sin(1/x + 1/y), deep=True) == sin((x + y)/(x*y))

    assert together(1/exp(x) + 1/(x*exp(x))) == (1 + x)/(x*exp(x))
    assert together(1/exp(2*x) + 1/(x*exp(3*x))) == (1 + exp(x)*x)/(x*exp(3*x))

    assert together(Integral(1/x + 1/y, x)) == Integral((x + y)/(x*y), x)
    assert together(Eq(1/x + 1/y, 1 + 1/z)) == Eq((x + y)/(x*y), (z + 1)/z)

    assert together((A*B)**-1 + (B*A)**-1) == (A*B)**-1 + (B*A)**-1


def test_thiele_interpolate():
    assert thiele([1, 2], [0, 0]) == 0
    assert thiele([2, 1], [0, 0]) == 0
    assert thiele([1, 2], [1, 1]) == 1
    assert thiele([2, 1], [1, 1]) == 1
    assert thiele([1, 2], [1, 2]) == x
    assert thiele([2, 1], [1, 2]) == 3 - x
    assert thiele([1, 2, 3], [1, 2, 3]) == x
    assert thiele([1, 2, 3], [3, 2, 1]) == 4 - x
    assert thiele([1, 2, 3], [2, 4, 6]) == 2*x
    assert thiele([2, 1, 3], [4, 2, 6]) == 2*x
    assert thiele([1, 2], [1, 4]) == 3*x - 2
    assert thiele([1, 2, 3, 4], [1, 4, 9, 16]) == x**2
    assert thiele([1, 2, 3, 4, 5], [1, 4, 9, 16, 25]) == x**2
    assert thiele([1, 2, 3, 4, 5, 6], [1, 8, 27, 64, 125, 216]) == x**3
    assert thiele([1, 2, 3, 4], [1, 2, 1, 2]) == (x**2 - 2*x - 2)/(2*x - 5)
    assert thiele([1, 2, 4, 9], [0, 1, 0, 1]) == (x**2 - 5*x + 4)/(6*x - 14)
