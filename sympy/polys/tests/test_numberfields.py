"""Tests for computational algebraic number field theory. """

from sympy import Poly, raises, sin, sqrt, I

from sympy.polys.numberfields import (
    minpoly, AlgebraicNumber,
)

from sympy.polys.polyerrors import (
    NotAlgebraic,
)

from sympy.abc import x, y, a

def test_minpoly():
    assert minpoly(-7, x) == x + 7
    assert minpoly(-1, x) == x + 1
    assert minpoly( 0, x) == x
    assert minpoly( 1, x) == x - 1
    assert minpoly( 7, x) == x - 7

    assert minpoly(sqrt(2), x) == x**2 - 2
    assert minpoly(sqrt(5), x) == x**2 - 5
    assert minpoly(sqrt(6), x) == x**2 - 6

    assert minpoly(2*sqrt(2), x) == x**2 - 8
    assert minpoly(3*sqrt(5), x) == x**2 - 45
    assert minpoly(4*sqrt(6), x) == x**2 - 96

    assert minpoly(2*sqrt(2) + 3, x) == x**2 -  6*x +  1
    assert minpoly(3*sqrt(5) + 6, x) == x**2 - 12*x -  9
    assert minpoly(4*sqrt(6) + 7, x) == x**2 - 14*x - 47

    assert minpoly(2*sqrt(2) - 3, x) == x**2 +  6*x +  1
    assert minpoly(3*sqrt(5) - 6, x) == x**2 + 12*x -  9
    assert minpoly(4*sqrt(6) - 7, x) == x**2 + 14*x - 47

    assert minpoly(sqrt(1 + sqrt(6)), x) == x**4 -  2*x**2 -  5
    assert minpoly(sqrt(I + sqrt(6)), x) == x**8 - 10*x**4 + 49

    assert minpoly(2*I + sqrt(2 + I), x) == x**4 + 4*x**2 + 8*x + 37

    assert minpoly(sqrt(2) + sqrt(3), x) == x**4 - 10*x**2 + 1
    assert minpoly(sqrt(2) + sqrt(3) + sqrt(6), x) == x**4 - 22*x**2 - 48*x - 23

    assert minpoly(1/(1 - 9*sqrt(2) + 7*sqrt(3)), x) == 392*x**4 - 1232*x**3 + 612*x**2 + 4*x - 1

    raises(NotAlgebraic, "minpoly(a, x)")
    raises(NotAlgebraic, "minpoly(2**y, x)")
    raises(NotAlgebraic, "minpoly(sin(1), x)")

    assert minpoly(sqrt(2), polys=True).is_Poly == True
    assert minpoly(sqrt(2), x, polys=True) == Poly(x**2 - 2)

