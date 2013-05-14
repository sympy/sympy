from sympy import exp, sqrt, Rational
from sympy.polys import Poly
from sympy.polys.functionfields import minpoly
from sympy.polys.polyerrors import NotAlgebraic, GeneratorsError
from sympy.utilities.pytest import raises
from sympy.abc import x, y, z

def test_minpoly():
    assert minpoly(1 / x, y) == -x*y + 1
    assert minpoly(1 / (x + 1), y) == (-x - 1)*y + 1

    assert minpoly(sqrt(x), y) == y**2 - x
    assert minpoly(sqrt(x + 1), y) == y**2 - x - 1
    assert minpoly(sqrt(x) / x, y) == x*y**2 - 1
    assert minpoly(sqrt(2) * sqrt(x), y) == y**2 - 2 * x
    assert minpoly(sqrt(2) + sqrt(x), y) == y**4 + (-2*x - 4)*y**2 + x**2 - 4*x + 4

    assert minpoly(x**Rational(1,3), y) == y**3 - x
    assert minpoly(x**Rational(1,3)+sqrt(x), y) == \
        y**6 - 3*x*y**4 - 2*x*y**3 + 3*x**2*y**2 - 6*x**2*y - x**3 + x**2

    assert minpoly(sqrt(x) / z, y) == z**2*y**2 - x
    assert minpoly(sqrt(x) / (z+1), y) == (z**2 + 2*z + 1)*y**2 - x

    assert minpoly(1 / x, y, polys=True) == Poly(-x*y + 1, y)
    assert minpoly(1 / (x + 1), y, polys=True) == Poly((-x - 1)*y + 1, y)
    assert minpoly(sqrt(x), y, polys=True) == Poly(y**2 - x, y)
    assert minpoly(sqrt(x) / z, y, polys=True) == Poly(z**2*y**2 - x, y)

    raises(NotAlgebraic, lambda: minpoly(exp(x), y))
    raises(GeneratorsError, lambda: minpoly(sqrt(x), x))
