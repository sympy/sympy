from sympy.core.numbers import Rational
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import cos
from sympy.geometry.point import Point2D
from sympy.geometry import (Circle, GeometryError, Hyperbola, Point)
from sympy.testing.pytest import raises
from sympy.functions.elementary.miscellaneous import Max


def test_hyperbola_equation_using_slope():
    from sympy.abc import x, y

    h1 = Hyperbola(Point(1, 0), 3, 2)
    assert str(h1.equation(_slope=1)) == str(-(-x + y + 1)**2/8 + (x + y - 1)**2/18 - 1)

    h2 = Hyperbola(Point(0, 0), 4, 1)
    assert str(h2.equation(_slope=1)) == str(-(-x + y)**2/2 + (x + y)**2/32 - 1)

    h3 = Hyperbola(Point(1, 5), 6, 2)
    assert str(h3.equation(_slope=2)) == str(-(-2*x + y - 3)**2/20 + (x + 2*y - 11)**2/180 - 1)


def test_construction():
    h1 = Hyperbola(hradius=2, vradius=1, eccentricity=None)
    assert h1.eccentricity == sqrt(5)/2

    h2 = Hyperbola(hradius=2, vradius=None, eccentricity=sqrt(5)/2)
    assert h2.vradius == 1

    h3 = Hyperbola(hradius=None, vradius=1, eccentricity=sqrt(5)/2)
    assert h3.hradius == 2

    #tests for eccentricity < 1
    raises(GeometryError, lambda: Hyperbola(Point(3, 1), hradius=3, eccentricity = S(1)/2))
    raises(GeometryError, lambda: Hyperbola(Point(3, 1), hradius=3, eccentricity = cos(5)))
    raises(GeometryError, lambda: Hyperbola(Point(3, 1), hradius=3, eccentricity = S.Pi - S(3)))
    raises(GeometryError, lambda: Hyperbola(Point(3, 1), hradius=3, eccentricity = -3))

    #tests for eccentricity = 1
    #if hradius, vradius is not defined
    raises(GeometryError, lambda: Hyperbola(None, None, 1, eccentricity = 1))


def test_auxiliary_circle():
    x, y, a, b = symbols('x y a b')
    h = Hyperbola((x, y), a, b)
    # the general result
    assert h.auxiliary_circle() == Circle((x, y), Max(a, b))


def test_director_circle():
    x, y, a, b = symbols('x y a b')
    h = Hyperbola((x, y), a, b)
    # the general result
    assert h.director_circle() == Circle(Point2D(x, y), sqrt(a**2 - b**2))


def test_evolute():
    #Hyperbola centered at h, k
    x, y, h, k = symbols('x y h k',real = True)
    a, b = symbols('a b')
    h1 = Hyperbola(Point(h, k), a, b)
    t1 = (h1.hradius*(x - h1.center.x))**Rational(2, 3)
    t2 = (h1.vradius*(y - h1.center.y))**Rational(2, 3)
    H = t1 - t2 - (h1.hradius**2 + h1.vradius**2)**Rational(2, 3)
    assert h1.evolute() == H
    #Numerical Example
    h1 = Hyperbola(Point(1, 1), 6, 3)
    t1 = (6*(x - 1))**Rational(2, 3)
    t2 = (3*(y - 1))**Rational(2, 3)
    H = t1 - t2 - (45)**Rational(2, 3)
    assert h1.evolute() == H