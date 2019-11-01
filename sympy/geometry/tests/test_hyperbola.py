from sympy import Rational, S, Symbol, symbols, pi, sqrt, oo, Point2D, Segment2D, Abs
from sympy.core.compatibility import range
from sympy.geometry import (Circle, Ellipse, GeometryError, Line, Point,
                            Polygon, Ray, RegularPolygon, Segment,
                            Triangle, intersection)
from sympy.geometry.hyperbola import Hyperbola
from sympy.utilities.pytest import raises, slow
from sympy import integrate
from sympy.functions.special.elliptic_integrals import elliptic_e
from sympy.functions.elementary.miscellaneous import Max, Min


def test_hyperbola_geom():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    t = Symbol('t', real=True)
    a = Symbol('a', real=True)
    b = Symbol('b', real=True)
    m = Symbol('m', real=True)
    p1 = Point(0, 0)
    p2 = Point(1, 1)
    p3 = Point(2, 3)
    p4 = Point(0, 1)
    h1 = Hyperbola(p1, 2, 1)
    h2 = Hyperbola(p2, 3, 4)

    assert Hyperbola(p3, 3, 4) == Hyperbola(Point2D(2, 3), 3, 4)

    raises(ValueError, lambda: Hyperbola(None, None, None, 1))

    assert h1 != h2

    #Center
    assert Hyperbola(None, 1, 1).center == p1
    assert h1.center == p1
    assert h2.center == p2

    #hradius and vradius
    assert h1.hradius == 2
    assert h1.vradius == 1

    #eccentricity
    assert h1.eccentricity == sqrt(5) / 2
    assert h2.eccentricity == Rational(5, 3)

    #foci
    assert h1.foci == (Point2D(-sqrt(5), 0), Point2D(sqrt(5), 0))
    assert h2.foci == (Point2D(-4, 1), Point2D(6, 1))

    #major and minor
    assert h1.major == 2
    assert h1.minor == 1
    assert Hyperbola(p1, a, b).major == a
    assert Hyperbola(p1, a, b).minor == b
    assert Hyperbola(p2, m, m+1).major == m+1
    assert Hyperbola(p2, m, m+1).minor == m

    #focus distance
    assert h1.focus_distance == sqrt(5)
    assert h2.focus_distance == 5

    #semi-latusrectum
    assert h1.semilatus_rectum == Rational(1, 2)
    assert h2.semilatus_rectum == Rational(64, 9)

