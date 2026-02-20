from sympy.core.numbers import Rational
from sympy.core.singleton import S
from sympy.core.symbol import symbols, Symbol
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import cos
from sympy.geometry.point import Point2D
from sympy.geometry import (Circle, Ellipse, GeometryError, Hyperbola, Line, Point,
                            Polygon, Ray, RegularPolygon, Segment, intersection)
from sympy.testing.pytest import raises, slow
from sympy.functions.elementary.miscellaneous import Max


def test_hyperbola_equation_using_slope_and_branch():
    from sympy.abc import x, y

    h1 = Hyperbola(Point(1, 0), 3, 2)
    assert str(h1.equation()) == str((x/3 - Rational(1, 3))**2 - y**2/4 - 1)
    assert str(h1.equation(branch='left')) == str(x + 3*sqrt(y**2 + 4)/2 - 1)
    assert str(h1.equation(branch='right')) == str(x - 3*sqrt(y**2 + 4)/2 - 1)
    assert str(h1.equation(_slope=1)) == str(-(-x + y + 1)**2/8 + (x + y - 1)**2/18 - 1)

    h2 = Hyperbola(Point(0, 0), 4, 1)
    assert str(h2.equation()) == str(x**2/16 - y**2 - 1)
    assert str(h2.equation(branch='left')) == str(x + 4*sqrt(y**2 + 1))
    assert str(h2.equation(branch='right')) == str(x - 4*sqrt(y**2 + 1))
    assert str(h2.equation(_slope=1)) == str(-(-x + y)**2/2 + (x + y)**2/32 - 1)

    h3 = Hyperbola(Point(1, 5), 6, 2)
    assert str(h3.equation()) == str((x/6 - Rational(1, 6))**2 - (y/2 - Rational(5, 2))**2 - 1)
    assert str(h3.equation(branch='left')) == str(x + 3*sqrt((y - 5)**2 + 4) - 1)
    assert str(h3.equation(branch='right')) == str(x - 3*sqrt((y - 5)**2 + 4) - 1)
    assert str(h3.equation(_slope=2)) == str(-(-2*x + y - 3)**2/20 + (x + 2*y - 11)**2/180 - 1)


@slow
def test_hyperbola_geom():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    t = Symbol('t', real=True)
    a = Symbol('a', real=True)
    b = Symbol('b', real=True)
    half = S.Half
    p1 = Point(0, 0)
    p2 = Point(1, 1)
    p3 = Point(1, 0)
    h1 = Hyperbola(p1, 1, 1)
    h2 = Hyperbola(p2, half, 1)
    c1 = Circle(p1, 1)
    l1 = Line(p1, p2)

    raises(ValueError, lambda: Hyperbola(None, None, None, 1))
    raises(ValueError, lambda: Ellipse())

    # Basic Stuff
    assert h1 != c1
    assert h1 != h2
    assert h1 != l1
    assert p1 not in h1
    assert p3 in h1
    assert h1 in h1
    assert h2 in h2

    # Center
    assert Hyperbola(None, 1, 1).center == Point(0, 0)
    assert h1.center == p1
    assert h2.center == p2

    # Hyperbola with major/minor axis as 0
    assert Hyperbola((1, 1), 0, 0) == Point(1, 1)
    assert Hyperbola((1, 1), 1, 0) == Segment(Point(0, 1), Point(2, 1))
    assert Hyperbola((1, 1), 0, 1) == Segment(Point(1, 0), Point(1, 2))

    # Minor and Major axis
    assert h1.major == 1
    assert h2.major == 1
    assert Hyperbola(p1, a, b).major == a
    assert Hyperbola(p1, a, b).minor == b
    assert Hyperbola(p2, t, t + 1).major == t + 1
    assert Hyperbola(p2, t, t + 1).minor == t

    # Foci and Focal distance
    f1, f2 = Point2D(-2*sqrt(5), 0), Point2D(2*sqrt(5), 0)
    h4 = Hyperbola(Point(0, 0), 4, 2)
    assert h4.focus_distance == 2*sqrt(5)
    assert h4.foci in [(f1, f2), (f2, f1)]

    # Checking Formulae for Properties
    major = 3
    minor = 1
    h5 = Hyperbola(p2, minor, major)
    assert h5.focus_distance == sqrt(major**2 + minor**2)
    ecc = h5.focus_distance / major
    assert h5.eccentricity == ecc
    assert h5.semilatus_rectum == major*(-1 + ecc ** 2)
    # independent of orientation
    h5 = Hyperbola(p2, major, minor)
    assert h5.focus_distance == sqrt(major**2 + minor**2)
    ecc = h5.focus_distance / major
    assert h5.eccentricity == ecc
    assert h5.semilatus_rectum == major*(-1 + ecc ** 2)

    # Encloses
    assert h1.encloses(Segment(Point(-0.5, -0.5), Point(0.5, 0.5))) is True
    assert h1.encloses(Line(p1, p2)) is False
    assert h1.encloses(Ray(p1, p2)) is False
    assert h1.encloses(
        Polygon(Point(-0.5, -0.5), Point(-0.5, 0.5), Point(0.5, 0.5))) is True
    assert h1.encloses(RegularPolygon(p1, 0.5, 3)) is True
    assert h1.encloses(RegularPolygon(p1, 5, 3)) is False
    assert h1.encloses(RegularPolygon(p2, 5, 3)) is False

    assert h2.arbitrary_point() in h2
    raises(ValueError, lambda: Hyperbola(Point(x, y), 1, 1).arbitrary_point(parameter='x'))

    # Tangents
    p1_2 = p2 + Point(half, 0)
    p1_3 = p2 + Point(0, 1)
    assert h2.is_tangent(Line(p1_2, p2 + Point(half, 1)))
    assert h2.is_tangent(Line(p1_3, p2 + Point(half, 1))) is False
    # Checking if asymptote is tangent
    assert h1.is_tangent(Line(Point(0, 0), Point(1, 1))) is False

    # Intersection
    l1 = Line(Point(1, -5), Point(1, 5))
    l4 = Line(Point(-10, 0), Point(0, 10))

    assert intersection(h2, l4) == [Point2D(13/3 - sqrt(403)/3, 43/3 - sqrt(403)/3), Point2D(13/3 + sqrt(403)/3, sqrt(403)/3 + 43/3)]
    assert h1.intersection(l1) == [Point(1, 0)]
    assert h2.intersection(l4) == [Point2D(13/3 - sqrt(403)/3, 43/3 - sqrt(403)/3), Point2D(13/3 + sqrt(403)/3, sqrt(403)/3 + 43/3)]
    assert h1.intersection(Circle(Point(0, 2), 1)) == []
    assert h1.intersection(Circle(Point(5, 0), 1)) == []
    assert h1.intersection(Ellipse(Point(2, 0), 1, 1)) == [Point(1, 0)]
    assert h1.intersection(Ellipse(Point(5, 0), 1, 1)) == []
    assert h1.intersection(Point(2, 0)) == []
    assert h1.intersection(h1) == h1
    assert intersection(Hyperbola(Point(0, 0), 2, 1), Hyperbola(Point(3, 0), 1, 2)) == [Point2D(2, 0), Point2D(22/5, -4*sqrt(6)/5), Point2D(22/5, 4*sqrt(6)/5)]
    assert intersection(Hyperbola(Point(0, 0), 5, 17), Ellipse(Point(4, 0), 1, 0.2)) == [Point(5, 0)]

    raises(TypeError, lambda: intersection(h2, Line((0, 0, 0), (0, 0, 1))))
    raises(TypeError, lambda: intersection(h2, Rational(12)))
    raises(TypeError, lambda: Hyperbola.intersection(h2, 1))

    # encloses_point
    h = Hyperbola((0, 0), 1, 2)
    assert h.encloses_point(h.center)
    assert h.encloses_point(h.center + Point(0, h.vradius - Rational(1, 10)))
    assert h.encloses_point(h.center + Point(h.hradius - Rational(1, 10), 0))
    assert h.encloses_point(h.center + Point(h.hradius, 0)) is False
    assert h.encloses_point(
        h.center + Point(h.hradius + Rational(1, 10), 0)) is False
    h = Hyperbola((0, 0), 2, 1)
    assert h.encloses_point(h.center)
    assert h.encloses_point(h.center + Point(0, h.vradius - Rational(1, 10)))
    assert h.encloses_point(h.center + Point(h.hradius - Rational(1, 10), 0))
    assert h.encloses_point(h.center + Point(h.hradius, 0)) is False
    assert h.encloses_point(
        h.center + Point(h.hradius + Rational(1, 10), 0)) is False


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
