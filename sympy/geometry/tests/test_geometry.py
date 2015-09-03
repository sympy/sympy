from __future__ import division
import warnings

from sympy import (Abs, I, Dummy, Rational, Float, S, Symbol, cos, oo, pi,
                   simplify, sin, sqrt, symbols, Derivative, asin, acos)
from sympy.core.compatibility import range
from sympy.functions.elementary.trigonometric import tan
from sympy.geometry import (Circle, Curve, Ellipse, GeometryError, Line, Point,
                            Polygon, Ray, RegularPolygon, Segment, Triangle,
                            are_similar, convex_hull, intersection,
                            Point3D, Line3D, Ray3D, Segment3D, Plane, centroid)
from sympy.geometry.line import Undecidable
from sympy.geometry.entity import rotate, scale, translate
from sympy.geometry.polygon import _asa as asa, rad, deg
from sympy.geometry.util import idiff, are_coplanar
from sympy.integrals.integrals import Integral
from sympy.matrices import Matrix
from sympy.solvers.solvers import solve
from sympy.utilities.iterables import cartes
from sympy.utilities.randtest import verify_numerically
from sympy.utilities.pytest import raises, slow

x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
t = Symbol('t', real=True)
k = Symbol('k', real=True)
x1 = Symbol('x1', real=True)
x2 = Symbol('x2', real=True)
x3 = Symbol('x3', real=True)
y1 = Symbol('y1', real=True)
y2 = Symbol('y2', real=True)
y3 = Symbol('y3', real=True)
z1 = Symbol('z1', real=True)
z2 = Symbol('z2', real=True)
z3 = Symbol('z3', real=True)
half = Rational(1, 2)


def feq(a, b):
    """Test if two floating point values are 'equal'."""
    t = Float("1.0E-10")
    return -t < a - b < t


def test_curve():
    s = Symbol('s')
    z = Symbol('z')

    # this curve is independent of the indicated parameter
    c = Curve([2*s, s**2], (z, 0, 2))

    assert c.parameter == z
    assert c.functions == (2*s, s**2)
    assert c.arbitrary_point() == Point(2*s, s**2)
    assert c.arbitrary_point(z) == Point(2*s, s**2)

    # this is how it is normally used
    c = Curve([2*s, s**2], (s, 0, 2))

    assert c.parameter == s
    assert c.functions == (2*s, s**2)
    t = Symbol('t')
    # the t returned as assumptions
    assert c.arbitrary_point() != Point(2*t, t**2)
    t = Symbol('t', real=True)
    # now t has the same assumptions so the test passes
    assert c.arbitrary_point() == Point(2*t, t**2)
    assert c.arbitrary_point(z) == Point(2*z, z**2)
    assert c.arbitrary_point(c.parameter) == Point(2*s, s**2)
    assert c.arbitrary_point(None) == Point(2*s, s**2)
    assert c.plot_interval() == [t, 0, 2]
    assert c.plot_interval(z) == [z, 0, 2]

    assert Curve([x, x], (x, 0, 1)).rotate(pi/2, (1, 2)).scale(2, 3).translate(
        1, 3).arbitrary_point(s) == \
        Line((0, 0), (1, 1)).rotate(pi/2, (1, 2)).scale(2, 3).translate(
            1, 3).arbitrary_point(s) == \
        Point(-2*s + 7, 3*s + 6)

    raises(ValueError, lambda: Curve((s), (s, 1, 2)))
    raises(ValueError, lambda: Curve((x, x * 2), (1, x)))

    raises(ValueError, lambda: Curve((s, s + t), (s, 1, 2)).arbitrary_point())
    raises(ValueError, lambda: Curve((s, s + t), (t, 1, 2)).arbitrary_point(s))


@slow
def test_ellipse_geom():
    p1 = Point(0, 0)
    p2 = Point(1, 1)
    p4 = Point(0, 1)

    e1 = Ellipse(p1, 1, 1)
    e2 = Ellipse(p2, half, 1)
    e3 = Ellipse(p1, y1, y1)
    c1 = Circle(p1, 1)
    c2 = Circle(p2, 1)
    c3 = Circle(Point(sqrt(2), sqrt(2)), 1)

    # Test creation with three points
    cen, rad = Point(3*half, 2), 5*half
    assert Circle(Point(0, 0), Point(3, 0), Point(0, 4)) == Circle(cen, rad)
    raises(
        GeometryError, lambda: Circle(Point(0, 0), Point(1, 1), Point(2, 2)))

    raises(ValueError, lambda: Ellipse(None, None, None, 1))
    raises(GeometryError, lambda: Circle(Point(0, 0)))

    # Basic Stuff
    assert Ellipse(None, 1, 1).center == Point(0, 0)
    assert e1 == c1
    assert e1 != e2
    assert p4 in e1
    assert p2 not in e2
    assert e1.area == pi
    assert e2.area == pi/2
    assert e3.area == pi*y1*abs(y1)
    assert c1.area == e1.area
    assert c1.circumference == e1.circumference
    assert e3.circumference == 2*pi*y1
    assert e1.plot_interval() == e2.plot_interval() == [t, -pi, pi]
    assert e1.plot_interval(x) == e2.plot_interval(x) == [x, -pi, pi]
    assert Ellipse(None, 1, None, 1).circumference == 2*pi
    assert c1.minor == 1
    assert c1.major == 1
    assert c1.hradius == 1
    assert c1.vradius == 1

    # Private Functions
    assert hash(c1) == hash(Circle(Point(1, 0), Point(0, 1), Point(0, -1)))
    assert c1 in e1
    assert (Line(p1, p2) in e1) is False
    assert e1.__cmp__(e1) == 0
    assert e1.__cmp__(Point(0, 0)) > 0

    # Encloses
    assert e1.encloses(Segment(Point(-0.5, -0.5), Point(0.5, 0.5))) is True
    assert e1.encloses(Line(p1, p2)) is False
    assert e1.encloses(Ray(p1, p2)) is False
    assert e1.encloses(e1) is False
    assert e1.encloses(
        Polygon(Point(-0.5, -0.5), Point(-0.5, 0.5), Point(0.5, 0.5))) is True
    assert e1.encloses(RegularPolygon(p1, 0.5, 3)) is True
    assert e1.encloses(RegularPolygon(p1, 5, 3)) is False
    assert e1.encloses(RegularPolygon(p2, 5, 3)) is False

    # with generic symbols, the hradius is assumed to contain the major radius
    M = Symbol('M')
    m = Symbol('m')
    c = Ellipse(p1, M, m).circumference
    _x = c.atoms(Dummy).pop()
    assert c == 4*M*Integral(
        sqrt((1 - _x**2*(M**2 - m**2)/M**2)/(1 - _x**2)), (_x, 0, 1))

    assert e2.arbitrary_point() in e2

    # Foci
    f1, f2 = Point(sqrt(12), 0), Point(-sqrt(12), 0)
    ef = Ellipse(Point(0, 0), 4, 2)
    assert ef.foci in [(f1, f2), (f2, f1)]

    # Tangents
    v = sqrt(2) / 2
    p1_1 = Point(v, v)
    p1_2 = p2 + Point(half, 0)
    p1_3 = p2 + Point(0, 1)
    assert e1.tangent_lines(p4) == c1.tangent_lines(p4)
    assert e2.tangent_lines(p1_2) == [Line(Point(3/2, 1), Point(3/2, 1/2))]
    assert e2.tangent_lines(p1_3) == [Line(Point(1, 2), Point(5/4, 2))]
    assert c1.tangent_lines(p1_1) != [Line(p1_1, Point(0, sqrt(2)))]
    assert c1.tangent_lines(p1) == []
    assert e2.is_tangent(Line(p1_2, p2 + Point(half, 1)))
    assert e2.is_tangent(Line(p1_3, p2 + Point(half, 1)))
    assert c1.is_tangent(Line(p1_1, Point(0, sqrt(2))))
    assert e1.is_tangent(Line(Point(0, 0), Point(1, 1))) is False
    assert c1.is_tangent(e1) is False
    assert c1.is_tangent(Ellipse(Point(2, 0), 1, 1)) is True
    assert c1.is_tangent(
        Polygon(Point(1, 1), Point(1, -1), Point(2, 0))) is True
    assert c1.is_tangent(
        Polygon(Point(1, 1), Point(1, 0), Point(2, 0))) is False
    assert Circle(Point(5, 5), 3).is_tangent(Circle(Point(0, 5), 1)) is False

    assert Ellipse(Point(5, 5), 2, 1).tangent_lines(Point(0, 0)) == \
        [Line(Point(0, 0), Point(77/25, 132/25)),
     Line(Point(0, 0), Point(33/5, 22/5))]
    assert Ellipse(Point(5, 5), 2, 1).tangent_lines(Point(3, 4)) == \
        [Line(Point(3, 4), Point(4, 4)), Line(Point(3, 4), Point(3, 5))]
    assert Circle(Point(5, 5), 2).tangent_lines(Point(3, 3)) == \
        [Line(Point(3, 3), Point(4, 3)), Line(Point(3, 3), Point(3, 4))]
    assert Circle(Point(5, 5), 2).tangent_lines(Point(5 - 2*sqrt(2), 5)) == \
        [Line(Point(5 - 2*sqrt(2), 5), Point(5 - sqrt(2), 5 - sqrt(2))),
     Line(Point(5 - 2*sqrt(2), 5), Point(5 - sqrt(2), 5 + sqrt(2))), ]

    e = Ellipse(Point(0, 0), 2, 1)
    assert e.normal_lines(Point(0, 0)) == \
        [Line(Point(0, 0), Point(0, 1)), Line(Point(0, 0), Point(1, 0))]
    assert e.normal_lines(Point(1, 0)) == \
        [Line(Point(0, 0), Point(1, 0))]
    assert e.normal_lines((0, 1)) == \
        [Line(Point(0, 0), Point(0, 1))]
    assert e.normal_lines(Point(1, 1), 2) == [
        Line(Point(-51/26, -1/5), Point(-25/26, 17/83)),
        Line(Point(28/29, -7/8), Point(57/29, -9/2))]
    # test the failure of Poly.intervals and checks a point on the boundary
    p = Point(sqrt(3), S.Half)
    assert p in e
    assert e.normal_lines(p, 2) == [
        Line(Point(-341/171, -1/13), Point(-170/171, 5/64)),
        Line(Point(26/15, -1/2), Point(41/15, -43/26))]
    # be sure to use the slope that isn't undefined on boundary
    e = Ellipse((0, 0), 2, 2*sqrt(3)/3)
    assert e.normal_lines((1, 1), 2) == [
        Line(Point(-64/33, -20/71), Point(-31/33, 2/13)),
        Line(Point(1, -1), Point(2, -4))]
    # general ellipse fails except under certain conditions
    e = Ellipse((0, 0), x, 1)
    assert e.normal_lines((x + 1, 0)) == [Line(Point(0, 0), Point(1, 0))]
    raises(NotImplementedError, lambda: e.normal_lines((x + 1, 1)))


    # Properties
    major = 3
    minor = 1
    e4 = Ellipse(p2, minor, major)
    assert e4.focus_distance == sqrt(major**2 - minor**2)
    ecc = e4.focus_distance / major
    assert e4.eccentricity == ecc
    assert e4.periapsis == major*(1 - ecc)
    assert e4.apoapsis == major*(1 + ecc)
    # independent of orientation
    e4 = Ellipse(p2, major, minor)
    assert e4.focus_distance == sqrt(major**2 - minor**2)
    ecc = e4.focus_distance / major
    assert e4.eccentricity == ecc
    assert e4.periapsis == major*(1 - ecc)
    assert e4.apoapsis == major*(1 + ecc)

    # Intersection
    l1 = Line(Point(1, -5), Point(1, 5))
    l2 = Line(Point(-5, -1), Point(5, -1))
    l3 = Line(Point(-1, -1), Point(1, 1))
    l4 = Line(Point(-10, 0), Point(0, 10))
    pts_c1_l3 = [Point(sqrt(2)/2, sqrt(2)/2), Point(-sqrt(2)/2, -sqrt(2)/2)]

    assert intersection(e2, l4) == []
    assert intersection(c1, Point(1, 0)) == [Point(1, 0)]
    assert intersection(c1, l1) == [Point(1, 0)]
    assert intersection(c1, l2) == [Point(0, -1)]
    assert intersection(c1, l3) in [pts_c1_l3, [pts_c1_l3[1], pts_c1_l3[0]]]
    assert intersection(c1, c2) == [Point(0, 1), Point(1, 0)]
    assert intersection(c1, c3) == [Point(sqrt(2)/2, sqrt(2)/2)]
    assert e1.intersection(l1) == [Point(1, 0)]
    assert e2.intersection(l4) == []
    assert e1.intersection(Circle(Point(0, 2), 1)) == [Point(0, 1)]
    assert e1.intersection(Circle(Point(5, 0), 1)) == []
    assert e1.intersection(Ellipse(Point(2, 0), 1, 1)) == [Point(1, 0)]
    assert e1.intersection(Ellipse(Point(5, 0), 1, 1,)) == []
    assert e1.intersection(Point(2, 0)) == []
    assert e1.intersection(e1) == e1

    # some special case intersections
    csmall = Circle(p1, 3)
    cbig = Circle(p1, 5)
    cout = Circle(Point(5, 5), 1)
    # one circle inside of another
    assert csmall.intersection(cbig) == []
    # separate circles
    assert csmall.intersection(cout) == []
    # coincident circles
    assert csmall.intersection(csmall) == csmall

    v = sqrt(2)
    t1 = Triangle(Point(0, v), Point(0, -v), Point(v, 0))
    points = intersection(t1, c1)
    assert len(points) == 4
    assert Point(0, 1) in points
    assert Point(0, -1) in points
    assert Point(v/2, v/2) in points
    assert Point(v/2, -v/2) in points

    circ = Circle(Point(0, 0), 5)
    elip = Ellipse(Point(0, 0), 5, 20)
    assert intersection(circ, elip) in \
        [[Point(5, 0), Point(-5, 0)], [Point(-5, 0), Point(5, 0)]]
    assert elip.tangent_lines(Point(0, 0)) == []
    elip = Ellipse(Point(0, 0), 3, 2)
    assert elip.tangent_lines(Point(3, 0)) == \
        [Line(Point(3, 0), Point(3, -12))]

    e1 = Ellipse(Point(0, 0), 5, 10)
    e2 = Ellipse(Point(2, 1), 4, 8)
    a = 53/17
    c = 2*sqrt(3991)/17
    ans = [Point(a - c/8, a/2 + c), Point(a + c/8, a/2 - c)]
    assert e1.intersection(e2) == ans
    e2 = Ellipse(Point(x, y), 4, 8)
    c = sqrt(3991)
    ans = [Point(-c/68 + a, 2*c/17 + a/2), Point(c/68 + a, -2*c/17 + a/2)]
    assert [p.subs({x: 2, y:1}) for p in e1.intersection(e2)] == ans

    # Combinations of above
    assert e3.is_tangent(e3.tangent_lines(p1 + Point(y1, 0))[0])

    e = Ellipse((1, 2), 3, 2)
    assert e.tangent_lines(Point(10, 0)) == \
        [Line(Point(10, 0), Point(1, 0)),
        Line(Point(10, 0), Point(14/5, 18/5))]

    # encloses_point
    e = Ellipse((0, 0), 1, 2)
    assert e.encloses_point(e.center)
    assert e.encloses_point(e.center + Point(0, e.vradius - Rational(1, 10)))
    assert e.encloses_point(e.center + Point(e.hradius - Rational(1, 10), 0))
    assert e.encloses_point(e.center + Point(e.hradius, 0)) is False
    assert e.encloses_point(
        e.center + Point(e.hradius + Rational(1, 10), 0)) is False
    e = Ellipse((0, 0), 2, 1)
    assert e.encloses_point(e.center)
    assert e.encloses_point(e.center + Point(0, e.vradius - Rational(1, 10)))
    assert e.encloses_point(e.center + Point(e.hradius - Rational(1, 10), 0))
    assert e.encloses_point(e.center + Point(e.hradius, 0)) is False
    assert e.encloses_point(
        e.center + Point(e.hradius + Rational(1, 10), 0)) is False
    assert c1.encloses_point(Point(1, 0)) is False
    assert c1.encloses_point(Point(0.3, 0.4)) is True

    assert e.scale(2, 3) == Ellipse((0, 0), 4, 3)
    assert e.scale(3, 6) == Ellipse((0, 0), 6, 6)
    assert e.rotate(pi) == e
    assert e.rotate(pi, (1, 2)) == Ellipse(Point(2, 4), 2, 1)
    raises(NotImplementedError, lambda: e.rotate(pi/3))

    # transformations
    c = Circle((1, 1), 2)
    assert c.scale(-1) == Circle((-1, 1), 2)
    assert c.scale(y=-1) == Circle((1, -1), 2)
    assert c.scale(2) == Ellipse((2, 1), 4, 2)


def test_ellipse_random_point():
    e3 = Ellipse(Point(0, 0), y1, y1)
    rx, ry = Symbol('rx'), Symbol('ry')
    for ind in range(0, 5):
        r = e3.random_point()
        # substitution should give zero*y1**2
        assert e3.equation(rx, ry).subs(zip((rx, ry), r.args)).equals(0)


def test_polygon():
    a, b, c = Point(0, 0), Point(2, 0), Point(3, 3)
    t = Triangle(a, b, c)
    assert Polygon(a, Point(1, 0), b, c) == t
    assert Polygon(Point(1, 0), b, c, a) == t
    assert Polygon(b, c, a, Point(1, 0)) == t
    # 2 "remove folded" tests
    assert Polygon(a, Point(3, 0), b, c) == t
    assert Polygon(a, b, Point(3, -1), b, c) == t
    raises(GeometryError, lambda: Polygon((0, 0), (1, 0), (0, 1), (1, 1)))
    # remove multiple collinear points
    assert Polygon(Point(-4, 15), Point(-11, 15), Point(-15, 15),
        Point(-15, 33/5), Point(-15, -87/10), Point(-15, -15),
        Point(-42/5, -15), Point(-2, -15), Point(7, -15), Point(15, -15),
        Point(15, -3), Point(15, 10), Point(15, 15)) == \
        Polygon(Point(-15,-15), Point(15,-15), Point(15,15), Point(-15,15))


    p1 = Polygon(
        Point(0, 0), Point(3, -1),
        Point(6, 0), Point(4, 5),
        Point(2, 3), Point(0, 3))
    p2 = Polygon(
        Point(6, 0), Point(3, -1),
        Point(0, 0), Point(0, 3),
        Point(2, 3), Point(4, 5))
    p3 = Polygon(
        Point(0, 0), Point(3, 0),
        Point(5, 2), Point(4, 4))
    p4 = Polygon(
        Point(0, 0), Point(4, 4),
        Point(5, 2), Point(3, 0))
    p5 = Polygon(
        Point(0, 0), Point(4, 4),
        Point(0, 4))
    p6 = Polygon(
        Point(-11, 1), Point(-9, 6.6),
        Point(-4, -3), Point(-8.4, -8.7))
    r = Ray(Point(-9,6.6), Point(-9,5.5))
    #
    # General polygon
    #
    assert p1 == p2
    assert len(p1.args) == 6
    assert len(p1.sides) == 6
    assert p1.perimeter == 5 + 2*sqrt(10) + sqrt(29) + sqrt(8)
    assert p1.area == 22
    assert not p1.is_convex()
    # ensure convex for both CW and CCW point specification
    assert p3.is_convex()
    assert p4.is_convex()
    dict5 = p5.angles
    assert dict5[Point(0, 0)] == pi / 4
    assert dict5[Point(0, 4)] == pi / 2
    assert p5.encloses_point(Point(x, y)) is None
    assert p5.encloses_point(Point(1, 3))
    assert p5.encloses_point(Point(0, 0)) is False
    assert p5.encloses_point(Point(4, 0)) is False
    p5.plot_interval('x') == [x, 0, 1]
    assert p5.distance(
        Polygon(Point(10, 10), Point(14, 14), Point(10, 14))) == 6 * sqrt(2)
    assert p5.distance(
        Polygon(Point(1, 8), Point(5, 8), Point(8, 12), Point(1, 12))) == 4
    warnings.filterwarnings(
        "error", message="Polygons may intersect producing erroneous output")
    raises(UserWarning,
           lambda: Polygon(Point(0, 0), Point(1, 0),
           Point(1, 1)).distance(
           Polygon(Point(0, 0), Point(0, 1), Point(1, 1))))
    warnings.filterwarnings(
        "ignore", message="Polygons may intersect producing erroneous output")
    assert hash(p5) == hash(Polygon(Point(0, 0), Point(4, 4), Point(0, 4)))
    assert p5 == Polygon(Point(4, 4), Point(0, 4), Point(0, 0))
    assert Polygon(Point(4, 4), Point(0, 4), Point(0, 0)) in p5
    assert p5 != Point(0, 4)
    assert Point(0, 1) in p5
    assert p5.arbitrary_point('t').subs(Symbol('t', real=True), 0) == \
        Point(0, 0)
    raises(ValueError, lambda: Polygon(
        Point(x, 0), Point(0, y), Point(x, y)).arbitrary_point('x'))
    assert p6.intersection(r) == [Point(-9, 33/5), Point(-9, -84/13)]
    #
    # Regular polygon
    #
    p1 = RegularPolygon(Point(0, 0), 10, 5)
    p2 = RegularPolygon(Point(0, 0), 5, 5)
    raises(GeometryError, lambda: RegularPolygon(Point(0, 0), Point(0,
           1), Point(1, 1)))
    raises(GeometryError, lambda: RegularPolygon(Point(0, 0), 1, 2))
    raises(ValueError, lambda: RegularPolygon(Point(0, 0), 1, 2.5))

    assert p1 != p2
    assert p1.interior_angle == 3*pi/5
    assert p1.exterior_angle == 2*pi/5
    assert p2.apothem == 5*cos(pi/5)
    assert p2.circumcenter == p1.circumcenter == Point(0, 0)
    assert p1.circumradius == p1.radius == 10
    assert p2.circumcircle == Circle(Point(0, 0), 5)
    assert p2.incircle == Circle(Point(0, 0), p2.apothem)
    assert p2.inradius == p2.apothem == (5 * (1 + sqrt(5)) / 4)
    p2.spin(pi / 10)
    dict1 = p2.angles
    assert dict1[Point(0, 5)] == 3 * pi / 5
    assert p1.is_convex()
    assert p1.rotation == 0
    assert p1.encloses_point(Point(0, 0))
    assert p1.encloses_point(Point(11, 0)) is False
    assert p2.encloses_point(Point(0, 4.9))
    p1.spin(pi/3)
    assert p1.rotation == pi/3
    assert p1.vertices[0] == Point(5, 5*sqrt(3))
    for var in p1.args:
        if isinstance(var, Point):
            assert var == Point(0, 0)
        else:
            assert var == 5 or var == 10 or var == pi / 3
    assert p1 != Point(0, 0)
    assert p1 != p5

    # while spin works in place (notice that rotation is 2pi/3 below)
    # rotate returns a new object
    p1_old = p1
    assert p1.rotate(pi/3) == RegularPolygon(Point(0, 0), 10, 5, 2*pi/3)
    assert p1 == p1_old

    assert p1.area == (-250*sqrt(5) + 1250)/(4*tan(pi/5))
    assert p1.length == 20*sqrt(-sqrt(5)/8 + 5/8)
    assert p1.scale(2, 2) == \
        RegularPolygon(p1.center, p1.radius*2, p1._n, p1.rotation)
    assert RegularPolygon((0, 0), 1, 4).scale(2, 3) == \
        Polygon(Point(2, 0), Point(0, 3), Point(-2, 0), Point(0, -3))

    assert repr(p1) == str(p1)

    #
    # Angles
    #
    angles = p4.angles
    assert feq(angles[Point(0, 0)].evalf(), Float("0.7853981633974483"))
    assert feq(angles[Point(4, 4)].evalf(), Float("1.2490457723982544"))
    assert feq(angles[Point(5, 2)].evalf(), Float("1.8925468811915388"))
    assert feq(angles[Point(3, 0)].evalf(), Float("2.3561944901923449"))

    angles = p3.angles
    assert feq(angles[Point(0, 0)].evalf(), Float("0.7853981633974483"))
    assert feq(angles[Point(4, 4)].evalf(), Float("1.2490457723982544"))
    assert feq(angles[Point(5, 2)].evalf(), Float("1.8925468811915388"))
    assert feq(angles[Point(3, 0)].evalf(), Float("2.3561944901923449"))

    #
    # Triangle
    #
    p1 = Point(0, 0)
    p2 = Point(5, 0)
    p3 = Point(0, 5)
    t1 = Triangle(p1, p2, p3)
    t2 = Triangle(p1, p2, Point(Rational(5, 2), sqrt(Rational(75, 4))))
    t3 = Triangle(p1, Point(x1, 0), Point(0, x1))
    s1 = t1.sides
    assert Triangle(p1, p2, p1) == Polygon(p1, p2, p1) == Segment(p1, p2)
    raises(GeometryError, lambda: Triangle(Point(0, 0)))

    # Basic stuff
    assert Triangle(p1, p1, p1) == p1
    assert Triangle(p2, p2*2, p2*3) == Segment(p2, p2*3)
    assert t1.area == Rational(25, 2)
    assert t1.is_right()
    assert t2.is_right() is False
    assert t3.is_right()
    assert p1 in t1
    assert t1.sides[0] in t1
    assert Segment((0, 0), (1, 0)) in t1
    assert Point(5, 5) not in t2
    assert t1.is_convex()
    assert feq(t1.angles[p1].evalf(), pi.evalf()/2)

    assert t1.is_equilateral() is False
    assert t2.is_equilateral()
    assert t3.is_equilateral() is False
    assert are_similar(t1, t2) is False
    assert are_similar(t1, t3)
    assert are_similar(t2, t3) is False
    assert t1.is_similar(Point(0, 0)) is False

    # Bisectors
    bisectors = t1.bisectors()
    assert bisectors[p1] == Segment(p1, Point(Rational(5, 2), Rational(5, 2)))
    ic = (250 - 125*sqrt(2)) / 50
    assert t1.incenter == Point(ic, ic)

    # Inradius
    assert t1.inradius == t1.incircle.radius == 5 - 5*sqrt(2)/2
    assert t2.inradius == t2.incircle.radius == 5*sqrt(3)/6
    assert t3.inradius == t3.incircle.radius == x1**2/((2 + sqrt(2))*Abs(x1))

    # Circumcircle
    assert t1.circumcircle.center == Point(2.5, 2.5)

    # Medians + Centroid
    m = t1.medians
    assert t1.centroid == Point(Rational(5, 3), Rational(5, 3))
    assert m[p1] == Segment(p1, Point(Rational(5, 2), Rational(5, 2)))
    assert t3.medians[p1] == Segment(p1, Point(x1/2, x1/2))
    assert intersection(m[p1], m[p2], m[p3]) == [t1.centroid]
    assert t1.medial == Triangle(Point(2.5, 0), Point(0, 2.5), Point(2.5, 2.5))

    # Perpendicular
    altitudes = t1.altitudes
    assert altitudes[p1] == Segment(p1, Point(Rational(5, 2), Rational(5, 2)))
    assert altitudes[p2] == s1[0]
    assert altitudes[p3] == s1[2]
    assert t1.orthocenter == p1
    t = S('''Triangle(
    Point(100080156402737/5000000000000, 79782624633431/500000000000),
    Point(39223884078253/2000000000000, 156345163124289/1000000000000),
    Point(31241359188437/1250000000000, 338338270939941/1000000000000000))''')
    assert t.orthocenter == S('''Point(-780660869050599840216997'''
    '''79471538701955848721853/80368430960602242240789074233100000000000000,'''
    '''20151573611150265741278060334545897615974257/16073686192120448448157'''
    '''8148466200000000000)''')

    # Ensure
    assert len(intersection(*bisectors.values())) == 1
    assert len(intersection(*altitudes.values())) == 1
    assert len(intersection(*m.values())) == 1

    # Distance
    p1 = Polygon(
        Point(0, 0), Point(1, 0),
        Point(1, 1), Point(0, 1))
    p2 = Polygon(
        Point(0, Rational(5)/4), Point(1, Rational(5)/4),
        Point(1, Rational(9)/4), Point(0, Rational(9)/4))
    p3 = Polygon(
        Point(1, 2), Point(2, 2),
        Point(2, 1))
    p4 = Polygon(
        Point(1, 1), Point(Rational(6)/5, 1),
        Point(1, Rational(6)/5))
    pt1 = Point(half, half)
    pt2 = Point(1, 1)

    '''Polygon to Point'''
    assert p1.distance(pt1) == half
    assert p1.distance(pt2) == 0
    assert p2.distance(pt1) == Rational(3)/4
    assert p3.distance(pt2) == sqrt(2)/2

    '''Polygon to Polygon'''
    # p1.distance(p2) emits a warning
    # First, test the warning
    warnings.filterwarnings("error",
        message="Polygons may intersect producing erroneous output")
    raises(UserWarning, lambda: p1.distance(p2))
    # now test the actual output
    warnings.filterwarnings("ignore",
        message="Polygons may intersect producing erroneous output")
    assert p1.distance(p2) == half/2

    assert p1.distance(p3) == sqrt(2)/2
    assert p3.distance(p4) == (sqrt(2)/2 - sqrt(Rational(2)/25)/2)


def test_convex_hull():
    p = [Point(-5, -1), Point(-2, 1), Point(-2, -1), Point(-1, -3),
        Point(0, 0), Point(1, 1), Point(2, 2), Point(2, -1), Point(3, 1),
        Point(4, -1), Point(6, 2)]
    ch = Polygon(p[0], p[3], p[9], p[10], p[6], p[1])
    #test handling of duplicate points
    p.append(p[3])

    #more than 3 collinear points
    another_p = [Point(-45, -85), Point(-45, 85), Point(-45, 26),
                 Point(-45, -24)]
    ch2 = Segment(another_p[0], another_p[1])

    assert convex_hull(*another_p) == ch2
    assert convex_hull(*p) == ch
    assert convex_hull(p[0]) == p[0]
    assert convex_hull(p[0], p[1]) == Segment(p[0], p[1])

    # no unique points
    assert convex_hull(*[p[-1]]*3) == p[-1]

    # collection of items
    assert convex_hull(*[Point(0, 0),
                        Segment(Point(1, 0), Point(1, 1)),
                        RegularPolygon(Point(2, 0), 2, 4)]) == \
        Polygon(Point(0, 0), Point(2, -2), Point(4, 0), Point(2, 2))


def test_concyclic_doctest_bug():
    p1, p2 = Point(-1, 0), Point(1, 0)
    p3, p4 = Point(0, 1), Point(-1, 2)
    assert Point.is_concyclic(p1, p2, p3)
    assert not Point.is_concyclic(p1, p2, p3, p4)


def test_subs():
    p = Point(x, 2)
    q = Point(1, 1)
    r = Point(3, 4)
    for o in [p,
              Segment(p, q),
              Ray(p, q),
              Line(p, q),
              Triangle(p, q, r),
              RegularPolygon(p, 3, 6),
              Polygon(p, q, r, Point(5, 4)),
              Circle(p, 3),
              Ellipse(p, 3, 4)]:
        assert 'y' in str(o.subs(x, y))
    assert p.subs({x: 1}) == Point(1, 2)
    assert Point(1, 2).subs(Point(1, 2), Point(3, 4)) == Point(3, 4)
    assert Point(1, 2).subs((1, 2), Point(3, 4)) == Point(3, 4)
    assert Point(1, 2).subs(Point(1, 2), Point(3, 4)) == Point(3, 4)
    assert Point(1, 2).subs(set([(1, 2)])) == Point(2, 2)
    raises(ValueError, lambda: Point(1, 2).subs(1))
    raises(ValueError, lambda: Point(1, 1).subs((Point(1, 1), Point(1,
           2)), 1, 2))


def test_encloses():
    # square with a dimpled left side
    s = Polygon(Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1),
        Point(S.Half, S.Half))
    # the following is True if the polygon isn't treated as closing on itself
    assert s.encloses(Point(0, S.Half)) is False
    assert s.encloses(Point(S.Half, S.Half)) is False  # it's a vertex
    assert s.encloses(Point(Rational(3, 4), S.Half)) is True


def test_free_symbols():
    a, b, c, d, e, f, s = symbols('a:f,s')
    assert Point(a, b).free_symbols == set([a, b])
    assert Line((a, b), (c, d)).free_symbols == set([a, b, c, d])
    assert Ray((a, b), (c, d)).free_symbols == set([a, b, c, d])
    assert Ray((a, b), angle=c).free_symbols == set([a, b, c])
    assert Segment((a, b), (c, d)).free_symbols == set([a, b, c, d])
    assert Line((a, b), slope=c).free_symbols == set([a, b, c])
    assert Curve((a*s, b*s), (s, c, d)).free_symbols == set([a, b, c, d])
    assert Ellipse((a, b), c, d).free_symbols == set([a, b, c, d])
    assert Ellipse((a, b), c, eccentricity=d).free_symbols == \
        set([a, b, c, d])
    assert Ellipse((a, b), vradius=c, eccentricity=d).free_symbols == \
        set([a, b, c, d])
    assert Circle((a, b), c).free_symbols == set([a, b, c])
    assert Circle((a, b), (c, d), (e, f)).free_symbols == \
        set([e, d, c, b, f, a])
    assert Polygon((a, b), (c, d), (e, f)).free_symbols == \
        set([e, b, d, f, a, c])
    assert RegularPolygon((a, b), c, d, e).free_symbols == set([e, a, b, c, d])


def test_util_centroid():
    p = Polygon((0, 0), (10, 0), (10, 10))
    q = p.translate(0, 20)
    assert centroid(p, q) == Point(20, 40)/3
    p = Segment((0, 0), (2, 0))
    q = Segment((0, 0), (2, 2))
    assert centroid(p, q) == Point(1, -sqrt(2) + 2)
    assert centroid(Point(0, 0), Point(2, 0)) == Point(2, 0)/2
    assert centroid(Point(0, 0), Point(0, 0), Point(2, 0)) == Point(2, 0)/3


def test_util():
    # coverage for some leftover functions in sympy.geometry.util
    assert intersection(Point(0, 0)) == []
    raises(ValueError, lambda: intersection(Point(0, 0), 3))
    raises(ValueError, lambda: convex_hull(Point(0, 0), 3))


def test_repr():
    assert repr(Circle((0, 1), 2)) == 'Circle(Point2D(0, 1), 2)'


def test_transform():
    p = Point(1, 1)
    assert p.transform(rotate(pi/2)) == Point(-1, 1)
    assert p.transform(scale(3, 2)) == Point(3, 2)
    assert p.transform(translate(1, 2)) == Point(2, 3)


def test_triangle_kwargs():
    assert Triangle(sss=(3, 4, 5)) == \
        Triangle(Point(0, 0), Point(3, 0), Point(3, 4))
    assert Triangle(asa=(30, 2, 30)) == \
        Triangle(Point(0, 0), Point(2, 0), Point(1, sqrt(3)/3))
    assert Triangle(sas=(1, 45, 2)) == \
        Triangle(Point(0, 0), Point(2, 0), Point(sqrt(2)/2, sqrt(2)/2))
    assert Triangle(sss=(1, 2, 5)) is None
    assert deg(rad(180)) == 180


def test_geometry_transforms():
    from sympy import Tuple
    c = Curve((x, x**2), (x, 0, 1))
    pts = [Point(0, 0), Point(1/2, 1/4), Point(1, 1)]
    cout = Curve((2*x - 4, 3*x**2 - 10), (x, 0, 1))
    pts_out = [Point(-4, -10), Point(-3, -37/4), Point(-2, -7)]
    assert c.scale(2, 3, (4, 5)) == cout
    assert [c.subs(x, xi/2) for xi in Tuple(0, 1, 2)] == pts
    assert [cout.subs(x, xi/2) for xi in Tuple(0, 1, 2)] == pts_out
    assert Triangle(*pts).scale(2, 3, (4, 5)) == Triangle(*pts_out)

    assert Ellipse((0, 0), 2, 3).scale(2, 3, (4, 5)) == \
        Ellipse(Point(-4, -10), 4, 9)
    assert Circle((0, 0), 2).scale(2, 3, (4, 5)) == \
        Ellipse(Point(-4, -10), 4, 6)
    assert Ellipse((0, 0), 2, 3).scale(3, 3, (4, 5)) == \
        Ellipse(Point(-8, -10), 6, 9)
    assert Circle((0, 0), 2).scale(3, 3, (4, 5)) == \
        Circle(Point(-8, -10), 6)
    assert Circle(Point(-8, -10), 6).scale(1/3, 1/3, (4, 5)) == \
        Circle((0, 0), 2)
    assert Curve((x + y, 3*x), (x, 0, 1)).subs(y, S.Half) == \
        Curve((x + 1/2, 3*x), (x, 0, 1))
    assert Curve((x, 3*x), (x, 0, 1)).translate(4, 5) == \
        Curve((x + 4, 3*x + 5), (x, 0, 1))
    assert Circle((0, 0), 2).translate(4, 5) == \
        Circle((4, 5), 2)
    assert Circle((0, 0), 2).scale(3, 3) == \
        Circle((0, 0), 6)
    assert Point(1, 1).scale(2, 3, (4, 5)) == \
        Point(-2, -7)
    assert Point(1, 1).translate(4, 5) == \
        Point(5, 6)
    assert scale(1, 2, (3, 4)).tolist() == \
        [[1, 0, 0], [0, 2, 0], [0, -4, 1]]
    assert RegularPolygon((0, 0), 1, 4).scale(2, 3, (4, 5)) == \
        Polygon(Point(-2, -10), Point(-4, -7), Point(-6, -10), Point(-4, -13))


def test_reflect():
    b = Symbol('b')
    m = Symbol('m')
    l = Line((0, b), slope=m)
    p = Point(x, y)
    r = p.reflect(l)
    dp = l.perpendicular_segment(p).length
    dr = l.perpendicular_segment(r).length
    assert verify_numerically(dp, dr)
    t = Triangle((0, 0), (1, 0), (2, 3))
    assert t.area == -t.reflect(l).area
    e = Ellipse((1, 0), 1, 2)
    assert e.area == -e.reflect(Line((1, 0), slope=0)).area
    assert e.area == -e.reflect(Line((1, 0), slope=oo)).area
    raises(NotImplementedError, lambda: e.reflect(Line((1, 0), slope=m)))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((3, 0), slope=oo)) \
        == Triangle(Point(5, 0), Point(4, 0), Point(4, 2))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((0, 3), slope=oo)) \
        == Triangle(Point(-1, 0), Point(-2, 0), Point(-2, 2))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((0, 3), slope=0)) \
        == Triangle(Point(1, 6), Point(2, 6), Point(2, 4))
    assert Polygon((1, 0), (2, 0), (2, 2)).reflect(Line((3, 0), slope=0)) \
        == Triangle(Point(1, 0), Point(2, 0), Point(2, -2))

    # test entity overrides
    c = Circle((x, y), 3)
    cr = c.reflect(l)
    assert cr == Circle(r, -3)
    assert c.area == -cr.area
    pent = RegularPolygon((1, 2), 1, 5)
    l = Line((0, pi), slope=sqrt(2))
    rpent = pent.reflect(l)
    poly_pent = Polygon(*pent.vertices)
    assert rpent.center == pent.center.reflect(l)
    assert str([w.n(3) for w in rpent.vertices]) == (
        '[Point2D(-0.586, 4.27), Point2D(-1.69, 4.66), '
        'Point2D(-2.41, 3.73), Point2D(-1.74, 2.76), '
        'Point2D(-0.616, 3.10)]')
    assert pent.area.equals(-rpent.area)


def test_idiff():
    # the use of idiff in ellipse also provides coverage
    circ = x**2 + y**2 - 4
    ans = -3*x*(x**2 + y**2)/y**5
    assert ans == idiff(circ, y, x, 3).simplify()
    assert ans == idiff(circ, [y], x, 3).simplify()
    assert idiff(circ, y, x, 3).simplify() == ans
    explicit  = 12*x/sqrt(-x**2 + 4)**5
    assert ans.subs(y, solve(circ, y)[0]).equals(explicit)
    assert True in [sol.diff(x, 3).equals(explicit) for sol in solve(circ, y)]
    assert idiff(x + t + y, [y, t], x) == -Derivative(t, x) - 1
