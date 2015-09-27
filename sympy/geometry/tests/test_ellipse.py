from __future__ import division

from sympy import Dummy, Rational, S, Symbol, pi, sqrt
from sympy.core.compatibility import range
from sympy.geometry import (Circle, Ellipse, GeometryError, Line, Point, Polygon, Ray, RegularPolygon, Segment,
                            Triangle, intersection)
from sympy.integrals.integrals import Integral
from sympy.utilities.pytest import raises, slow

x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
t = Symbol('t', real=True)
y1 = Symbol('y1', real=True)
half = Rational(1, 2)


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


def test_repr():
    assert repr(Circle((0, 1), 2)) == 'Circle(Point2D(0, 1), 2)'
