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


@slow
def test_plane():
    p1 = Point3D(0, 0, 0)
    p2 = Point3D(1, 1, 1)
    p3 = Point3D(1, 2, 3)
    p4 = Point3D(x, x, x)
    p5 = Point3D(y, y, y)

    pl3 = Plane(p1, p2, p3)
    pl4 = Plane(p1, normal_vector=(1, 1, 1))
    pl4b = Plane(p1, p2)
    pl5 = Plane(p3, normal_vector=(1, 2, 3))
    pl6 = Plane(Point3D(2, 3, 7), normal_vector=(2, 2, 2))
    pl7 = Plane(Point3D(1, -5, -6), normal_vector=(1, -2, 1))

    l1 = Line3D(Point3D(5, 0, 0), Point3D(1, -1, 1))
    l2 = Line3D(Point3D(0, -2, 0), Point3D(3, 1, 1))
    l3 = Line3D(Point3D(0, -1, 0), Point3D(5, -1, 9))

    assert Plane(p1, p2, p3) != Plane(p1, p3, p2)
    assert Plane(p1, p2, p3).is_coplanar(Plane(p1, p3, p2))
    assert pl3 == Plane(Point3D(0, 0, 0), normal_vector=(1, -2, 1))
    assert pl3 != pl4
    assert pl4 == pl4b
    assert pl5 == Plane(Point3D(1, 2, 3), normal_vector=(1, 2, 3))

    assert pl5.equation(x, y, z) == x + 2*y + 3*z - 14
    assert pl3.equation(x, y, z) == x - 2*y + z

    assert pl3.p1 == p1
    assert pl4.p1 == p1
    assert pl5.p1 == p3

    assert pl4.normal_vector == (1, 1, 1)
    assert pl5.normal_vector == (1, 2, 3)

    assert p1 in pl3
    assert p1 in pl4
    assert p3 in pl5

    assert pl3.projection(Point(0, 0)) == p1
    p = pl3.projection(Point3D(1, 1, 0))
    assert p == Point3D(7/6, 2/3, 1/6)
    assert p in pl3

    l = pl3.projection_line(Line(Point(0, 0), Point(1, 1)))
    assert l == Line3D(Point3D(0, 0, 0), Point3D(7/6, 2/3, 1/6))
    assert l in pl3
    # get a segment that does not intersect the plane which is also
    # parallel to pl3's normal veector
    t = Dummy()
    r = pl3.random_point()
    a = pl3.perpendicular_line(r).arbitrary_point(t)
    s = Segment3D(a.subs(t, 1), a.subs(t, 2))
    assert s.p1 not in pl3 and s.p2 not in pl3
    assert pl3.projection_line(s).equals(r)
    assert pl3.projection_line(Segment(Point(1, 0), Point(1, 1))) == \
               Segment3D(Point3D(5/6, 1/3, -1/6), Point3D(7/6, 2/3, 1/6))
    assert pl6.projection_line(Ray(Point(1, 0), Point(1, 1))) == \
               Ray3D(Point3D(14/3, 11/3, 11/3), Point3D(13/3, 13/3, 10/3))
    assert pl3.perpendicular_line(r.args) == pl3.perpendicular_line(r)


    assert pl3.is_parallel(pl6) is False
    assert pl4.is_parallel(pl6)
    assert pl6.is_parallel(l1) is False

    assert pl3.is_perpendicular(pl6)
    assert pl4.is_perpendicular(pl7)
    assert pl6.is_perpendicular(pl7)
    assert pl6.is_perpendicular(l1) is False

    assert pl7.distance(Point3D(1, 3, 5)) == 5*sqrt(6)/6
    assert pl6.distance(Point3D(0, 0, 0)) == 4*sqrt(3)
    assert pl6.distance(pl6.p1) == 0
    assert pl7.distance(pl6) == 0
    assert pl7.distance(l1) == 0
    assert pl6.distance(Segment3D(Point3D(2, 3, 1), Point3D(1, 3, 4))) == 0
    pl6.distance(Plane(Point3D(5, 5, 5), normal_vector=(8, 8, 8))) == sqrt(3)

    assert pl6.angle_between(pl3) == pi/2
    assert pl6.angle_between(pl6) == 0
    assert pl6.angle_between(pl4) == 0
    assert pl7.angle_between(Line3D(Point3D(2, 3, 5), Point3D(2, 4, 6))) == \
        -asin(sqrt(3)/6)
    assert pl6.angle_between(Ray3D(Point3D(2, 4, 1), Point3D(6, 5, 3))) == \
        asin(sqrt(7)/3)
    assert pl7.angle_between(Segment3D(Point3D(5, 6, 1), Point3D(1, 2, 4))) == \
        -asin(7*sqrt(246)/246)

    assert are_coplanar(l1, l2, l3) is False
    assert are_coplanar(l1) is False
    assert are_coplanar(Point3D(2, 7, 2), Point3D(0, 0, 2),
        Point3D(1, 1, 2), Point3D(1, 2, 2))
    assert are_coplanar(Plane(p1, p2, p3), Plane(p1, p3, p2))
    assert Plane.are_concurrent(pl3, pl4, pl5) is False
    assert Plane.are_concurrent(pl6) is False
    raises(ValueError, lambda: Plane.are_concurrent(Point3D(0, 0, 0)))

    assert pl3.parallel_plane(Point3D(1, 2, 5)) == Plane(Point3D(1, 2, 5), \
                                                      normal_vector=(1, -2, 1))

    # perpendicular_plane
    p = Plane((0, 0, 0), (1, 0, 0))
    # default
    assert p.perpendicular_plane() == Plane(Point3D(0, 0, 0), (0, 1, 0))
    # 1 pt
    assert p.perpendicular_plane(Point3D(1, 0, 1)) == \
        Plane(Point3D(1, 0, 1), (0, 1, 0))
    # pts as tuples
    assert p.perpendicular_plane((1, 0, 1), (1, 1, 1)) == \
        Plane(Point3D(1, 0, 1), (0, 0, -1))

    a, b = Point3D(0, 0, 0), Point3D(0, 1, 0)
    Z = (0, 0, 1)
    p = Plane(a, normal_vector=Z)
    # case 4
    assert p.perpendicular_plane(a, b) == Plane(a, (1, 0, 0))
    n = Point3D(*Z)
    # case 1
    assert p.perpendicular_plane(a, n) == Plane(a, (-1, 0, 0))
    # case 2
    assert Plane(a, normal_vector=b.args).perpendicular_plane(a, a + b) == \
        Plane(Point3D(0, 0, 0), (1, 0, 0))
    # case 1&3
    assert Plane(b, normal_vector=Z).perpendicular_plane(b, b + n) == \
        Plane(Point3D(0, 1, 0), (-1, 0, 0))
    # case 2&3
    assert Plane(b, normal_vector=b.args).perpendicular_plane(n, n + b) == \
        Plane(Point3D(0, 0, 1), (1, 0, 0))

    assert pl6.intersection(pl6) == [pl6]
    assert pl4.intersection(pl4.p1) == [pl4.p1]
    assert pl3.intersection(pl6) == [
        Line3D(Point3D(8, 4, 0), Point3D(2, 4, 6))]
    assert pl3.intersection(Line3D(Point3D(1,2,4), Point3D(4,4,2))) == [
        Point3D(2, 8/3, 10/3)]
    assert pl3.intersection(Plane(Point3D(6, 0, 0), normal_vector=(2, -5, 3))
        ) == [Line3D(Point3D(-24, -12, 0), Point3D(-25, -13, -1))]
    assert pl6.intersection(Ray3D(Point3D(2, 3, 1), Point3D(1, 3, 4))) == [
        Point3D(-1, 3, 10)]
    assert pl6.intersection(Segment3D(Point3D(2, 3, 1), Point3D(1, 3, 4))) == [
        Point3D(-1, 3, 10)]
    assert pl7.intersection(Line(Point(2, 3), Point(4, 2))) == [
        Point3D(13/2, 3/4, 0)]
    r = Ray(Point(2, 3), Point(4, 2))
    assert Plane((1,2,0), normal_vector=(0,0,1)).intersection(r) == [
        Ray3D(Point(2, 3), Point(4, 2))]

    assert pl3.random_point() in pl3

    # issue 8570
    l2 = Line3D(Point3D(S(50000004459633)/5000000000000,
                        -S(891926590718643)/1000000000000000,
                        S(231800966893633)/100000000000000),
                Point3D(S(50000004459633)/50000000000000,
                        -S(222981647679771)/250000000000000,
                        S(231800966893633)/100000000000000))

    p2 = Plane(Point3D(S(402775636372767)/100000000000000,
                       -S(97224357654973)/100000000000000,
                       S(216793600814789)/100000000000000),
               (-S('9.00000087501922'), -S('4.81170658872543e-13'),
                S('0.0')))

    assert str([i.n(2) for i in p2.intersection(l2)]) == \
           '[Point3D(4.0, -0.89, 2.3)]'
