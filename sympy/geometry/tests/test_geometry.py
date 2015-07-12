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


def test_line_geom():
    p1 = Point(0, 0)
    p2 = Point(1, 1)
    p3 = Point(x1, x1)
    p4 = Point(y1, y1)
    p5 = Point(x1, 1 + x1)
    p6 = Point(1, 0)
    p7 = Point(0, 1)
    p8 = Point(2, 0)
    p9 = Point(2, 1)

    l1 = Line(p1, p2)
    l2 = Line(p3, p4)
    l3 = Line(p3, p5)
    l4 = Line(p1, p6)
    l5 = Line(p1, p7)
    l6 = Line(p8, p9)
    l7 = Line(p2, p9)
    raises(ValueError, lambda: Line(Point(0, 0), Point(0, 0)))

    # Basic stuff
    assert Line((1, 1), slope=1) == Line((1, 1), (2, 2))
    assert Line((1, 1), slope=oo) == Line((1, 1), (1, 2))
    assert Line((1, 1), slope=-oo) == Line((1, 1), (1, 2))
    raises(ValueError, lambda: Line((1, 1), 1))
    assert Line(p1, p2) == Line(p1, p2)
    assert Line(p1, p2) != Line(p2, p1)
    assert l1 != l2
    assert l1 != l3
    assert l1.slope == 1
    assert l1.length == oo
    assert l3.slope == oo
    assert l4.slope == 0
    assert l4.coefficients == (0, 1, 0)
    assert l4.equation(x=x, y=y) == y
    assert l5.slope == oo
    assert l5.coefficients == (1, 0, 0)
    assert l5.equation() == x
    assert l6.equation() == x - 2
    assert l7.equation() == y - 1
    assert p1 in l1  # is p1 on the line l1?
    assert p1 not in l3
    assert Line((-x, x), (-x + 1, x - 1)).coefficients == (1, 1, 0)

    assert simplify(l1.equation()) in (x - y, y - x)
    assert simplify(l3.equation()) in (x - x1, x1 - x)

    assert Line(p1, p2).scale(2, 1) == Line(p1, p9)

    assert l2.arbitrary_point() in l2
    for ind in range(0, 5):
        assert l3.random_point() in l3

    # Orthogonality
    p1_1 = Point(-x1, x1)
    l1_1 = Line(p1, p1_1)
    assert l1.perpendicular_line(p1.args) == Line(Point(0, 0), Point(1, -1))
    assert l1.perpendicular_line(p1) == Line(Point(0, 0), Point(1, -1))
    assert Line.is_perpendicular(l1, l1_1)
    assert Line.is_perpendicular(l1, l2) is False
    p = l1.random_point()
    assert l1.perpendicular_segment(p) == p

    # Parallelity
    l2_1 = Line(p3, p5)
    assert l2.parallel_line(p1_1) == Line(Point(-x1, x1), Point(-y1, 2*x1 - y1))
    assert l2_1.parallel_line(p1.args) == Line(Point(0, 0), Point(0, -1))
    assert l2_1.parallel_line(p1) == Line(Point(0, 0), Point(0, -1))
    assert Line.is_parallel(l1, l2)
    assert Line.is_parallel(l2, l3) is False
    assert Line.is_parallel(l2, l2.parallel_line(p1_1))
    assert Line.is_parallel(l2_1, l2_1.parallel_line(p1))

    # Intersection
    assert intersection(l1, p1) == [p1]
    assert intersection(l1, p5) == []
    assert intersection(l1, l2) in [[l1], [l2]]
    assert intersection(l1, l1.parallel_line(p5)) == []

    # Concurrency
    l3_1 = Line(Point(5, x1), Point(-Rational(3, 5), x1))
    assert Line.are_concurrent(l1) is False
    assert Line.are_concurrent(l1, l3)
    assert Line.are_concurrent(l1, l3, l3_1)
    assert Line.are_concurrent(l1, l1_1, l3) is False

    # Projection
    assert l2.projection(p4) == p4
    assert l1.projection(p1_1) == p1
    assert l3.projection(p2) == Point(x1, 1)
    raises(GeometryError, lambda: Line(Point(0, 0), Point(1, 0))
           .projection(Circle(Point(0, 0), 1)))

    # Finding angles
    l1_1 = Line(p1, Point(5, 0))
    assert feq(Line.angle_between(l1, l1_1).evalf(), pi.evalf()/4)

    # Testing Rays and Segments (very similar to Lines)
    assert Ray((1, 1), angle=pi/4) == Ray((1, 1), (2, 2))
    assert Ray((1, 1), angle=pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=-pi/2) == Ray((1, 1), (1, 0))
    assert Ray((1, 1), angle=-3*pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=5*pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=5.0*pi/2) == Ray((1, 1), (1, 2))
    assert Ray((1, 1), angle=pi) == Ray((1, 1), (0, 1))
    assert Ray((1, 1), angle=3.0*pi) == Ray((1, 1), (0, 1))
    assert Ray((1, 1), angle=4.0*pi) == Ray((1, 1), (2, 1))
    assert Ray((1, 1), angle=0) == Ray((1, 1), (2, 1))
    assert Ray((1, 1), angle=4.05*pi) == Ray(Point(1, 1),
               Point(2, -sqrt(5)*sqrt(2*sqrt(5) + 10)/4 - sqrt(2*sqrt(5) + 10)/4 + 2 + sqrt(5)))
    assert Ray((1, 1), angle=4.02*pi) == Ray(Point(1, 1),
               Point(2, 1 + tan(4.02*pi)))
    assert Ray((1, 1), angle=5) == Ray((1, 1), (2, 1 + tan(5)))
    raises(ValueError, lambda: Ray((1, 1), 1))

    # issue 7963
    r = Ray((0, 0), angle=x)
    assert r.subs(x, 3*pi/4) == Ray((0, 0), (-1, 1))
    assert r.subs(x, 5*pi/4) == Ray((0, 0), (-1, -1))
    assert r.subs(x, -pi/4) == Ray((0, 0), (1, -1))
    assert r.subs(x, pi/2) == Ray((0, 0), (0, 1))
    assert r.subs(x, -pi/2) == Ray((0, 0), (0, -1))

    r1 = Ray(p1, Point(-1, 5))
    r2 = Ray(p1, Point(-1, 1))
    r3 = Ray(p3, p5)
    r4 = Ray(p1, p2)
    r5 = Ray(p2, p1)
    r6 = Ray(Point(0, 1), Point(1, 2))
    r7 = Ray(Point(0.5, 0.5), Point(1, 1))
    assert l1.projection(r1) == Ray(Point(0, 0), Point(2, 2))
    assert l1.projection(r2) == p1
    assert r3 != r1
    t = Symbol('t', real=True)
    assert Ray((1, 1), angle=pi/4).arbitrary_point() == \
        Point(t + 1, t + 1)
    r8 = Ray(Point(0, 0), Point(0, 4))
    r9 = Ray(Point(0, 1), Point(0, -1))
    assert r8.intersection(r9) == [Segment(Point(0, 0), Point(0, 1))]

    s1 = Segment(p1, p2)
    s2 = Segment(p1, p1_1)
    assert s1.midpoint == Point(Rational(1, 2), Rational(1, 2))
    assert s2.length == sqrt( 2*(x1**2) )
    assert Segment((1, 1), (2, 3)).arbitrary_point() == Point(1 + t, 1 + 2*t)
    assert s1.perpendicular_bisector() == \
        Line(Point(1/2, 1/2), Point(3/2, -1/2))
    # intersections
    assert s1.intersection(Line(p6, p9)) == []
    s3 = Segment(Point(0.25, 0.25), Point(0.5, 0.5))
    assert s1.intersection(s3) == [s1]
    assert s3.intersection(s1) == [s3]
    assert r4.intersection(s3) == [s3]
    assert r4.intersection(Segment(Point(2, 3), Point(3, 4))) == []
    assert r4.intersection(Segment(Point(-1, -1), Point(0.5, 0.5))) == \
        [Segment(p1, Point(0.5, 0.5))]
    s3 = Segment(Point(1, 1), Point(2, 2))
    assert s1.intersection(s3) == [Point(1, 1)]
    s3 = Segment(Point(0.5, 0.5), Point(1.5, 1.5))
    assert s1.intersection(s3) == [Segment(Point(0.5, 0.5), p2)]
    assert s1.intersection(Segment(Point(4, 4), Point(5, 5))) == []
    assert s1.intersection(Segment(Point(-1, -1), p1)) == [p1]
    assert s1.intersection(Segment(Point(-1, -1), Point(0.5, 0.5))) == \
        [Segment(p1, Point(0.5, 0.5))]
    assert r4.intersection(r5) == [s1]
    assert r5.intersection(r6) == []
    assert r4.intersection(r7) == r7.intersection(r4) == [r7]

    # Segment contains
    a, b = symbols('a,b')
    s = Segment((0, a), (0, b))
    assert Point(0, (a + b)/2) in s
    s = Segment((a, 0), (b, 0))
    assert Point((a + b)/2, 0) in s

    raises(Undecidable, lambda: Point(2*a, 0) in s)

    # Testing distance from a Segment to an object
    s1 = Segment(Point(0, 0), Point(1, 1))
    s2 = Segment(Point(half, half), Point(1, 0))
    pt1 = Point(0, 0)
    pt2 = Point(Rational(3)/2, Rational(3)/2)
    assert s1.distance(pt1) == 0
    assert s1.distance((0, 0)) == 0
    assert s2.distance(pt1) == 2**(half)/2
    assert s2.distance(pt2) == 2**(half)
    # Line to point
    p1, p2 = Point(0, 0), Point(1, 1)
    s = Line(p1, p2)
    assert s.distance(Point(-1, 1)) == sqrt(2)
    assert s.distance(Point(1, -1)) == sqrt(2)
    assert s.distance(Point(2, 2)) == 0
    assert s.distance((-1, 1)) == sqrt(2)
    assert Line((0, 0), (0, 1)).distance(p1) == 0
    assert Line((0, 0), (0, 1)).distance(p2) == 1
    assert Line((0, 0), (1, 0)).distance(p1) == 0
    assert Line((0, 0), (1, 0)).distance(p2) == 1
    m = symbols('m')
    l = Line((0, 5), slope=m)
    p = Point(2, 3)
    assert l.distance(p) == 2*abs(m + 1)/sqrt(m**2 + 1)
    # Ray to point
    r = Ray(p1, p2)
    assert r.distance(Point(-1, -1)) == sqrt(2)
    assert r.distance(Point(1, 1)) == 0
    assert r.distance(Point(-1, 1)) == sqrt(2)
    assert Ray((1, 1), (2, 2)).distance(Point(1.5, 3)) == 3*sqrt(2)/4
    assert r.distance((1, 1)) == 0

    #Line contains
    p1, p2 = Point(0, 1), Point(3, 4)
    l = Line(p1, p2)
    assert l.contains(p1) is True
    assert l.contains((0, 1)) is True
    assert l.contains((0, 0)) is False

    #Ray contains
    p1, p2 = Point(0, 0), Point(4, 4)
    r = Ray(p1, p2)
    assert r.contains(p1) is True
    assert r.contains((1, 1)) is True
    assert r.contains((1, 3)) is False
    s = Segment((1, 1), (2, 2))
    assert r.contains(s) is True
    s = Segment((1, 2), (2, 5))
    assert r.contains(s) is False
    r1 = Ray((2, 2), (3, 3))
    assert r.contains(r1) is True
    r1 = Ray((2, 2), (3, 5))
    assert r.contains(r1) is False


    # Special cases of projection and intersection
    r1 = Ray(Point(1, 1), Point(2, 2))
    r2 = Ray(Point(2, 2), Point(0, 0))
    r3 = Ray(Point(1, 1), Point(-1, -1))
    r4 = Ray(Point(0, 4), Point(-1, -5))
    r5 = Ray(Point(2, 2), Point(3, 3))
    assert intersection(r1, r2) == [Segment(Point(1, 1), Point(2, 2))]
    assert intersection(r1, r3) == [Point(1, 1)]
    assert r1.projection(r3) == Point(1, 1)
    assert r1.projection(r4) == Segment(Point(1, 1), Point(2, 2))

    r5 = Ray(Point(0, 0), Point(0, 1))
    r6 = Ray(Point(0, 0), Point(0, 2))
    assert r5 in r6
    assert r6 in r5

    s1 = Segment(Point(0, 0), Point(2, 2))
    s2 = Segment(Point(-1, 5), Point(-5, -10))
    s3 = Segment(Point(0, 4), Point(-2, 2))
    assert intersection(r1, s1) == [Segment(Point(1, 1), Point(2, 2))]
    assert r1.projection(s2) == Segment(Point(1, 1), Point(2, 2))
    assert s3.projection(r1) == Segment(Point(0, 4), Point(-1, 3))

    l1 = Line(Point(0, 0), Point(3, 4))
    r1 = Ray(Point(0, 0), Point(3, 4))
    s1 = Segment(Point(0, 0), Point(3, 4))
    assert intersection(l1, l1) == [l1]
    assert intersection(l1, r1) == [r1]
    assert intersection(l1, s1) == [s1]
    assert intersection(r1, l1) == [r1]
    assert intersection(s1, l1) == [s1]

    entity1 = Segment(Point(-10, 10), Point(10, 10))
    entity2 = Segment(Point(-5, -5), Point(-5, 5))
    assert intersection(entity1, entity2) == []

    r1 = Ray(p1, Point(0, 1))
    r2 = Ray(Point(0, 1), p1)
    r3 = Ray(p1, p2)
    r4 = Ray(p2, p1)
    s1 = Segment(p1, Point(0, 1))
    assert Line(r1.source, r1.random_point()).slope == r1.slope
    assert Line(r2.source, r2.random_point()).slope == r2.slope
    assert Segment(Point(0, -1), s1.random_point()).slope == s1.slope
    p_r3 = r3.random_point()
    p_r4 = r4.random_point()
    assert p_r3.x >= p1.x and p_r3.y >= p1.y
    assert p_r4.x <= p2.x and p_r4.y <= p2.y
    p10 = Point(2000, 2000)
    s1 = Segment(p1, p10)
    p_s1 = s1.random_point()
    assert p1.x <= p_s1.x and p_s1.x <= p10.x and \
        p1.y <= p_s1.y and p_s1.y <= p10.y
    s2 = Segment(p10, p1)
    assert hash(s1) == hash(s2)
    p11 = p10.scale(2, 2)
    assert s1.is_similar(Segment(p10, p11))
    assert s1.is_similar(r1) is False
    assert (r1 in s1) is False
    assert Segment(p1, p2) in s1
    assert s1.plot_interval() == [t, 0, 1]
    assert s1 in Line(p1, p10)
    assert Line(p1, p10) != Line(p10, p1)
    assert Line(p1, p10) != p1
    assert Line(p1, p10).plot_interval() == [t, -5, 5]
    assert Ray((0, 0), angle=pi/4).plot_interval() == \
        [t, 0, 10]

def test_line3d():
    p1 = Point3D(0, 0, 0)
    p2 = Point3D(1, 1, 1)
    p3 = Point3D(x1, x1, x1)
    p4 = Point3D(y1, y1, y1)
    p5 = Point3D(x1, 1 + x1, 1)
    p6 = Point3D(1, 0, 1)
    p7 = Point3D(0, 1, 1)
    p8 = Point3D(2, 0, 3)
    p9 = Point3D(2, 1, 4)

    l1 = Line3D(p1, p2)
    l2 = Line3D(p3, p4)
    l3 = Line3D(p3, p5)
    l4 = Line3D(p1, p6)
    l5 = Line3D(p1, p7)
    l6 = Line3D(p8, p9)
    l7 = Line3D(p2, p9)
    raises(ValueError, lambda: Line3D(Point3D(0, 0, 0), Point3D(0, 0, 0)))

    assert Line3D((1, 1, 1), direction_ratio=[2, 3, 4]) == \
        Line3D(Point3D(1, 1, 1), Point3D(3, 4, 5))
    assert Line3D((1, 1, 1), direction_ratio=[1, 5, 7 ]) == \
        Line3D(Point3D(1, 1, 1), Point3D(2, 6, 8))
    assert Line3D((1, 1, 1), direction_ratio=[1, 2, 3]) == \
        Line3D(Point3D(1, 1, 1), Point3D(2, 3, 4))
    raises(TypeError, lambda: Line3D((1, 1), 1))
    assert Line3D(p1, p2) != Line3D(p2, p1)
    assert l1 != l3
    assert l1.is_parallel(l1)  # same as in 2D
    assert l1 != l2
    assert l1.direction_ratio == [1, 1, 1]
    assert l1.length == oo
    assert l1.equation() == (x, y, z, k)
    assert l2.equation() == \
        ((x - x1)/(-x1 + y1), (-x1 + y)/(-x1 + y1), (-x1 + z)/(-x1 + y1), k)
    assert p1 in l1
    assert p1 not in l3

    # Orthogonality
    p1_1 = Point3D(x1, x1, x1)
    l1_1 = Line3D(p1, p1_1)
    assert Line3D.is_perpendicular(l1, l2) is False
    p = l1.arbitrary_point()
    raises(NotImplementedError , lambda: l1.perpendicular_segment(p))

    # Parallelity
    assert l1.parallel_line(p1_1) == Line3D(Point3D(x1, x1, x1),
        Point3D(x1 + 1, x1 + 1, x1 + 1))
    assert l1.parallel_line(p1_1.args) == \
        Line3D(Point3D(x1, x1, x1), Point3D(x1 + 1, x1 + 1, x1 + 1))

    # Intersection
    assert intersection(l1, p1) == [p1]
    assert intersection(l1, p5) == []
    assert intersection(l1, l1.parallel_line(p1)) == [
        Line3D(Point3D(0, 0, 0), Point3D(1, 1, 1))]
    # issue 8517
    line3 = Line3D(Point3D(4, 0, 1), Point3D(0, 4, 1))
    line4 = Line3D(Point3D(0, 0, 1), Point3D(4, 4, 1))
    assert line3.intersection(line4) == [Point3D(2, 2, 1)]
    assert line3.is_parallel(line4) is False
    assert Line3D((0, 1, 2), (0, 2, 3)).intersection(
        Line3D((0, 1, 2), (0, 1, 1))) == []
    ray0 = Ray3D((0, 0), (3, 0))
    ray1 = Ray3D((1, 0), (3, 0))
    assert ray0.intersection(ray1) == [ray1]
    assert ray1.intersection(ray0) == [ray1]
    assert Segment3D((0, 0), (3, 0)).intersection(
        Segment3D((1, 0), (2, 0))) == [Segment3D((1, 0), (2, 0))]
    assert Segment3D((1, 0), (2, 0)).intersection(
        Segment3D((0, 0), (3, 0))) == [Segment3D((1, 0), (2, 0))]
    assert Segment3D((0, 0), (3, 0)).intersection(
        Segment3D((3, 0), (4, 0))) == [Point3D((3, 0))]
    assert Segment3D((0, 0), (3, 0)).intersection(
        Segment3D((2, 0), (5, 0))) == [Segment3D((3, 0), (2, 0))]
    assert Segment3D((0, 0), (3, 0)).intersection(
        Segment3D((-2, 0), (1, 0))) == [Segment3D((0, 0), (1, 0))]
    assert Segment3D((0, 0), (3, 0)).intersection(
        Segment3D((-2, 0), (0, 0))) == [Point3D(0, 0, 0)]
    # issue 7757
    p = Ray3D(Point3D(1, 0, 0), Point3D(-1, 0, 0))
    q = Ray3D(Point3D(0, 1, 0), Point3D(0, -1, 0))
    assert intersection(p, q) == [Point3D(0, 0, 0)]

    # Concurrency
    assert Line3D.are_concurrent(l1) is False
    assert Line3D.are_concurrent(l1, l2)
    assert Line3D.are_concurrent(l1, l1_1, l3) is False
    parallel_1 = Line3D(Point3D(0, 0, 0), Point3D(1, 0, 0))
    parallel_2 = Line3D(Point3D(0, 1, 0), Point3D(1, 1, 0))
    assert Line3D.are_concurrent(parallel_1, parallel_2) == False

    # Finding angles
    l1_1 = Line3D(p1, Point3D(5, 0, 0))
    assert Line3D.angle_between(l1, l1_1), acos(sqrt(3)/3)

    # Testing Rays and Segments (very similar to Lines)
    assert Ray3D((1, 1, 1), direction_ratio=[4, 4, 4]) == \
        Ray3D(Point3D(1, 1, 1), Point3D(5, 5, 5))
    assert Ray3D((1, 1, 1), direction_ratio=[1, 2, 3]) == \
        Ray3D(Point3D(1, 1, 1), Point3D(2, 3, 4))
    assert Ray3D((1, 1, 1), direction_ratio=[1, 1, 1]) == \
        Ray3D(Point3D(1, 1, 1), Point3D(2, 2, 2))

    r1 = Ray3D(p1, Point3D(-1, 5, 0))
    r2 = Ray3D(p1, Point3D(-1, 1, 1))
    r3 = Ray3D(p1, p2)
    r4 = Ray3D(p2, p1)
    r5 = Ray3D(Point3D(0, 1, 1), Point3D(1, 2, 0))
    assert l1.projection(r1) == [
        Ray3D(Point3D(0, 0, 0), Point3D(4/3, 4/3, 4/3))]
    assert l1.projection(r2) == [
        Ray3D(Point3D(0, 0, 0), Point3D(1/3, 1/3, 1/3))]
    assert r3 != r1
    t = Symbol('t', real=True)
    assert Ray3D((1, 1, 1), direction_ratio=[1, 2, 3]).arbitrary_point() == \
        Point3D(t + 1, 2*t + 1, 3*t + 1)
    r6 = Ray3D(Point3D(0, 0, 0), Point3D(0, 4, 0))
    r7 = Ray3D(Point3D(0, 1, 1), Point3D(0, -1, 1))
    assert r6.intersection(r7) == []

    s1 = Segment3D(p1, p2)
    s2 = Segment3D(p3, p4)
    assert s1.midpoint == \
        Point3D(Rational(1, 2), Rational(1, 2), Rational(1, 2))
    assert s2.length == sqrt(3)*sqrt((x1 - y1)**2)
    assert Segment3D((1, 1, 1), (2, 3, 4)).arbitrary_point() == \
        Point3D(t + 1, 2*t + 1, 3*t + 1)

    # Segment contains
    s = Segment3D((0, 1, 0), (0, 1, 0))
    assert Point3D(0, 1, 0) in s
    s = Segment3D((1, 0, 0), (1, 0, 0))
    assert Point3D(1, 0, 0) in s

    # Testing distance from a Segment to an object
    s1 = Segment3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
    s2 = Segment3D(Point3D(1/2, 1/2, 1/2), Point3D(1, 0, 1))
    pt1 = Point3D(0, 0, 0)
    pt2 = Point3D(Rational(3)/2, Rational(3)/2, Rational(3)/2)
    assert s1.distance(pt1) == 0
    assert s2.distance(pt1) == sqrt(3)/2
    assert s2.distance(pt2) == 2
    assert s1.distance((0,0,0)) == 0
    assert s2.distance((0,0,0)) == sqrt(3)/2
    # Line to point
    p1, p2 = Point3D(0, 0, 0), Point3D(1, 1, 1)
    s = Line3D(p1, p2)
    assert s.distance(Point3D(-1, 1, 1)) == 2*sqrt(6)/3
    assert s.distance(Point3D(1, -1, 1)) == 2*sqrt(6)/3
    assert s.distance(Point3D(2, 2, 2)) == 0
    assert s.distance((2, 2, 2)) == 0
    assert s.distance((1, -1, 1)) == 2*sqrt(6)/3
    assert Line3D((0, 0, 0), (0, 1, 0)).distance(p1) == 0
    assert Line3D((0, 0, 0), (0, 1, 0)).distance(p2) == sqrt(2)
    assert Line3D((0, 0, 0), (1, 0, 0)).distance(p1) == 0
    assert Line3D((0, 0, 0), (1, 0, 0)).distance(p2) == sqrt(2)
    # Ray to point
    r = Ray3D(p1, p2)
    assert r.distance(Point3D(-1, -1, -1)) == sqrt(3)
    assert r.distance(Point3D(1, 1, 1)) == 0
    assert r.distance((-1, -1, -1)) == sqrt(3)
    assert r.distance((1, 1, 1)) == 0
    assert Ray3D((1, 1, 1), (2, 2, 2)).distance(Point3D(1.5, 3, 1)) == \
        sqrt(17)/2


    # Special cases of projection and intersection
    r1 = Ray3D(Point3D(1, 1, 1), Point3D(2, 2, 2))
    r2 = Ray3D(Point3D(2, 2, 2), Point3D(0, 0, 0))
    r3 = Ray3D(Point3D(1, 1, 1), Point3D(-1, -1, -1))
    r4 = Ray3D(Point3D(0, 4, 2), Point3D(-1, -5, -1))
    r5 = Ray3D(Point3D(2, 2, 2), Point3D(3, 3, 3))
    assert intersection(r1, r2) == \
        [Segment3D(Point3D(1, 1, 1), Point3D(2, 2, 2))]
    assert intersection(r1, r3) == [Point3D(1, 1, 1)]

    r5 = Ray3D(Point3D(0, 0, 0), Point3D(1, 1, 1))
    r6 = Ray3D(Point3D(0, 0, 0), Point3D(2, 2, 2))
    assert r5 in r6
    assert r6 in r5

    s1 = Segment3D(Point3D(0, 0, 0), Point3D(2, 2, 2))
    s2 = Segment3D(Point3D(-1, 5, 2), Point3D(-5, -10, 0))
    assert intersection(r1, s1) == [
        Segment3D(Point3D(1, 1, 1), Point3D(2, 2, 2))]

    l1 = Line3D(Point3D(0, 0, 0), Point3D(3, 4, 0))
    r1 = Ray3D(Point3D(0, 0, 0), Point3D(3, 4, 0))
    s1 = Segment3D(Point3D(0, 0, 0), Point3D(3, 4, 0))
    assert intersection(l1, r1) == [r1]
    assert intersection(l1, s1) == [s1]
    assert intersection(r1, l1) == [r1]
    assert intersection(s1, r1) == [s1]

    # check that temporary symbol is Dummy
    assert Line3D((0, 0), (t, t)).perpendicular_line((0, 1)) == \
        Line3D(Point3D(0, 1, 0), Point3D(1/2, 1/2, 0))
    assert Line3D((0, 0), (t, t)).perpendicular_segment((0, 1)) == \
        Segment3D(Point3D(0, 1, 0), Point3D(1/2, 1/2, 0))
    assert Line3D((0, 0), (t, t)).intersection(Line3D((0, 1), (t, t))) == \
        [Point3D(t, t, 0)]
    assert Line3D((0, 0, 0), (x, y, z)).contains((2*x, 2*y, 2*z))

    # Test is_perpendicular
    perp_1 = Line3D(p1, Point3D(0, 1, 0))
    assert Line3D.is_perpendicular(parallel_1, perp_1) is True
    assert Line3D.is_perpendicular(parallel_1, parallel_2) is False

    # Test projection
    assert parallel_1.projection(Point3D(5, 5, 0)) == Point3D(5, 0, 0)
    assert parallel_1.projection(parallel_2) == [parallel_1]
    raises(GeometryError, lambda: parallel_1.projection(Plane(p1, p2, p6)))

    # Test __new__
    assert Line3D(perp_1) == perp_1
    raises(ValueError, lambda: Line3D(p1))

    # Test contains
    pt2d = Point(1.0, 1.0)
    assert perp_1.contains(pt2d) is False

    # Test equals
    assert perp_1.equals(pt2d) is False
    col1 = Line3D(Point3D(0, 0, 0), Point3D(1, 0, 0))
    col2 = Line3D(Point3D(-5, 0, 0), Point3D(-1, 0, 0))
    assert col1.equals(col2) is True
    assert col1.equals(perp_1) is False

    # Begin ray
    # Test __new__
    assert Ray3D(col1) == Ray3D(p1, Point3D(1, 0, 0))
    raises(ValueError, lambda: Ray3D(pt2d))

    # Test zdirection
    negz = Ray3D(p1, Point3D(0, 0, -1))
    assert negz.zdirection == S.NegativeInfinity

    # Test contains
    assert negz.contains(Segment3D(p1, Point3D(0, 0, -10))) is True
    assert negz.contains(Segment3D(Point3D(1, 1, 1), Point3D(2, 2, 2))) is False
    posy = Ray3D(p1, Point3D(0, 1, 0))
    posz = Ray3D(p1, Point3D(0, 0, 1))
    assert posy.contains(p1) is True
    assert posz.contains(p1) is True
    assert posz.contains(pt2d) is False
    ray1 = Ray3D(Point3D(1, 1, 1), Point3D(1, 0, 0))
    raises(TypeError, lambda: ray1.contains([]))

    # Test equals
    assert negz.equals(pt2d) is False
    assert negz.equals(negz) is True

    assert ray1.is_similar(Line3D(Point3D(1, 1, 1), Point3D(1, 0, 0))) is True
    assert ray1.is_similar(perp_1) is False
    raises(NotImplementedError, lambda: ray1.is_similar(ray1))

    # Begin Segment
    seg1 = Segment3D(p1, Point3D(1, 0, 0))
    raises(TypeError, lambda: seg1.contains([]))
    seg2= Segment3D(Point3D(2, 2, 2), Point3D(3, 2, 2))
    assert seg1.contains(seg2) is False

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


def test_line_intersection():
    assert asa(120, 8, 52) == \
        Triangle(
            Point(0, 0),
            Point(8, 0),
            Point(-4*cos(19*pi/90)/sin(2*pi/45),
            4*sqrt(3)*cos(19*pi/90)/sin(2*pi/45)))
    assert Line((0, 0), (1, 1)).intersection(Ray((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert Line((0, 0), (1, 1)).intersection(Segment((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert Ray((0, 0), (1, 1)).intersection(Ray((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert Ray((0, 0), (1, 1)).intersection(Segment((1, 0), (1, 2))) == \
        [Point(1, 1)]
    assert Ray((0, 0), (10, 10)).contains(Segment((1, 1), (2, 2))) is True
    assert Segment((1, 1), (2, 2)) in Line((0, 0), (10, 10))
    x = 8*tan(13*pi/45)/(tan(13*pi/45) + sqrt(3))
    y = (-8*sqrt(3)*tan(13*pi/45)**2 + 24*tan(13*pi/45))/ \
        (-3 + tan(13*pi/45)**2)
    assert Line(Point(0, 0), Point(1, -sqrt(3))).contains(Point(x, y)) is True


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


def test_issue_2941():
    def _check():
        for f, g in cartes(*[(Line, Ray, Segment)]*2):
            l1 = f(a, b)
            l2 = g(c, d)
            assert l1.intersection(l2) == l2.intersection(l1)
    # intersect at end point
    c, d = (-2, -2), (-2, 0)
    a, b = (0, 0), (1, 1)
    _check()
    # midline intersection
    c, d = (-2, -3), (-2, 0)
    a, b = (0, 0), (1, 1)
    _check()


def test_symbolic_intersect():
    # Issue 7814.
    circle = Circle(Point(x, 0), y)
    line = Line(Point(k, z), slope=0)
    assert line.intersection(circle) == [
        Point(x - sqrt(y**2 - z**2), z), Point(x + sqrt(y**2 - z**2), z)]
