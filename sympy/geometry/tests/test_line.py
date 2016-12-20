from __future__ import division

from sympy import Rational, Float, S, Symbol, cos, oo, pi, simplify, sin, sqrt, symbols, acos
from sympy.core.compatibility import range
from sympy.functions.elementary.trigonometric import tan
from sympy.geometry import (Circle, GeometryError, Line, Point, Ray, Segment, Triangle, intersection, Point3D, Line3D,
                            Ray3D, Segment3D, Plane, Point2D, Line2D, Ray2D, Segment2D)
from sympy.geometry.line import Undecidable
from sympy.geometry.polygon import _asa as asa
from sympy.utilities.iterables import cartes
from sympy.utilities.pytest import raises, slow

import traceback
import warnings
import sys

# make warnings show tracebacks
def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    traceback.print_stack()
    log = file if hasattr(file,'write') else sys.stderr
    log.write(warnings.formatwarning(message, category, filename, lineno, line))
warnings.showwarning = warn_with_traceback
warnings.simplefilter('always', UserWarning)     # make sure to show warnings every time they occurr

def feq(a, b):
    """Test if two floating point values are 'equal'."""
    t_float = Float("1.0E-10")
    return -t_float < a - b < t_float


def test_line_geom():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    x1 = Symbol('x1', real=True)
    y1 = Symbol('y1', real=True)
    half = Rational(1, 2)
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
    raises(TypeError, lambda: Line((1, 1), 1))
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
    assert l1.perpendicular_line(p1.args).equals( Line(Point(0, 0), Point(1, -1)) )
    assert l1.perpendicular_line(p1).equals( Line(Point(0, 0), Point(1, -1)) )
    assert Line.is_perpendicular(l1, l1_1)
    assert Line.is_perpendicular(l1, l2) is False
    p = l1.random_point()
    assert l1.perpendicular_segment(p) == p

    # Parallelity
    l2_1 = Line(p3, p5)
    assert l2.parallel_line(p1_1).equals( Line(Point(-x1, x1), Point(-y1, 2*x1 - y1)) )
    assert l2_1.parallel_line(p1.args).equals( Line(Point(0, 0), Point(0, -1)) )
    assert l2_1.parallel_line(p1).equals( Line(Point(0, 0), Point(0, -1)) )
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
    assert Line.are_concurrent(l1, l1, l1, l3)
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
    a = Point(1, 2, 3, 4)
    b = a.orthogonal_direction
    o = a.origin
    assert Line(a, o).angle_between(Line(b, o)) == pi/2

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
    raises(TypeError, lambda: Ray((1, 1), 1))

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
    aline = Line(Point(1/2, 1/2), Point(3/2, -1/2))
    assert s1.perpendicular_bisector().equals(aline)
    on_line = Segment(Point(1/2, 1/2), Point(3/2, -1/2)).midpoint
    assert s1.perpendicular_bisector(on_line) == Segment(s1.midpoint, on_line)
    assert s1.perpendicular_bisector(on_line + (1, 0)).equals(aline)
    # intersections
    assert s1.intersection(Line(p6, p9)) == []
    s3 = Segment(Point(0.25, 0.25), Point(0.5, 0.5))
    assert s1.intersection(s3) == [s3]
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
    a, b = symbols('a,b', real=True)
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
    m = symbols('m', real=True)
    l = Line((0, 5), slope=m)
    p = Point(2, 3)
    assert (l.distance(p) - 2*abs(m + 1)/sqrt(m**2 + 1)).equals(0)
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
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    z = Symbol('z', real=True)
    k = Symbol('k', real=True)
    x1 = Symbol('x1', real=True)
    y1 = Symbol('y1', real=True)
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
    assert l1.perpendicular_segment(p) == p

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
        Line3D((0, 1, 2), (0, 1, 1))) == [Point3D(0,1,2)]
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
        Segment3D((-2, 0), (0, 0))) == [Point3D(0, 0)]
    # issue 7757
    p = Ray3D(Point3D(1, 0, 0), Point3D(-1, 0, 0))
    q = Ray3D(Point3D(0, 1, 0), Point3D(0, -1, 0))
    assert intersection(p, q) == [Point3D(0, 0, 0)]

    # Concurrency
    assert Line3D.are_concurrent(l1) is False
    assert Line3D.are_concurrent(l1, l2) is False
    assert Line3D.are_concurrent(l1, l1_1, l3) is True
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
    assert l1.projection(r1) == Ray3D(Point3D(0, 0, 0), Point3D(4/3, 4/3, 4/3))
    assert l1.projection(r2) == Ray3D(Point3D(0, 0, 0), Point3D(1/3, 1/3, 1/3))
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
    assert s2.distance(pt2) == 2*sqrt(6)/3
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
    assert Ray3D((0,0,0),(1,1,2)).distance((-1,-1,2)) == 4*sqrt(3)/3
    assert Ray3D((1, 1, 1), (2, 2, 2)).distance(Point3D(1.5, -3, -1)) == \
        Rational(9)/2
    assert Ray3D((1, 1, 1), (2, 2, 2)).distance(Point3D(1.5, 3, 1)) == \
        sqrt(78)/6


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
    assert Line3D((0, 0), (t, t)).perpendicular_line((0, 1, 0)).equals( \
        Line3D(Point3D(0, 1, 0), Point3D(1/2, 1/2, 0)))
    assert Line3D((0, 0), (t, t)).perpendicular_segment((0, 1, 0)).equals( \
        Segment3D((0, 1), (1/2, 1/2)))
    assert Line3D((0, 0), (t, t)).intersection(Line3D((0, 1), (t, t))) == \
        [Point3D(t, t)]
    assert Line3D((0, 0, 0), (x, y, z)).contains((2*x, 2*y, 2*z))

    # Test is_perpendicular
    perp_1 = Line3D(p1, Point3D(0, 1, 0))
    assert Line3D.is_perpendicular(parallel_1, perp_1) is True
    assert Line3D.is_perpendicular(parallel_1, parallel_2) is False

    # Test projection
    assert parallel_1.projection(Point3D(5, 5, 0)) == Point3D(5, 0)
    assert parallel_1.projection(parallel_2).equals( parallel_1 )
    raises(GeometryError, lambda: parallel_1.projection(Plane(p1, p2, p6)))

    # Test __new__
    assert Line3D(perp_1) == perp_1
    raises(ValueError, lambda: Line3D(p1))

    # Test contains
    pt2d = Point(1.0, 1.0)
    with warnings.catch_warnings(record=True) as w:
        assert perp_1.contains(pt2d) is False
        assert len(w) == 1

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
    with warnings.catch_warnings(record=True) as w:
        assert posz.contains(pt2d) is False
        assert len(w) == 1
    ray1 = Ray3D(Point3D(1, 1, 1), Point3D(1, 0, 0))
    assert ray1.contains([]) is False

    # Test equals
    assert negz.equals(pt2d) is False
    assert negz.equals(negz) is True

    assert ray1.is_similar(Line3D(Point3D(1, 1, 1), Point3D(1, 0, 0))) is True
    assert ray1.is_similar(perp_1) is False
    assert ray1.is_similar(ray1) is True

    # Begin Segment
    seg1 = Segment3D(p1, Point3D(1, 0, 0))
    assert seg1.contains([]) is True
    seg2= Segment3D(Point3D(2, 2, 2), Point3D(3, 2, 2))
    assert seg1.contains(seg2) is False

@slow
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
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    z = Symbol('z', real=True)
    k = Symbol('k', real=True)
    # Issue 7814.
    circle = Circle(Point(x, 0), y)
    line = Line(Point(k, z), slope=0)
    assert line.intersection(circle) == [
        Point(x - sqrt(y**2 - z**2), z), Point(x + sqrt(y**2 - z**2), z)]


def test_arguments():
    """Functions accepting `Point` objects in `geometry`
    should also accept tuples, lists, and generators and
    automatically convert them to points."""
    from sympy import subsets

    singles2d = ((1,2), [1,3], Point(1,5))
    doubles2d = subsets(singles2d, 2)
    l2d = Line(Point2D(1,2), Point2D(2,3))
    singles3d = ((1,2,3), [1,2,4], Point(1,2,6))
    doubles3d = subsets(singles3d, 2)
    l3d = Line(Point3D(1,2,3), Point3D(1, 1, 2))
    singles4d = ((1,2,3,4), [1,2,3,5], Point(1,2,3,7))
    doubles4d = subsets(singles4d, 2)
    l4d = Line(Point(1,2,3,4), Point(2,2,2,2))

    # test 2D
    test_single = ['contains', 'distance', 'equals', 'parallel_line', 'perpendicular_line', 'perpendicular_segment', 'projection', 'intersection']
    for p in doubles2d:
        Line2D(*p)
    for func in test_single:
        for p in singles2d:
            getattr(l2d, func)(p)

    # test 3D
    for p in doubles3d:
        Line3D(*p)
    for func in test_single:
        for p in singles3d:
            getattr(l3d, func)(p)

    # test 4D
    for p in doubles4d:
        Line(*p)
    for func in test_single:
        for p in singles4d:
            getattr(l4d, func)(p)

def test_dimension_normalization():
    with warnings.catch_warnings(record=True) as w:
        assert Ray((1, 1), (2, 1, 2)) == Ray((1, 1, 0), (2, 1, 2))
        assert len(w) == 1
