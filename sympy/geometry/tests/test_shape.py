from __future__ import division

from sympy import Rational, oo, sqrt
from sympy import Line, Point, Point2D, Parabola, Segment2D, Ray2D, Line2D
from sympy import Circle, Ellipse, Parabola, Shape
from sympy.utilities.pytest import raises
from sympy import S

def test_shape_geom():
    e = Ellipse((0, 0), 4, 2)
    p_h = Parabola((1, 0), Line((-3, 0), (-3, 2)))
    c = Circle((0, 0), 4)
    l_v = Line((0, 0), slope=oo)
    l_h = Line((0, 0), slope=0)

    s1 = Shape(e, l_h)
    s2 = Shape(e, l_v)
    s3 = Shape(p_h, l_v)
    s4 = Shape(c, l_h)
    s5 = Shape(c, l_v)
    raises(TypeError, lambda:
           Shape(p_h, l_h))

    assert s1 != s2
    assert s4 != s5
    assert s1.conic == Ellipse(Point2D(0, 0), 4, 2)
    assert s3.conic == Parabola(Point2D(1, 0), Line2D(Point2D(-3, 0), Point2D(-3, 2)))
    assert s4.conic == Circle(Point2D(0, 0), 4)
    assert s1.area == S.Pi*4*2/2
    assert s2.area == S.Pi*4*2/2
    assert s4.area == S.Pi*4*4/2
    assert s5.area == S.Pi*4*4/2
    assert s1.area == s2.area
    assert s4.area == s5.area
    assert s1.centroid == Point2D(0, 8/(3*S.Pi))
    assert s2.centroid == Point2D(16/(3*S.Pi), 0)
    assert s3.centroid == Point2D(-2/5, 0)
    assert s4.centroid == Point2D(0, 16/(3*S.Pi))
    assert s5.centroid == Point2D(16/(3*S.Pi), 0)
