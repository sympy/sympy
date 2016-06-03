from __future__ import division

from sympy import Dummy, Rational, S, Symbol, pi, sqrt, oo
from sympy.core.compatibility import range
from sympy.geometry import (Circle, Ellipse, GeometryError, Line, Point,
                            Point2D, Polygon, Ray, RegularPolygon,
                            Segment, Triangle, intersection, Parabola)
from sympy.integrals.integrals import Integral
from sympy.utilities.pytest import raises, slow


def test_parabola_geom():
    p1 = Point(0, 0)
    p2 = Point(3, 7)
    p3 = Point(0, 4)
    p4 = Point(6, 0)
    d1 = Line(Point(4, 0), Point(4, 9))
    d2 = Line(Point(7, 6), Point(3, 6))
    d3 = Line(Point(4, 0), slope=oo)
    d4 = Line(Point(7, 6), slope=0)

    pa1 = Parabola(None, d2)
    pa2 = Parabola(directrix=d1)
    pa3 = Parabola(p1, d1)
    pa4 = Parabola(p2, d2)
    pa5 = Parabola(p2, d4)
    pa6 = Parabola(p3, d2)
    pa7 = Parabola(p2, d1)
    pa8 = Parabola(p4, d1)
    pa9 = Parabola(p4, d3)

    # Basic Stuff
    assert pa1.focus == Point(0, 0)
    assert pa2 == pa3
    assert pa4 != pa7
    assert pa6 != pa7
    assert pa6.focus == Point2D(0, 4)
    assert pa6.focal_length == 1
    assert pa6.p_parameter == -1
    assert pa6.vertex == Point2D(0, 5)
    assert pa6.eccentricity == 1
    assert pa7.focus == Point2D(3, 7)
    assert pa7.focal_length == 1/2
    assert pa7.p_parameter == -1/2
    assert pa7.vertex == Point2D(7/2, 7)
    assert pa4.focal_length == 1/2
    assert pa4.p_parameter == 1/2
    assert pa4.vertex == Point2D(3, 13/2)
    assert pa8.focal_length == 1
    assert pa8.p_parameter == 1
    assert pa8.vertex == Point2D(5, 0)
    assert pa4.focal_length == pa5.focal_length
    assert pa4.p_parameter == pa5.p_parameter
    assert pa4.vertex == pa5.vertex
    assert pa4.equation() == pa5.equation()
    assert pa8.focal_length == pa9.focal_length
    assert pa8.p_parameter == pa9.p_parameter
    assert pa8.vertex == pa9.vertex
    assert pa8.equation() == pa9.equation()
