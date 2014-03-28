from __future__ import division
from sympy import sqrt
from sympy.physics.optics.point import Point
from sympy.abc import x, y

def test_point_3d():
    p1 = Point(0, 0, 0)
    p1_dash = Point(0, 0, 0)
    p2 = Point(1, 1, 1)
    p3 = Point(2, 2, 2)
    p4 = Point(2, 4, 6)
    p5 = Point(x, y, 4)
    p6 = Point(x, x, x)
    p7 = Point(y, y, y)
    assert p1 == p1_dash
    assert p1 + p2 == p2
    assert p3 - p2 == p2
    assert Point.is_collinear(p1, p2, p3) is True
    assert Point.is_collinear(p1, p2, p3, p4) is False
    assert Point.is_collinear(p1, p2, p3, p6, p7) is True
    assert p2.distance(p4) == sqrt(35)
    assert p4.midpoint(p2) == Point(3/2, 5/2, 7/2)
    p4.translate(2, 2, 7)
    assert p4 == Point(4, 6, 13)


def test_point_2d():
    p1 = Point(0, 0)
    p1_dash = Point(0, 0)
    p2 = Point(1, 1)
    p3 = Point(2, 2)
    p4 = Point(2, 4)
    p5 = Point(x, y)
    p6 = Point(x, x)
    p7 = Point(y, y)
    assert p1 == p1_dash
    assert p1 + p2 == p2
    assert p3 - p2 == p2
    assert Point.is_collinear(p1, p2, p3) is True
    assert Point.is_collinear(p1, p2, p3, p4) is False
    assert Point.is_collinear(p1, p2, p3, p6, p7) is True
    assert p2.distance(p4) == sqrt(10)
    assert p4.midpoint(p2) == Point(3/2, 5/2)
    p4.translate(2, 2)
    assert p4 == Point(4, 6)
    p4.reflect(axis = 'x')
    assert p4 == Point(4, -6)
    p4.reflect(axis = 'y')
    assert p4 == Point(-4, -6)
