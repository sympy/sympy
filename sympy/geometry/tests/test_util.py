from __future__ import division

from sympy import Symbol, sqrt, Derivative
from sympy.geometry import Point, Point2D, Polygon, Segment, convex_hull, intersection, centroid
from sympy.geometry.util import idiff, closest_points, farthest_points
from sympy.solvers.solvers import solve
from sympy.utilities.pytest import raises


def test_idiff():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)
    t = Symbol('t', real=True)
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


def test_util():
    # coverage for some leftover functions in sympy.geometry.util
    assert intersection(Point(0, 0)) == []
    raises(ValueError, lambda: intersection(Point(0, 0), 3))
    raises(ValueError, lambda: convex_hull(Point(0, 0), 3))


def test_util_centroid():
    p = Polygon((0, 0), (10, 0), (10, 10))
    q = p.translate(0, 20)
    assert centroid(p, q) == Point(20, 40)/3
    p = Segment((0, 0), (2, 0))
    q = Segment((0, 0), (2, 2))
    assert centroid(p, q) == Point(1, -sqrt(2) + 2)
    assert centroid(Point(0, 0), Point(2, 0)) == Point(2, 0)/2
    assert centroid(Point(0, 0), Point(0, 0), Point(2, 0)) == Point(2, 0)/3


def test_closest_points():
    points = [Point2D(1, 1), Point2D(2, 2), Point2D(11, 11)]
    assert closest_points(*points) == (Point2D(1, 1), Point2D(2, 2))
    points = [Point2D(38, 40), Point2D(83, 89), Point2D(74, 46), Point2D(76, 73), Point2D(75, 13), Point2D(98, 87), Point2D(5, 73), Point2D(75, 18), Point2D(30, 7)]
    assert closest_points(*points) == (Point2D(75, 18), Point2D(75, 13))
    points = [Point2D(34, 54), Point2D(92, 10), Point2D(99, 2), Point2D(29, 4), Point2D(4, 51), Point2D(2, 40), Point2D(18, 60), Point2D(30, 69), Point2D(10, 32)]
    assert closest_points(*points) == (Point2D(99, 2), Point2D(92, 10))


def test_farthest_points():
    points = [Point2D(1, 1), Point2D(2, 2), Point2D(11, 11)]
    assert farthest_points(*points) == (Point2D(1, 1), Point2D(11, 11))
    points = [Point2D(38, 40), Point2D(83, 89), Point2D(74, 46), Point2D(76, 73), Point2D(75, 13), Point2D(98, 87), Point2D(5, 73), Point2D(75, 18), Point2D(30, 7)]
    assert farthest_points(*points) == (Point2D(98, 87), Point2D(30, 7))
    points = [Point2D(34, 54), Point2D(92, 10), Point2D(99, 2), Point2D(29, 4), Point2D(4, 51), Point2D(2, 40), Point2D(18, 60), Point2D(30, 69), Point2D(10, 32)]
    assert farthest_points(*points) == (Point2D(4, 51), Point2D(99, 2))
