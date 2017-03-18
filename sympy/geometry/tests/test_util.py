from __future__ import division

from sympy import Symbol, sqrt, Derivative, symbols, S, Mul
from sympy.geometry import (Line, Point, Point2D, Polygon,
    Segment, convex_hull, intersection, centroid, Circle, Ellipse)
from sympy.geometry.util import (idiff, closest_points,
    farthest_points, _ordered_points, gsolve)
from sympy.solvers.solvers import solve
from sympy.polys.rootoftools import CRootOf
from sympy.utilities.pytest import raises, slow


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
    raises(TypeError, lambda: intersection(Point(0, 0), 3))


def test_convex_hull():
    raises(TypeError, lambda: convex_hull(Point(0, 0), 3))
    points = [(1, -1), (1, -2), (3, -1), (-5, -2), (15, -4)]
    assert convex_hull(*points, **dict(polygon=False)) == (
        [Point2D(-5, -2), Point2D(1, -1), Point2D(3, -1), Point2D(15, -4)],
        [Point2D(-5, -2), Point2D(15, -4)])


def test_util_centroid():
    p = Polygon((0, 0), (10, 0), (10, 10))
    q = p.translate(0, 20)
    assert centroid(p, q) == Point(20, 40)/3
    p = Segment((0, 0), (2, 0))
    q = Segment((0, 0), (2, 2))
    assert centroid(p, q) == Point(1, -sqrt(2) + 2)
    assert centroid(Point(0, 0), Point(2, 0)) == Point(2, 0)/2
    assert centroid(Point(0, 0), Point(0, 0), Point(2, 0)) == Point(2, 0)/3


def test_farthest_points_closest_points():
    from random import randint
    from sympy.utilities.iterables import subsets

    for how in (min, max):
        if how is min:
            func = closest_points
        else:
            func = farthest_points

        raises(ValueError, lambda: func(Point2D(0, 0), Point2D(0, 0)))

        # 3rd pt dx is close and pt is closer to 1st pt
        p1 = [Point2D(0, 0), Point2D(3, 0), Point2D(1, 1)]
        # 3rd pt dx is close and pt is closer to 2nd pt
        p2 = [Point2D(0, 0), Point2D(3, 0), Point2D(2, 1)]
        # 3rd pt dx is close and but pt is not closer
        p3 = [Point2D(0, 0), Point2D(3, 0), Point2D(1, 10)]
        # 3rd pt dx is not closer and it's closer to 2nd pt
        p4 = [Point2D(0, 0), Point2D(3, 0), Point2D(4, 0)]
        # 3rd pt dx is not closer and it's closer to 1st pt
        p5 = [Point2D(0, 0), Point2D(3, 0), Point2D(-1, 0)]
        # duplicate point doesn't affect outcome
        dup = [Point2D(0, 0), Point2D(3, 0), Point2D(3, 0), Point2D(-1, 0)]
        # symbolic
        x = Symbol('x', positive=True)
        s = [Point2D(a) for a in ((x, 1), (x + 3, 2), (x + 2, 2))]

        for points in (p1, p2, p3, p4, p5, s, dup):
            d = how(i.distance(j) for i, j in subsets(points, 2))
            ans = a, b = list(func(*points))[0]
            a.distance(b) == d
            assert ans == _ordered_points(ans)

        # if the following ever fails, the above tests were not sufficient
        # and the logical error in the routine should be fixed
        points = set()
        while len(points) != 7:
            points.add(Point2D(randint(1, 100), randint(1, 100)))
        points = list(points)
        d = how(i.distance(j) for i, j in subsets(points, 2))
        ans = a, b = list(func(*points))[0]
        a.distance(b) == d
        assert ans == _ordered_points(ans)

    # equidistant points
    a, b, c = (
        Point2D(0, 0), Point2D(1, 0), Point2D(1/2, sqrt(3)/2))
    ans = set([_ordered_points((i, j))
        for i, j in subsets((a, b, c), 2)])
    assert closest_points(b, c, a) == ans
    assert farthest_points(b, c, a) == ans

    # unique to farthest
    points = [(1, 1), (1, 2), (3, 1), (-5, 2), (15, 4)]
    assert farthest_points(*points) == set(
        [(Point2D(-5, 2), Point2D(15, 4))])
    points = [(1, -1), (1, -2), (3, -1), (-5, -2), (15, -4)]
    assert farthest_points(*points) == set(
        [(Point2D(-5, -2), Point2D(15, -4))])
    assert farthest_points((1, 1), (0, 0)) == set(
        [(Point2D(0, 0), Point2D(1, 1))])
    raises(ValueError, lambda: farthest_points((1, 1)))


@slow
def test_gsolve():
    a, x, y = symbols('a x y')
    raises(AssertionError, lambda: gsolve(1, 2, x, y))
    raises(AssertionError, lambda: gsolve(a + x, a - y, x, y))
    raises(AssertionError, lambda: gsolve(x + 1/x, x - 2, x, y))
    raises(AssertionError, lambda: gsolve(x - 2, y - 2))
    raises(AttributeError, lambda: gsolve(Segment((0, 1), (1, 2)), x - 2, x, y))
    assert gsolve(Line((0, 1), slope=1), Line((0, 1), slope=-2)
        ) == set([(0, 1)])
    assert gsolve(Circle((0, 0), 2), Circle((0, 0), 1)) == set()
    assert gsolve(x - 2, x - 3, x, y) == set([])
    assert gsolve(y - 1, y - 2, x, y) == set([])
    assert gsolve(x**2 - 4, x - 2, x, y) == set([(2, y)])
    assert gsolve(y**2 - 4, y - 2, x, y) == set([(x, 2)])
    assert gsolve(x - 1, y - 2, x, y) == set([(1, 2)])
    assert gsolve(y - 2, x - 1, x, y) == set([(1, 2)])
    assert gsolve(x - 2, y**2 - x + 4, x, y) == set([])
    assert gsolve(y**2 - x + 4, x - 2, x, y) == set([])
    assert gsolve(y - 2, x**2 - y + 4, x, y) == set([])
    assert gsolve(x**2 - y + 4, y - 2, x, y) == set([])
    assert gsolve(x**2 + y - 2, x - 1, x, y) == set([(1, 1)])
    assert gsolve(x - 1, x**2 + y - 2, x, y) == set([(1, 1)])
    assert gsolve(x**2 + y - 2, y - 1, x, y) == set([(-1, 1), (1, 1)])
    assert gsolve(y - 1, x**2 + y - 2, x, y) == set([(-1, 1), (1, 1)])
    assert gsolve(
        x**2 + 3*x*y + 3*x + y**2 + 3*y + 1, 
        x**2 + 2*x*y + x + y**2 + 2*y + 1, x, y) == set([(1, -1)])
    eqs = (x-2)**2-x*y*3-y**2,x**2-3*x*y-y**2
    assert gsolve(eqs[0], eqs[1], x, y, check=False) == set([
        (-sqrt(-S(3)/2 + sqrt(13)/2)*sqrt(S(9)/2 + 13*sqrt(13)/2)/2 -
        S(1)/4 + 3*sqrt(13)/4, -S(3)/2 + sqrt(13)/2), (-S(1)/4 + sqrt(-S(3)/2 +
        sqrt(13)/2)*sqrt(S(9)/2 + 13*sqrt(13)/2)/2 + 3*sqrt(13)/4, -S(3)/2
        + sqrt(13)/2), (-sqrt(-S(9)/2 + 13*sqrt(13)/2)*sqrt(S(3)/2 +
        sqrt(13)/2)/2 - 3*sqrt(13)/4 - S(1)/4, -sqrt(13)/2 - S(3)/2),
        (-3*sqrt(13)/4 - S(1)/4 + sqrt(-S(9)/2 + 13*sqrt(13)/2)*sqrt(S(3)/2 +
        sqrt(13)/2)/2, -sqrt(13)/2 - S(3)/2)])
    assert gsolve(eqs[0], eqs[1], x, y) == set([
        (-sqrt(-S(3)/2 + sqrt(13)/2)*sqrt(S(9)/2 + 13*sqrt(13)/2)/2 - S(1)/4 +
        3*sqrt(13)/4, -S(3)/2 + sqrt(13)/2), (-3*sqrt(13)/4 - S(1)/4 + sqrt(-S(9)/2 +
        13*sqrt(13)/2)*sqrt(S(3)/2 + sqrt(13)/2)/2, -sqrt(13)/2 - S(3)/2)])
    assert gsolve(y-(x-3), x**2/9 + (y/4 - S(1)/4)**2 - 1,
        x, y, check=False) == set([(0, -3), (S(72)/25, -S(3)/25)])
    assert gsolve(
        -8*x*y - 6*x - 7*y**2 - y - 7,
        -9*x**2 - 10*x*y - x - 4*y**2 + 5*y, x, y, check=False) == set([
        (Mul(-1, CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 0) +
        7*CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 0)**2 +
        7, evaluate=False)/(8*CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 0) + 6),
        CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 0)),
        (Mul(-1, CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 1) + 7 +
        7*CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 1)**2, evaluate=False)/(6 +
        8*CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 1)),
        CRootOf(137*y**4 - 366*y**3 - 115*y**2 - 536*y + 399, 1))])
    assert gsolve(Circle(Point2D(4, -8), 2), Ellipse(Point2D(4, 5), 5, 7)) == set()
