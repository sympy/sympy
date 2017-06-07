from __future__ import print_function, division

from sympy import sqrt

from sympy.core import S

from sympy.integrals.intpoly import (is_vertex, intersection, norm,
                                     decompose, best_origin,
                                     hyperplane_parameters,
                                     integration_reduction, polytope_integrate)

from sympy.geometry.line import Segment2D
from sympy.geometry.polygon import Polygon
from sympy.geometry.point import Point
from sympy.abc import x, y


def test_is_vertex():
    pass


def test_intersection():
    pass


def test_norm():
    pass


def test_decompose():
    assert decompose(x) == {1: x}
    assert decompose(x**2) == {2: x**2}
    assert decompose(x*y) == {2: x*y}
    assert decompose(x + y) == {1: x + y}
    assert decompose(x**2 + y) == {1: y, 2: x**2}
    assert decompose(8*x**2 + 4*y) == {1: 4*y, 2: 8*x**2}
    assert decompose(x**2 + 3*y*x) == {2: x**2 + 3*x*y}
    assert decompose(9*x**2 + y + 4*x + x**3 + y**2*x) ==\
        {1: 4*x + y, 2: 9*x**2, 3: x**3 + x*y**2}


def test_best_origin():
    expr1 = y ** 2 * x ** 5 + y ** 5 * x ** 7 + 7 * x + x ** 12 + y ** 7 * x

    l1 = Segment2D(Point(0, 3), Point(1, 1))
    l2 = Segment2D(Point(1.5, 0), Point(1.5, 3))
    l3 = Segment2D(Point(0, 1.5), Point(3, 1.5))
    l4 = Segment2D(Point(0, 2), Point(2, 0))
    l5 = Segment2D(Point(0, 2), Point(1, 1))
    l6 = Segment2D(Point(2, 0), Point(1, 1))

    assert best_origin((2, 1), 3, l1, expr1) == (0, 3.0)
    assert best_origin((2, 0), 3, l2, x ** 7) == (1.5, 0)
    assert best_origin((0, 2), 3, l3, x ** 7) == (0, 1.5)
    assert best_origin((1, 1), 2, l4, x ** 7 * y ** 3) == (0, 2.0)
    assert best_origin((1, 1), 2, l4, x ** 3 * y ** 7) == (2.0, 0)
    assert best_origin((1, 1), 2, l5, x ** 2 * y ** 9) == (0, 2.0)
    assert best_origin((1, 1), 2, l6, x ** 9 * y ** 2) == (2.0, 0)


def test_hyperplane_parameters():
    pass


def test_integration_reduction():
    pass


def test_polytope_integrate():
    #  Convex 2-Polytopes
    assert polytope_integrate(Polygon(Point(0, 0), Point(0, 2),
                                      Point(4, 0)), 1, dims=(x, y)) == 4
    assert polytope_integrate(Polygon(Point(0, 0), Point(0, 1),
                                      Point(1, 1), Point(1, 0)), x * y) == 1/4
    assert polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)),
                              6*x**2 - 40*y) == float(-935)/3

    assert polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(3)),
                                      Point(sqrt(3), sqrt(3)),
                                      Point(sqrt(3), 0)), 1) == 3

    hexagon = Polygon(Point(0, 0), Point(-sqrt(3) / 2, 0.5),
                      Point(-sqrt(3) / 2, 3 / 2), Point(0, 2),
                      Point(sqrt(3) / 2, 3 / 2), Point(sqrt(3) / 2, 0.5))

    assert polytope_integrate(hexagon, 1) == S(3*sqrt(3)) / 2

    #  Non-convex polytopes
    assert polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1),
                                      Point(1, 1), Point(0, 0),
                                      Point(1, -1)), 1) == 3
    assert polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1),
                                      Point(0, 0), Point(1, 1),
                                      Point(1, -1), Point(0, 0)), 1) == 2
