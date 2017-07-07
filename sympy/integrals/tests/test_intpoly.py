from __future__ import print_function, division

from sympy import sqrt

from sympy.core import S

from sympy.integrals.intpoly import (decompose, best_origin,
                                     polytope_integrate)

from sympy.geometry.line import Segment2D
from sympy.geometry.polygon import Polygon
from sympy.geometry.point import Point
from sympy.abc import x, y

from sympy.utilities.pytest import raises, XFAIL


def test_decompose():
    assert decompose(x) == {1: x}
    assert decompose(x**2) == {2: x**2}
    assert decompose(x*y) == {2: x*y}
    assert decompose(x + y) == {1: x + y}
    assert decompose(x**2 + y) == {1: y, 2: x**2}
    assert decompose(8*x**2 + 4*y + 7) == {0: 7, 1: 4*y, 2: 8*x**2}
    assert decompose(x**2 + 3*y*x) == {2: x**2 + 3*x*y}
    assert decompose(9*x**2 + y + 4*x + x**3 + y**2*x + 3) ==\
        {0: 3, 1: 4*x + y, 2: 9*x**2, 3: x**3 + x*y**2}

    assert decompose(x, True) == [x]
    assert decompose(x ** 2, True) == [x ** 2]
    assert decompose(x * y, True) == [x * y]
    assert decompose(x + y, True) == [x, y]
    assert decompose(x ** 2 + y, True) == [y, x ** 2]
    assert decompose(8 * x ** 2 + 4 * y + 7, True) == [7, 4*y, 8*x**2]
    assert decompose(x ** 2 + 3 * y * x, True) == [x ** 2, 3 * x * y]
    assert decompose(9 * x ** 2 + y + 4 * x + x ** 3 + y ** 2 * x + 3, True) == \
           [3, y, x**3, 4*x, 9*x**2, x*y**2]


def test_best_origin():
    expr1 = y ** 2 * x ** 5 + y ** 5 * x ** 7 + 7 * x + x ** 12 + y ** 7 * x

    l1 = Segment2D(Point(0, 3), Point(1, 1))
    l2 = Segment2D(Point(S(3) / 2, 0), Point(S(3) / 2, 3))
    l3 = Segment2D(Point(0, S(3) / 2), Point(3, S(3) / 2))
    l4 = Segment2D(Point(0, 2), Point(2, 0))
    l5 = Segment2D(Point(0, 2), Point(1, 1))
    l6 = Segment2D(Point(2, 0), Point(1, 1))

    assert best_origin((2, 1), 3, l1, expr1) == (0, 3)
    assert best_origin((2, 0), 3, l2, x ** 7) == (S(3) / 2, 0)
    assert best_origin((0, 2), 3, l3, x ** 7) == (0, S(3) / 2)
    assert best_origin((1, 1), 2, l4, x ** 7 * y ** 3) == (0, 2)
    assert best_origin((1, 1), 2, l4, x ** 3 * y ** 7) == (2, 0)
    assert best_origin((1, 1), 2, l5, x ** 2 * y ** 9) == (0, 2)
    assert best_origin((1, 1), 2, l6, x ** 9 * y ** 2) == (2, 0)


def test_polytope_integrate():
    #  Convex 2-Polytopes
    #  Vertex representation
    assert polytope_integrate(Polygon(Point(0, 0), Point(0, 2),
                                      Point(4, 0)), 1, dims=(x, y)) == 4
    assert polytope_integrate(Polygon(Point(0, 0), Point(0, 1),
                                      Point(1, 1), Point(1, 0)), x * y) ==\
                                      S(1)/4
    assert polytope_integrate(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)),
                              6*x**2 - 40*y) == S(-935)/3

    assert polytope_integrate(Polygon(Point(0, 0), Point(0, sqrt(3)),
                                      Point(sqrt(3), sqrt(3)),
                                      Point(sqrt(3), 0)), 1) == 3

    hexagon = Polygon(Point(0, 0), Point(-sqrt(3) / 2, S(1)/2),
                      Point(-sqrt(3) / 2, 3 / 2), Point(0, 2),
                      Point(sqrt(3) / 2, 3 / 2), Point(sqrt(3) / 2, S(1)/2))

    assert polytope_integrate(hexagon, 1) == S(3*sqrt(3)) / 2

    #  Hyperplane representation
    assert polytope_integrate([((-1, 0), 0), ((1, 2), 4),
                               ((0, -1), 0)], 1, dims=(x, y)) == 4
    assert polytope_integrate([((-1, 0), 0), ((0, 1), 1),
                               ((1, 0), 1), ((0, -1), 0)], x * y) == S(1)/4
    assert polytope_integrate([((0, 1), 3), ((1, -2), -1),
                               ((-2, -1), -3)], 6*x**2 - 40*y) == S(-935)/3
    assert polytope_integrate([((-1, 0), 0), ((0, sqrt(3)), 3),
                               ((sqrt(3), 0), 3), ((0, -1), 0)], 1) == 3

    hexagon = [((-1 / 2, -sqrt(3) / 2), 0),
               ((-1, 0), sqrt(3) / 2),
               ((-1 / 2, sqrt(3) / 2), sqrt(3)),
               ((1 / 2, sqrt(3) / 2), sqrt(3)),
               ((1, 0), sqrt(3) / 2),
               ((1 / 2, -sqrt(3) / 2), 0)]
    assert polytope_integrate(hexagon, 1) == S(3*sqrt(3)) / 2

    #  Non-convex polytopes
    #  Vertex representation
    assert polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1),
                                      Point(1, 1), Point(0, 0),
                                      Point(1, -1)), 1) == 3
    assert polytope_integrate(Polygon(Point(-1, -1), Point(-1, 1),
                                      Point(0, 0), Point(1, 1),
                                      Point(1, -1), Point(0, 0)), 1) == 2
    #  Hyperplane representation
    assert polytope_integrate([((-1, 0), 1), ((0, 1), 1), ((1, -1), 0),
                               ((1, 1), 0), ((0, -1), 1)], 1) == 3
    assert polytope_integrate([((-1, 0), 1), ((1, 1), 0), ((-1, 1), 0),
                               ((1, 0), 1), ((-1, -1), 0),
                               ((1, -1), 0)], 1) == 2

    #  Tests for 2D polytopes mentioned in Chin et al(Page 10):
    #  http://dilbert.engr.ucdavis.edu/~suku/quadrature/cls-integration.pdf
    fig1 = Polygon(Point(1.220, -0.827), Point(-1.490, -4.503),
                   Point(-3.766, -1.622), Point(-4.240, -0.091),
                   Point(-3.160, 4), Point(-0.981, 4.447),
                   Point(0.132, 4.027))
    assert polytope_integrate(fig1, x**2 + x*y + y**2) ==\
        S(2031627344735367)/(8*10**12)

    fig2 = Polygon(Point(4.561, 2.317), Point(1.491, -1.315),
                   Point(-3.310, -3.164), Point(-4.845, -3.110),
                   Point(-4.569, 1.867))
    assert polytope_integrate(fig2, x**2 + x*y + y**2) ==\
        S(517091313866043)/(16*10**11)

    fig3 = Polygon(Point(-2.740, -1.888), Point(-3.292, 4.233),
                   Point(-2.723, -0.697), Point(-0.643, -3.151))
    assert polytope_integrate(fig3, x**2 + x*y + y**2) ==\
        S(147449361647041)/(8*10**12)

    fig4 = Polygon(Point(0.211, -4.622), Point(-2.684, 3.851),
                   Point(0.468, 4.879), Point(4.630, -1.325),
                   Point(-0.411, -1.044))
    assert polytope_integrate(fig4, x**2 + x*y + y**2) ==\
        S(180742845225803)/(10**12)

    #  Tests for many polynomials with maximum degree given.
    tri = Polygon(Point(0, 3), Point(5, 3), Point(1, 1))
    polys = []
    expr1 = x**9*y + x**7*y**3 + 2*x**2*y**8
    expr2 = x**6*y**4 + x**5*y**5 + 2*y**10
    expr3 = x**10 + x**9*y + x**8*y**2 + x**5*y**5
    polys.extend((expr1, expr2, expr3))
    result_dict = polytope_integrate(tri, polys, max_degree=10)
    assert result_dict[expr1] == 615780107/594
    assert result_dict[expr2] == 13062161/27
    assert result_dict[expr3] == 1946257153/924


@XFAIL
def test_polytopes_intersecting_sides():
    #  Intersecting polygons not implemented yet in SymPy. Will be implemented
    #  soon. As of now, the intersection point will have to be manually
    #  supplied by user.
    fig5 = Polygon(Point(-4.165, -0.832), Point(-3.668, 1.568),
                   Point(-3.266, 1.279), Point(-1.090, -2.080),
                   Point(3.313, -0.683), Point(3.033, -4.845),
                   Point(-4.395, 4.840), Point(-1.007, -3.328))
    assert polytope_integrate(fig5, x**2 + x*y + y**2) ==\
        S(1633405224899363)/(24*10**12)

    fig6 = Polygon(Point(-3.018, -4.473), Point(-0.103, 2.378),
                   Point(-1.605, -2.308), Point(4.516, -0.771),
                   Point(4.203, 0.478))
    assert polytope_integrate(fig6, x**2 + x*y + y**2) ==\
        S(88161333955921)/(3*10**12)
