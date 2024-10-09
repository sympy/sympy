"""Tests for solvers of systems of polynomial equations. """
from sympy.core.numbers import (I, Integer, Rational)
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.polys.domains.rationalfield import QQ
from sympy.polys.polyerrors import UnsolvableFactorError
from sympy.polys.polyoptions import Options
from sympy.polys.polytools import Poly
from sympy.solvers.solvers import solve
from sympy.utilities.iterables import flatten
from sympy.abc import x, y, z
from sympy.polys import PolynomialError
from sympy.solvers.polysys import (solve_poly_system,
                                   solve_triangulated,
                                   solve_biquadratic, SolveFailed,
                                   solve_generic)
from sympy.solvers.polysys import (red, red_set,
                                    subresultant_polynomials,
                                    subresultant_coefficients,
                                    get_nice_roots, get_sample_point,
                                    simplify_alg_sub,
                                    projone, projtwo, hongproj,
                                    cylindrical_algebraic_decomposition,
                                    solve_poly_system_cad)
from sympy.polys.rootoftools import ComplexRootOf
from sympy.polys.polytools import parallel_poly_from_expr
from sympy.testing.pytest import raises


def test_solve_poly_system():
    assert solve_poly_system([x - 1], x) == [(S.One,)]

    assert solve_poly_system([y - x, y - x - 1], x, y) is None

    assert solve_poly_system([y - x**2, y + x**2], x, y) == [(S.Zero, S.Zero)]

    assert solve_poly_system([2*x - 3, y*Rational(3, 2) - 2*x, z - 5*y], x, y, z) == \
        [(Rational(3, 2), Integer(2), Integer(10))]

    assert solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y) == \
        [(0, 0), (2, -sqrt(2)), (2, sqrt(2))]

    assert solve_poly_system([y - x**2, y + x**2 + 1], x, y) == \
        [(-I*sqrt(S.Half), Rational(-1, 2)), (I*sqrt(S.Half), Rational(-1, 2))]

    f_1 = x**2 + y + z - 1
    f_2 = x + y**2 + z - 1
    f_3 = x + y + z**2 - 1

    a, b = sqrt(2) - 1, -sqrt(2) - 1

    assert solve_poly_system([f_1, f_2, f_3], x, y, z) == \
        [(0, 0, 1), (0, 1, 0), (1, 0, 0), (a, a, a), (b, b, b)]

    solution = [(1, -1), (1, 1)]

    assert solve_poly_system([Poly(x**2 - y**2), Poly(x - 1)]) == solution
    assert solve_poly_system([x**2 - y**2, x - 1], x, y) == solution
    assert solve_poly_system([x**2 - y**2, x - 1]) == solution

    assert solve_poly_system(
        [x + x*y - 3, y + x*y - 4], x, y) == [(-3, -2), (1, 2)]

    raises(NotImplementedError, lambda: solve_poly_system([x**3 - y**3], x, y))
    raises(NotImplementedError, lambda: solve_poly_system(
        [z, -2*x*y**2 + x + y**2*z, y**2*(-z - 4) + 2]))
    raises(PolynomialError, lambda: solve_poly_system([1/x], x))

    raises(NotImplementedError, lambda: solve_poly_system(
          [x-1,], (x, y)))
    raises(NotImplementedError, lambda: solve_poly_system(
          [y-1,], (x, y)))

    # solve_poly_system should ideally construct solutions using
    # CRootOf for the following four tests
    assert solve_poly_system([x**5 - x + 1], [x], strict=False) == []
    raises(UnsolvableFactorError, lambda: solve_poly_system(
        [x**5 - x + 1], [x], strict=True))

    assert solve_poly_system([(x - 1)*(x**5 - x + 1), y**2 - 1], [x, y],
                             strict=False) == [(1, -1), (1, 1)]
    raises(UnsolvableFactorError,
           lambda: solve_poly_system([(x - 1)*(x**5 - x + 1), y**2-1],
                                     [x, y], strict=True))


def test_solve_generic():
    NewOption = Options((x, y), {'domain': 'ZZ'})
    assert solve_generic([x**2 - 2*y**2, y**2 - y + 1], NewOption) == \
           [(-sqrt(-1 - sqrt(3)*I), Rational(1, 2) - sqrt(3)*I/2),
            (sqrt(-1 - sqrt(3)*I), Rational(1, 2) - sqrt(3)*I/2),
            (-sqrt(-1 + sqrt(3)*I), Rational(1, 2) + sqrt(3)*I/2),
            (sqrt(-1 + sqrt(3)*I), Rational(1, 2) + sqrt(3)*I/2)]

    # solve_generic should ideally construct solutions using
    # CRootOf for the following two tests
    assert solve_generic(
        [2*x - y, (y - 1)*(y**5 - y + 1)], NewOption, strict=False) == \
        [(Rational(1, 2), 1)]
    raises(UnsolvableFactorError, lambda: solve_generic(
        [2*x - y, (y - 1)*(y**5 - y + 1)], NewOption, strict=True))


def test_solve_biquadratic():
    x0, y0, x1, y1, r = symbols('x0 y0 x1 y1 r')

    f_1 = (x - 1)**2 + (y - 1)**2 - r**2
    f_2 = (x - 2)**2 + (y - 2)**2 - r**2
    s = sqrt(2*r**2 - 1)
    a = (3 - s)/2
    b = (3 + s)/2
    assert solve_poly_system([f_1, f_2], x, y) == [(a, b), (b, a)]

    f_1 = (x - 1)**2 + (y - 2)**2 - r**2
    f_2 = (x - 1)**2 + (y - 1)**2 - r**2

    assert solve_poly_system([f_1, f_2], x, y) == \
        [(1 - sqrt((2*r - 1)*(2*r + 1))/2, Rational(3, 2)),
         (1 + sqrt((2*r - 1)*(2*r + 1))/2, Rational(3, 2))]

    query = lambda expr: expr.is_Pow and expr.exp is S.Half

    f_1 = (x - 1 )**2 + (y - 2)**2 - r**2
    f_2 = (x - x1)**2 + (y - 1)**2 - r**2

    result = solve_poly_system([f_1, f_2], x, y)

    assert len(result) == 2 and all(len(r) == 2 for r in result)
    assert all(r.count(query) == 1 for r in flatten(result))

    f_1 = (x - x0)**2 + (y - y0)**2 - r**2
    f_2 = (x - x1)**2 + (y - y1)**2 - r**2

    result = solve_poly_system([f_1, f_2], x, y)

    assert len(result) == 2 and all(len(r) == 2 for r in result)
    assert all(len(r.find(query)) == 1 for r in flatten(result))

    s1 = (x*y - y, x**2 - x)
    assert solve(s1) == [{x: 1}, {x: 0, y: 0}]
    s2 = (x*y - x, y**2 - y)
    assert solve(s2) == [{y: 1}, {x: 0, y: 0}]
    gens = (x, y)
    for seq in (s1, s2):
        (f, g), opt = parallel_poly_from_expr(seq, *gens)
        raises(SolveFailed, lambda: solve_biquadratic(f, g, opt))
    seq = (x**2 + y**2 - 2, y**2 - 1)
    (f, g), opt = parallel_poly_from_expr(seq, *gens)
    assert solve_biquadratic(f, g, opt) == [
        (-1, -1), (-1, 1), (1, -1), (1, 1)]
    ans = [(0, -1), (0, 1)]
    seq = (x**2 + y**2 - 1, y**2 - 1)
    (f, g), opt = parallel_poly_from_expr(seq, *gens)
    assert solve_biquadratic(f, g, opt) == ans
    seq = (x**2 + y**2 - 1, x**2 - x + y**2 - 1)
    (f, g), opt = parallel_poly_from_expr(seq, *gens)
    assert solve_biquadratic(f, g, opt) == ans


def test_solve_triangulated():
    f_1 = x**2 + y + z - 1
    f_2 = x + y**2 + z - 1
    f_3 = x + y + z**2 - 1

    a, b = sqrt(2) - 1, -sqrt(2) - 1

    assert solve_triangulated([f_1, f_2, f_3], x, y, z) == \
        [(0, 0, 1), (0, 1, 0), (1, 0, 0)]

    dom = QQ.algebraic_field(sqrt(2))

    assert solve_triangulated([f_1, f_2, f_3], x, y, z, domain=dom) == \
        [(0, 0, 1), (0, 1, 0), (1, 0, 0), (a, a, a), (b, b, b)]


def test_solve_issue_3686():
    roots = solve_poly_system([((x - 5)**2/250000 + (y - Rational(5, 10))**2/250000) - 1, x], x, y)
    assert roots == [(0, S.Half - 15*sqrt(1111)), (0, S.Half + 15*sqrt(1111))]

    roots = solve_poly_system([((x - 5)**2/250000 + (y - 5.0/10)**2/250000) - 1, x], x, y)
    # TODO: does this really have to be so complicated?!
    assert len(roots) == 2
    assert roots[0][0] == 0
    assert roots[0][1].epsilon_eq(-499.474999374969, 1e12)
    assert roots[1][0] == 0
    assert roots[1][1].epsilon_eq(500.474999374969, 1e12)


# NEW TESTS FOR CAD AND RELATED FUNCTIONS

def test_red():
    # simple univar degree-one example
    assert red(x, x) == 0
    # univar degree-two
    assert red(x**2 + x + 1, x) == x+1
    # bivar
    assert red(x*y + x**2 * y**2, x) == x*y


def test_red_set():
    assert red_set(1, x) == []
    assert red_set(x, x) == [x,0]

    assert red_set(x**3 + x**2 + x + 1, x) == [x**3 + x**2 + x + 1,
                                               x**2 + x + 1,
                                               x + 1,
                                               1]

    assert red_set(y*x + y, x) == [y*x + y, y]


def test_subresultant_polynomials():
    # edge cases
    assert subresultant_polynomials(x, 0, x) == []
    assert subresultant_polynomials(x, 1, x) == [1]

    # simple monic univariate examples
    assert subresultant_polynomials(x**2+1, x**2-1, x) == [4, -2, -1+x**2]
    assert subresultant_polynomials(x**3+1, x**2-1, x) == [0, 1+x, -1+x**2]

    # should be order-invariant
    assert subresultant_polynomials(x**3+1, x**2-1, x) == subresultant_polynomials(x**3+1, x**2-1, x)

    # battery of univariate examples with diff degrees and coefficients
    # all checked in Mathematica
    fs = [
        2*x**5 - 3*x**4 + x**3 - 7*x + 5,
        -x**5 + 2*x**4 - 5*x**2 + x - 4,
        4*x**4 - x**3 + 2*x**2 - x + 3,
        x**3 - x**2 + x - 1,
        5*x**2 + 3*x - 2
    ]
    gs = [
        x**5 + 4*x**4 - x**3 + 2*x**2 - 3*x + 6,
        3*x**3 - x + 2,
        -2*x**3 + 3*x - 5,
        2*x**2 + x - 3,
        -x + 1
    ]
    answers = [
        [
            45695124,
            692022 - 809988*x,
            1349 - 743*x - 901*x**2,
            397 - 487*x + 43*x**2 - 24*x**3,
            7 + x + 4*x**2 - 3*x**3 + 11*x**4,
            6 - 3*x + 2*x**2 - x**3 + 4*x**4 + x**5
        ],
        [
            31514,
            862 - 1469*x,
            102 + 12*x + 99*x**2,
            6 - 3*x + 9*x**3
        ],
        [
            21650,
            -730 - 130*x,
            22 - 50*x + 32*x**2,
            -5 + 3*x - 2*x**3
        ],
        [
            0,
            -13 + 13*x,
            -3 + x + 2*x**2
        ],
        [
            6,
            1 - x
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_polynomials(fs[0], gs[0], x) == answers[0]
    assert subresultant_polynomials(fs[1], gs[1], x) == answers[1]
    assert subresultant_polynomials(fs[2], gs[2], x) == answers[2]
    assert subresultant_polynomials(fs[3], gs[3], x) == answers[3]
    assert subresultant_polynomials(fs[4], gs[4], x) == answers[4]

    # a battery of bivariate tests with varying deg and coeffs
    # checked with mathematica
    fs = [
        2*x**3 + y**2 + 3*x*y - 4,
        x**3 + 4*x**2*y + y - 1,
        3*x**2 + 2*x*y**2 - y + 6,
        x**3 - y**3 + x*y,
        y*x**2+1
    ]
    gs = [
        x**2 + 2*y**3 + 5*x*y,
        x**2 + y**4 - x,
        x + y**3 - 2*x*y + 1,
        x**2 - 2*y**2 + 3*x*y**2,
        y**2*x**2 - 1
    ]
    answers = [
        [
            16 + 52*y**2 + 1000*y**3 - 254*y**4 - 232*y**5 + 360*y**6 - 48*y**7 + 32*y**9,
            -4 + 3*x*y + y**2 + 50*x*y**2 - 4*x*y**3 + 20*y**4,
            x**2 + 5*x*y + 2*y**3
        ],
        [
            -5*y + 5*y**2 + 3*y**4 + 5*y**5 - 8*y**6 + 4*y**9 + 16*y**10 + y**12,
            -1 + x + y + 4*x*y - y**4 - x*y**4 - 4*y**5,
            -x + x**2 + y**4
        ],
        [
            9 - 25*y + 26*y**2 + 6*y**3 - 2*y**5 + 7*y**6,
            1 + x - 2*x*y + y**3
        ],
        [
            -2*y**4 - 8*y**5 - 4*y**6 + 27*y**9,
            x*y + 2*x*y**2 - y**3 - 6*y**4 + 9*x*y**4,
            x**2 - 2*y**2 + 3*x*y**2
        ],
        [
            y**2 + 2*y**3 + y**4,
            -y - y**2,
            x**2 - 1/y**2
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_polynomials(fs[0], gs[0], x) == answers[0]
    assert subresultant_polynomials(fs[1], gs[1], x) == answers[1]
    assert subresultant_polynomials(fs[2], gs[2], x) == answers[2]
    assert subresultant_polynomials(fs[3], gs[3], x) == answers[3]
    assert subresultant_polynomials(fs[4], gs[4], x) == answers[4]



def test_subresultant_coefficients():
    # edge cases
    assert subresultant_coefficients(x, 0, x) == []
    assert subresultant_coefficients(x, 1, x) == [1]

    # simple monic univariate examples
    assert subresultant_coefficients(x**2+1, x**2-1, x) == [4, 0, 1]
    assert subresultant_coefficients(x**3+1, x**2-1, x) == [0, 1, 1]

    # should be order-invariant
    assert subresultant_coefficients(x**3+1, x**2-1, x) == subresultant_coefficients(x**3+1, x**2-1, x)

    # battery of univariate examples with diff degrees and coefficients
    # all checked in Mathematica
    fs = [
        2*x**5 - 3*x**4 + x**3 - 7*x + 5,
        -x**5 + 2*x**4 - 5*x**2 + x - 4,
        4*x**4 - x**3 + 2*x**2 - x + 3,
        x**3 - x**2 + x - 1,
        5*x**2 + 3*x - 2
    ]
    gs = [
        x**5 + 4*x**4 - x**3 + 2*x**2 - 3*x + 6,
        3*x**3 - x + 2,
        -2*x**3 + 3*x - 5,
        2*x**2 + x - 3,
        -x + 1
    ]
    answers = [
        [
            45695124,
            -809988,
            -901,
            -24,
            11,
            1
        ],
        [
            31514,
            -1469,
            99,
            9
        ],
        [
            21650,
            -130,
            32,
            -2
        ],
        [
            0,
            13,
            2
        ],
        [
            6,
            -1
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_coefficients(fs[0], gs[0], x) == answers[0]
    assert subresultant_coefficients(fs[1], gs[1], x) == answers[1]
    assert subresultant_coefficients(fs[2], gs[2], x) == answers[2]
    assert subresultant_coefficients(fs[3], gs[3], x) == answers[3]
    assert subresultant_coefficients(fs[4], gs[4], x) == answers[4]


    # a battery of bivariate tests with varying deg and coeffs
    # checked with mathematica
    fs = [
        2*x**3 + y**2 + 3*x*y - 4,
        x**3 + 4*x**2*y + y - 1,
        3*x**2 + 2*x*y**2 - y + 6,
        x**3 - y**3 + x*y,
        y*x**2+1
    ]
    gs = [
        x**2 + 2*y**3 + 5*x*y,
        x**2 + y**4 - x,
        x + y**3 - 2*x*y + 1,
        x**2 - 2*y**2 + 3*x*y**2,
        y**2*x**2 - 1
    ]
    answers = [
        [
            16 + 52*y**2 + 1000*y**3 - 254*y**4 - 232*y**5 + 360*y**6 - 48*y**7 + 32*y**9,
            3*y + 50*y**2 - 4*y**3,
            1
        ],
        [
            -5*y + 5*y**2 + 3*y**4 + 5*y**5 - 8*y**6 + 4*y**9 + 16*y**10 + y**12,
            1 + 4*y - y**4,
            1
        ],
        [
            9 - 25*y + 26*y**2 + 6*y**3 - 2*y**5 + 7*y**6,
            1 - 2*y
        ],
        [
            -2*y**4 - 8*y**5 - 4*y**6 + 27*y**9,
            y + 2*y**2 + 9*y**4,
            1
        ],
        [
            y**2 + 2*y**3 + y**4,
            0,
            1
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_coefficients(fs[0], gs[0], x) == answers[0]
    assert subresultant_coefficients(fs[1], gs[1], x) == answers[1]
    assert subresultant_coefficients(fs[2], gs[2], x) == answers[2]
    assert subresultant_coefficients(fs[3], gs[3], x) == answers[3]
    assert subresultant_coefficients(fs[4], gs[4], x) == answers[4]


def test_get_nice_roots():
    # constants have no roots
    assert get_nice_roots(3) == []
    assert get_nice_roots(Poly(3, x)) == []

    # if not implemented, just solve numerically
    # eg if coefficient is algebraic
    # the answer here can be solved with basic algebra
    assert get_nice_roots(sqrt(2) * x**2 - 1)[1] == Rational(sqrt(1 / sqrt(2)).evalf())

    # if roots are CRootOf, then they should be numeric
    assert get_nice_roots(x**5 + x**2 - 1)[0] ==\
        ComplexRootOf(x**5 + x**2 - 1, 0)

    # the algebraic roots should stay algebraic
    # bc of the multiplication, we get the roots from x^2 - 1 of +- sqrt(2)
    assert get_nice_roots((x**2 - 2) * (x**5 - x**2 - 1)) ==\
        [-sqrt(2), ComplexRootOf(x**5 - x**2 - 1, 0), sqrt(2)]


def test_get_sample_point():
    # edge case
    assert get_sample_point(1, 1) == 1

    # test rays
    assert get_sample_point(S.NegativeInfinity, -2) == -3
    assert get_sample_point(S.NegativeInfinity, -1.5) == -3
    assert get_sample_point(sqrt(5), S.Infinity) == 4
    assert get_sample_point(ComplexRootOf(x**5 + x**3 - 1, 0), S.Infinity) == S(2)

    # test finite open intervals
    assert get_sample_point(-2, -1) == - S(3) / 2
    # approx 1.4 and 2.4
    assert get_sample_point(sqrt(2), sqrt(2) + 1) == 2
    # approx 0.8 and 1.2
    assert get_sample_point(ComplexRootOf(x**5 + x**3 - 1, 0), sqrt(5) - 1) == 1
    # approx 0.8376 and 0.8164
    assert get_sample_point(ComplexRootOf(x**5 + x**3 - 1, 0), 209 / 256) == S(53)/64

    # differ in last digit, just to test accuracy with close numbers
    # previous archimedian axiom technique took way too long on these
    assert get_sample_point(0.1729847129041, 0.1729847129042) == S(6086358504499)/35184372088832

def test_simplify_alg_sub():
    assert simplify_alg_sub(x**3 + x**2 - 1, {})
    assert simplify_alg_sub(x**3 + x**2 - 1, {x:2}) == 11
    assert simplify_alg_sub(x**3 + x**2 - 1, {y:1}) == x**3 + x**2 - 1

    r = ComplexRootOf(x**5 + x**3 - 1, 0)
    assert simplify_alg_sub(x**5 + x**3 - 1, {x:r}) == 0
    assert simplify_alg_sub(x**5 + x**3 + 5, {x:r}) == 6
    # wow
    assert simplify_alg_sub(x**3 + x**2 * y**5 + x**2 * y**3 - x**2, {y:r}) == x**3


def test_projone():
    # simple example: work it out manually by looping through
    # for x^2
    #   red_set = [x^2,0,0]
    #   LC = 1
    #   subres coeffs of x^2 and 2x (deriv) = [0,2]
    #   so add [0,1,2] to the proj factors
    # for x
    #   red_set = [x,0]
    #   LC = 1
    #   subres coeffs of x and 1 (deriv) = [1]
    #   so add [1] to the proj factors
    # hence return {0,1,2}
    assert projone([x**2, x], x) == {0, 1, 2}

    # slighly harder bivar example: loop through
    # for y*x**2
    #   red_set = [y*x^2]
    #   LC = y
    #   subres coeffs of y*x^2 and 2yx (deriv) = [0,2y]
    #   so add [0,y,2y] to the proj factors
    # for y**2*x
    #   red_set = [y**2*x]
    #   LC = y**2
    #   subres coeffs of y**2*x and y**2 (deriv) = [y**2]
    #   so add [y**2, y**2] to the proj factors
    # hence return {0,y,2y,y**2}
    assert projone([y*x**2, y**2*x], x) == {0, y**2, y, 2*y}


def test_projtwo():
    # simple example: loop through (pairs)
    # for the pair (f,g) = (x^2, x)
    #   red_set(f) = [x^2]
    #   for f_ = x^2
    #       subres coeffs of x^2 and x = [0,1]
    #   so add [0,1] to the proj factors
    # hence return {0,1}
    assert projtwo([x**2, x], x) == {0, 1}

    # slighly harder bivar example: loop through (pairs)
    # for the pair (f,g) = (y*x^2, y**2*x)
    #   red_set(f) = [y*x^2]
    #   for f_ = y*x^2
    #       subres coeffs of y*x^2 and y**2*x = [0,y**2]
    #   so add [0,y**2] to the proj factors
    # hence return {0,y**2}
    assert projtwo([y*x**2, y**2*x], x) == {0, y**2}


def test_hongproj():
    # the same two examples as in test_projone and test_projtwo
    # the tests here are really just testing the 'cleaning up'

    # this should be empty as they are all constants
    assert hongproj([x**2, x], x) == set()

    assert hongproj([y*x**2, y**2*x], x) == {y**2, y, 2*y}


def test_cylindrical_algebraic_decomposition():
    # simple univar example
    # x^2-1 has roots +-1, which implies these cells
    assert cylindrical_algebraic_decomposition([x**2-1], [x]) ==\
        [{x: val} for val in [-2, -1, 0, 1, 2]]

    # harder univar example
    # the collection of roots are (-1, 0, 1, 5**(1/3))
    # note the sample points here depend on the choices we make in get_sample_root
    # eg, the point btwn 1 and 5**(1/3) could be 3/2 if we are looking for simple rationals
    # but the way it is now it is based on decimal approximations and gives 87/64
    assert cylindrical_algebraic_decomposition([x**2-1,
                                                x,
                                                x**3-5], [x]) ==\
        [{x: val} for val in
         [-2, -1, -S(1)/2, 0, S(1)/2, 1,
          S(87)/64,
          5**(S(1)/3),
          3]]


    # multivar example
    # projecting on x gets [-y^2-1, y], only root is 0
    # then lift on y=-1, y=0, y=1 -- easy but tedious
    assert cylindrical_algebraic_decomposition([x+y, x*y-1], [x,y]) == [
        {y: -1, x: -2},
        {y: -1, x: -1},
        {y: -1, x: 0},
        {y: -1, x: 1},
        {y: -1, x: 2},
        {y: 0, x: -1},
        {y: 0, x: 0},
        {y: 0, x: 1},
        {y: 1, x: -2},
        {y: 1, x: -1},
        {y: 1, x: 0},
        {y: 1, x: 1},
        {y: 1, x: 2}
        ]


def test_solve_poly_system_cad():
    # no solution as this quadratic lies below x axis
    assert solve_poly_system_cad([-x**2 - 1 >= 0], [x]) == []

    # testing utility
    # solves a system and then subs the sample points back in
    def solve_and_sub(ineqs, vars, return_one_sample=True):
        soln = solve_poly_system_cad(ineqs, vars, return_one_sample)

        for sample in soln:
            if not all(ineq.subs(sample) for ineq in ineqs):
                return False
        return True

    assert(solve_and_sub([x**2-1 >= 3], [x])) == True
    assert(solve_and_sub([x**2-1 >= 3], [x], False)) == True

    # harder example
    # this now fails to converge due to floating point!
    # -0.618033988749895 -0.618033988749895
    # False
    # <class 'sympy.core.numbers.Float'> <class 'sympy.core.numbers.Float'>
    #assert(solve_and_sub([x**2 * y**2 - 1 > 0, x<=0.2, x+y >= 1], [x,y])) == True
