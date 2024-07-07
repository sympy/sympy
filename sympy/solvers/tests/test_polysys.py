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
from sympy.solvers.polysys import (red, red_level, red_set, 
                                    subresultant_polynomials,
                                    subresultant_coefficients,
                                    get_nice_roots, projone, projtwo,
                                    hongproj, cylindrical_algebraic_decomposition,
                                    solve_poly_system_cad)
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


def test_red_level():
    assert red_level(x, x, 0) == x
    assert red_level(x, x, 1) == 0
    assert red_level(x, x, 100) == 0 # level > degree

    assert red_level(x**3 + x**2 + x + 1, x, 2) == x + 1

    assert red_level(y*x**3 + y**2 * x**2 + 1, x, 1) == y**2 * x**2 + 1


def test_red_set():
    assert red_set(1, x) == []
    assert red_set(x, x) == [x]

    assert red_set(x**3 + x**2 + x + 1, x) == [x**3 + x**2 + x + 1,
                                               x**2 + x + 1,
                                               x + 1,
                                               1]
    
    assert red_set(y*x + y, x) == [y*x + y, y]


def test_subresultant_polynomials():
    # simple monic univariate examples
    assert subresultant_polynomials(x**2+1, x**2-1, x) == [4, -2, -1+x**2]
    assert subresultant_polynomials(x**3+1, x**2-1, x) == [0, 1+x, -1+x**2]

    # should be order-invariant
    assert subresultant_polynomials(x**3+1, x**2-1, x) == subresultant_polynomials(x**3+1, x**2-1, x)

    # battery of univariate examples with diff degrees and coefficients
    # all checked in Mathematica
    # note scaling: 2nd example, the PRS element is 2-x+3x^3
    #   but the corresponding subresultant is 6-3x+9x^3
    #   ie times LC(g)^(deg(f)-deg(g)-1) = 3^1 = 3
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
            1349-743*x - 901*x**2,  
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
            13 + 13*x,
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

    # a multivariate example
    # note the division by y^2!
    assert subresultant_polynomials(y*x**2+1, y**2*x**2-1, x) ==\
        [y**4 + 2*y**3 + y**2, -y**2 - y, (x**2*y**2 - 1)/y**2]
    

def test_subresultant_coefficients():
    # simple monic univariate examples
    assert subresultant_coefficients(x**2+1, x**2-1, x) == [4, 0, 1]
    assert subresultant_polynomials(x**3+1, x**2-1, x) == [0, 1, 1]
    
    # should be order-invariant
    assert subresultant_coefficients(x**3+1, x**2-1, x) == subresultant_coefficients(x**3+1, x**2-1, x)

    # battery of univariate examples with diff degrees and coefficients
    # all checked in Mathematica
    # note scaling: 2nd example, the PRS element is 2-x+3x^3
    #   but the corresponding subresultant is 6-3x+9x^3
    #   ie times LC(g)^(deg(f)-deg(g)-1) = 3^1 = 3
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

    # a multivariate example
    # note the division by y^2!
    assert subresultant_polynomials(y*x**2+1, y**2*x**2-1, x) ==\
        [y**4 + 2*y**3 + y**2, 0, 1]
    








