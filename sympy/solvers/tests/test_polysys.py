"""Tests for solvers of systems of polynomial equations. """
from sympy.polys.domains import  ZZ, QQ_I
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
from sympy.abc import a, b, c, x, y, z
from sympy.polys import PolynomialError
from sympy.solvers.polysys import (solve_poly_system,
                                   solve_triangulated,
                                   solve_biquadratic, SolveFailed,
                                   solve_generic, factor_system_bool,
                                   factor_system_cond, factor_system_poly,
                                   factor_system)
from sympy.polys.polytools import parallel_poly_from_expr
from sympy.testing.pytest import raises
from sympy.core.relational import Eq
from sympy.functions.elementary.trigonometric import sin, cos

from sympy.functions.elementary.exponential import exp


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


def test_factor_system():

    assert factor_system([x**2 + 2*x + 1]) ==  ([[x + 1]], True)
    assert factor_system([x**2 + 2*x + 1, y**2 + 2*y + 1]) ==  ([[x + 1, y + 1]], True)
    assert factor_system([x**2 + 1]) ==  ([[x**2 + 1]], True)
    assert factor_system([]) == ([[]], True)

    result = factor_system([x**2 + y**2 + 2*x*y, x**2 - 2], extension=sqrt(2))
    expected = ([[x + y, x - sqrt(2)], [x + y, x + sqrt(2)]], True)
    assert result == expected

    result = factor_system([x ** 2 + 1, y ** 2 + 1], gaussian=True)
    expcted = ([[x - I, y - I], [x - I, y + I], [x + I, y - I], [x + I, y + I]], True)
    assert result == expcted

    result = factor_system([x ** 2 + 1, y ** 2 + 1], domain=QQ_I)
    expected = ([[x - I, y - I], [x - I, y + I], [x + I, y - I], [x + I, y + I]], True)
    assert result == expected

    assert factor_system([0]) == ([[]], True)
    assert factor_system([1]) == ([], True)
    assert factor_system([0 , x]) == ([[x]], True)
    assert factor_system([1, 0, x]) == ([], True)

    result = factor_system([x**4 - 1, y**6 - 1])
    expected = (
        [
            [x - 1, y - 1],
            [x - 1, y**2 - y + 1],
            [x - 1, y + 1],
            [x - 1, y**2 + y + 1],
            [x + 1, y - 1],
            [x + 1, y**2 - y + 1],
            [x + 1, y + 1],
            [x + 1, y**2 + y + 1],
            [x**2 + 1, y - 1],
            [x**2 + 1, y**2 - y + 1],
            [x**2 + 1, y + 1],
            [x**2 + 1, y**2 + y + 1]
        ],
        True
    )
    assert result == expected

    result = factor_system([(x - 1)*(y - 2), (y  - 2)*(z - 3)])
    expected = ([[x - 1, z - 3], [y - 2]], True)
    assert result == expected

    result = factor_system([sin(x) ** 2 + cos(x) ** 2 - 1, x])
    expected = ([[sin(x)**2 + cos(x)**2 - 1, x]], True)
    assert result == expected

    result = factor_system([sin(x)**2 + cos(x)**2 - 1])
    expected = ([[sin(x)**2 + cos(x)**2 - 1]], True)
    assert result == expected

    result = factor_system([sin(x) ** 2 + cos(x) ** 2])
    expected = ([[sin(x)**2 + cos(x)**2]], True)
    assert result == expected

    result = factor_system([a*x*(x-1), b*y, c], [x, y])
    expected = ([[x, y], [x - 1, y]], Eq(c, 0))
    assert result == expected

    result = factor_system([x**2 - 2], [y])
    expected = ([[]], Eq(x**2 - 2, 0))
    assert result == expected

    # Generators not symbols
    expr1 = cos(x) ** 2 - sin(x) ** 2
    expr2 = cos(x) ** 2 + sin(x) ** 2 - 1
    result = factor_system([expr1, expr2])
    expected = (
        [
            [sin(x) + cos(x), sin(x)**2 + cos(x)**2 - 1],
            [-sin(x) + cos(x), sin(x)**2 + cos(x)**2 - 1]],
        True
    )
    assert result == expected

    expr3 = (cos(x) + sin(x)) ** 2 - 1
    expr4 = cos(x) ** 2 - sin(x) ** 2 - cos(2 * x)
    result = factor_system([expr3, expr4])
    expected = (
        [
            [sin(x) + cos(x) + 1, sin(x)**2 - cos(x)**2 + cos(2*x)],
            [sin(x) + cos(x) - 1, sin(x)**2 - cos(x)**2 + cos(2*x)]],
        True
    )
    assert result == expected

    expr5 = (cos(x) + sin(x)) * exp(y) - 1
    expr6 = (cos(x) - sin(x)) * exp(y) - 1
    result = factor_system([expr5, expr6])
    expected = ([[exp(y)*sin(x) + exp(y)*cos(x) - 1, -exp(y)*sin(x) + exp(y)*cos(x) - 1]], True)
    assert result == expected


def test_factor_system_poly():

    p1 = Poly(x ** 2 - 1, x)
    p2 = Poly(x ** 2 - 4, x)
    result = factor_system_poly([p1, p2])
    expected = ([[(Poly(x + 1, x, domain='ZZ'), ()), (Poly(x - 2, x, domain='ZZ'), ())],
                 [(Poly(x + 1, x, domain='ZZ'), ()), (Poly(x + 2, x, domain='ZZ'), ())],
                 [(Poly(x - 1, x, domain='ZZ'), ()), (Poly(x - 2, x, domain='ZZ'), ())],
                 [(Poly(x - 1, x, domain='ZZ'), ()), (Poly(x + 2, x, domain='ZZ'), ())]], [])
    assert result == expected

    p = Poly(x ** 2 - 1, x)
    result = factor_system_poly([p])
    expected = ([[(Poly(x - 1, x, domain='ZZ'), ())], [(Poly(x + 1, x, domain='ZZ'), ())]], [])
    assert result == expected

    p1 = Poly(x ** 2 * y - y, x, y, z)
    p2 = Poly(x ** 2 * z - z, x, y, z)
    result = factor_system_poly([p1, p2])
    expected = ([[(Poly(x + 1, x, y, z, domain='ZZ'), ())],
                 [(Poly(x - 1, x, y, z, domain='ZZ'), ())],
                 [(Poly(y, x, y, z, domain='ZZ'), ()), (Poly(z, x, y, z, domain='ZZ'), ())]], [])
    assert result == expected

    p1 = Poly(x ** 2 * (x - 1) ** 2, x)
    p2 = Poly(x * (x - 1), x)
    result = factor_system_poly([p1, p2])
    expected = ([[(Poly(x, x, domain='ZZ'), ())], [(Poly(x - 1, x, domain='ZZ'), ())]], [])
    assert result == expected

    p1 = Poly(x ** 2 + y * x, x, y, z)
    p2 = Poly(x ** 2 + z * x, x, y, z)
    result = factor_system_poly([p1, p2])
    expected = ([[(Poly(x + y, x, y, z, domain='ZZ'), ()), (Poly(x + z, x, y, z, domain='ZZ'), ())],
                 [(Poly(x, x, y, z, domain='ZZ'), ())]], [])
    assert result == expected

    p1 = Poly((a - 1) * (x - 2), x, domain=ZZ[a, b])
    p2 = Poly((b - 3) * (x - 2), x, domain=ZZ[a, b])
    result = factor_system_poly([p1, p2])
    expected = (
        [[
            (
                Poly(x - 2, x, domain='ZZ[a,b]'),
                (Poly(a - 1, x, domain='ZZ[a,b]'), Poly(b - 3, x, domain='ZZ[a,b]'))
            )
        ]],
        []
    )
    assert result == expected

    result = factor_system_poly([Poly(x ** 2 + 1, x, domain=QQ_I)])
    expected = ([[(Poly(x - I, x, domain='QQ_I'), ())], [(Poly(x + I, x, domain='QQ_I'), ())]], [])
    assert result == expected

    assert factor_system_poly([]) == ([[]], [])
    assert factor_system_poly([Poly(1, x), Poly(x, x)]) == ([], [])
    assert factor_system_poly([Poly(0, x), Poly(x, x)]) == ([[(Poly(x, x, domain='ZZ'), ())]], [])


def test_factor_system_cond():

    result = factor_system_cond([x ** 2 - 1, x ** 2 - 4])
    expected = (
        [
            [(x + 1, []), (x - 2, [])],
            [(x + 1, []), (x + 2, [])],
            [(x - 1, []), (x - 2, [])],
            [(x - 1, []), (x + 2, [])]
        ],
        []
    )
    assert result == expected

    assert factor_system_cond([1]) == ([], [])
    assert factor_system_cond([0]) == ([[]], [])
    assert factor_system_cond([1, x]) == ([], [])
    assert factor_system_cond([0, x]) == ([[(x, [])]], [])
    assert factor_system_cond([]) == ([[]], [])

    result = factor_system_cond([x ** 2 + y * x])
    expected = ([[(x, [])], [(x + y, [])]], [])
    assert result == expected

    result = factor_system_cond([(a - 1) * (x - 2), (b - 3) * (x - 2)], [x])
    expected = ([[(x - 2, [a - 1, b - 3])]], [])
    assert result == expected

    result = factor_system_cond([a * (x - 1), b], [x])
    expected = ([[(x - 1, [a])]], [b])
    assert result == expected

    result = factor_system_cond([a*x*(x-1), b*y, c], [x, y])
    expected = ([[(x, [a]), (y, [b])], [(x - 1, [a]), (y, [b])]], [c])
    assert result == expected

    result = factor_system_cond([x*(x-1), y], [x, y])
    expected = ([[(x, []), (y, [])], [(x - 1, []), (y, [])]], [])
    assert result == expected


def test_factor_system_bool():

    eqs = [a * (x - 1) * (y - 1), b * (x - 2) * (y - 1) * (y - 2)]
    expected = (Eq(y - 1, 0) |
                (Eq(a, 0) & Eq(b, 0)) |
                ((Eq(a, 0) | Eq(x - 1, 0)) & (Eq(b, 0) | Eq(x - 2, 0))) |
                ((Eq(a, 0) | Eq(x - 1, 0)) & (Eq(b, 0) | Eq(y - 2, 0))))
    assert factor_system_bool(eqs, [x, y]) == expected

    eqs = [x - 1]
    assert factor_system_bool(eqs, [x]) == Eq(x - 1, 0)

    eqs = [(x - 1) * (x - 2)]
    assert factor_system_bool(eqs, [x]) == Eq(x - 2, 0) | Eq(x - 1, 0)

    assert factor_system_bool([], [x]) == True
    assert factor_system_bool([0], [x]) == True
    assert factor_system_bool([1], [x]) == False
    assert factor_system_bool([a], [x]) == Eq(a, 0)

    eqs = [a * x * y, b * y * z]
    result = factor_system_bool(eqs, [x, y, z])
    expected = (
        Eq(y, 0)
        | (Eq(a, 0) & Eq(b, 0))
        | ((Eq(a, 0) | Eq(x, 0)) & (Eq(b, 0) | Eq(z, 0)))
    )
    assert result == expected

    eqs = [a * (x - 1), b]
    assert factor_system_bool(eqs, [x]) == Eq(b, 0) & (Eq(a, 0) | Eq(x - 1, 0))
