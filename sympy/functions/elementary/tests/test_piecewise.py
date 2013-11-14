from sympy import (
    adjoint, And, Basic, conjugate, diff, expand, Eq, Function, I, im,
    Integral, integrate, Interval, lambdify, log, Max, Min, oo, Or, pi,
    Piecewise, piecewise_fold, Rational, re, solve, symbols, transpose,
    cos, exp, Abs, Not
)
from sympy.utilities.pytest import XFAIL, raises

x, y = symbols('x y')
z = symbols('z', nonzero=True)


def test_piecewise():

    # Test canonization
    assert Piecewise((x, x < 1), (0, True)) == Piecewise((x, x < 1), (0, True))
    assert Piecewise((x, x < 1), (0, True), (1, True)) == \
        Piecewise((x, x < 1), (0, True))
    assert Piecewise((x, x < 1), (0, False), (-1, 1 > 2)) == \
        Piecewise((x, x < 1))
    assert Piecewise((x, x < 1), (0, x < 1), (0, True)) == \
        Piecewise((x, x < 1), (0, True))
    assert Piecewise((x, x < 1), (0, x < 2), (0, True)) == \
        Piecewise((x, x < 1), (0, True))
    assert Piecewise((x, x < 1), (x, x < 2), (0, True)) == \
        Piecewise((x, Or(x < 1, x < 2)), (0, True))
    assert Piecewise((x, x < 1), (x, x < 2), (x, True)) == x
    assert Piecewise((x, True)) == x
    raises(TypeError, lambda: Piecewise(x))
    raises(TypeError, lambda: Piecewise((x, x**2)))

    # Test subs
    p = Piecewise((-1, x < -1), (x**2, x < 0), (log(x), x >= 0))
    p_x2 = Piecewise((-1, x**2 < -1), (x**4, x**2 < 0), (log(x**2), x**2 >= 0))
    assert p.subs(x, x**2) == p_x2
    assert p.subs(x, -5) == -1
    assert p.subs(x, -1) == 1
    assert p.subs(x, 1) == log(1)

    # More subs tests
    p2 = Piecewise((1, x < pi), (-1, x < 2*pi), (0, x > 2*pi))
    p3 = Piecewise((1, Eq(x, 0)), (1/x, True))
    p4 = Piecewise((1, Eq(x, 0)), (2, 1/x>2))
    assert p2.subs(x, 2) == 1
    assert p2.subs(x, 4) == -1
    assert p2.subs(x, 10) == 0
    assert p3.subs(x, 0.0) == 1
    assert p4.subs(x, 0.0) == 1


    f, g, h = symbols('f,g,h', cls=Function)
    pf = Piecewise((f(x), x < -1), (f(x) + h(x) + 2, x <= 1))
    pg = Piecewise((g(x), x < -1), (g(x) + h(x) + 2, x <= 1))
    assert pg.subs(g, f) == pf

    assert Piecewise((1, Eq(x, 0)), (0, True)).subs(x, 0) == 1
    assert Piecewise((1, Eq(x, 0)), (0, True)).subs(x, 1) == 0
    assert Piecewise((1, Eq(x, y)), (0, True)).subs(x, y) == 1
    assert Piecewise((1, Eq(x, z)), (0, True)).subs(x, z) == 1
    assert Piecewise((1, Eq(exp(x), cos(z))), (0, True)).subs(x, z) == \
        Piecewise((1, Eq(exp(z), cos(z))), (0, True))
    assert Piecewise((1, Eq(x, y*(y + 1))), (0, True)).subs(x, y**2 + y) == 1

    p5 = Piecewise( (0, Eq(cos(x) + y, 0)), (1, True))
    assert p5.subs(y, 0) == Piecewise( (0, Eq(cos(x), 0)), (1, True))

    # Test evalf
    assert p.evalf() == p
    assert p.evalf(subs={x: -2}) == -1
    assert p.evalf(subs={x: -1}) == 1
    assert p.evalf(subs={x: 1}) == log(1)

    # Test doit
    f_int = Piecewise((Integral(x, (x, 0, 1)), x < 1))
    assert f_int.doit() == Piecewise( (1.0/2.0, x < 1) )

    # Test differentiation
    f = x
    fp = x*p
    dp = Piecewise((0, x < -1), (2*x, x < 0), (1/x, x >= 0))
    fp_dx = x*dp + p
    assert diff(p, x) == dp
    assert diff(f*p, x) == fp_dx

    # Test simple arithmetic
    assert x*p == fp
    assert x*p + p == p + x*p
    assert p + f == f + p
    assert p + dp == dp + p
    assert p - dp == -(dp - p)

    # Test power
    dp2 = Piecewise((0, x < -1), (4*x**2, x < 0), (1/x**2, x >= 0))
    assert dp**2 == dp2

    # Test _eval_interval
    f1 = x*y + 2
    f2 = x*y**2 + 3
    peval = Piecewise( (f1, x < 0), (f2, x > 0))
    peval_interval = f1.subs(
        x, 0) - f1.subs(x, -1) + f2.subs(x, 1) - f2.subs(x, 0)
    assert peval._eval_interval(x, 0, 0) == 0
    assert peval._eval_interval(x, -1, 1) == peval_interval
    peval2 = Piecewise((f1, x < 0), (f2, True))
    assert peval2._eval_interval(x, 0, 0) == 0
    assert peval2._eval_interval(x, 1, -1) == -peval_interval
    assert peval2._eval_interval(x, -1, -2) == f1.subs(x, -2) - f1.subs(x, -1)
    assert peval2._eval_interval(x, -1, 1) == peval_interval
    assert peval2._eval_interval(x, None, 0) == peval2.subs(x, 0)
    assert peval2._eval_interval(x, -1, None) == -peval2.subs(x, -1)

    # Test integration
    p_int = Piecewise((-x, x < -1), (x**3/3.0, x < 0), (-x + x*log(x), x >= 0))
    assert integrate(p, x) == p_int
    p = Piecewise((x, x < 1), (x**2, -1 <= x), (x, 3 < x))
    assert integrate(p, (x, -2, 2)) == 5.0/6.0
    assert integrate(p, (x, 2, -2)) == -5.0/6.0
    p = Piecewise((0, x < 0), (1, x < 1), (0, x < 2), (1, x < 3), (0, True))
    assert integrate(p, (x, -oo, oo)) == 2
    p = Piecewise((x, x < -10), (x**2, x <= -1), (x, 1 < x))
    raises(ValueError, lambda: integrate(p, (x, -2, 2)))

    # Test commutativity
    assert p.is_commutative is True


def test_piecewise_free_symbols():
    a = symbols('a')
    f = Piecewise((x, a < 0), (y, True))
    assert f.free_symbols == set([x, y, a])


def test_piecewise_integrate():
    # XXX Use '<=' here! '>=' is not yet implemented ..
    f = Piecewise(((x - 2)**2, 0 <= x), (1, True))
    assert integrate(f, (x, -2, 2)) == Rational(14, 3)

    g = Piecewise(((x - 5)**5, 4 <= x), (f, True))
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == Rational(43, 6)

    g = Piecewise(((x - 5)**5, 4 <= x), (f, x < 4))
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == Rational(43, 6)

    g = Piecewise(((x - 5)**5, 2 <= x), (f, x < 2))
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == -Rational(701, 6)

    g = Piecewise(((x - 5)**5, 2 <= x), (f, True))
    assert integrate(g, (x, -2, 2)) == Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == -Rational(701, 6)

    g = Piecewise(((x - 5)**5, 2 <= x), (2 * f, True))
    assert integrate(g, (x, -2, 2)) == 2 * Rational(14, 3)
    assert integrate(g, (x, -2, 5)) == -Rational(673, 6)

    g = Piecewise((1, x > 0), (0, Eq(x, 0)), (-1, x < 0))
    assert integrate(g, (x, -1, 1)) == 0

    g = Piecewise((1, x - y < 0), (0, True))
    assert integrate(g, (y, -oo, 0)) == -Min(0, x)
    assert integrate(g, (y, 0, oo)) == oo - Max(0, x)
    assert integrate(g, (y, -oo, oo)) == oo - x

    g = Piecewise((0, x < 0), (x, x <= 1), (1, True))
    assert integrate(g, (x, -5, 1)) == Rational(1, 2)
    assert integrate(g, (x, -5, y)).subs(y, 1) == Rational(1, 2)
    assert integrate(g, (x, y, 1)).subs(y, -5) == Rational(1, 2)
    assert integrate(g, (x, 1, -5)) == -Rational(1, 2)
    assert integrate(g, (x, 1, y)).subs(y, -5) == -Rational(1, 2)
    assert integrate(g, (x, y, -5)).subs(y, 1) == -Rational(1, 2)
    assert integrate(g, (x, -5, y)) == Piecewise((0, y < 0),
        (y**2/2, y <= 1), (y - 0.5, True))
    assert integrate(g, (x, y, 1)) == Piecewise((0.5, y < 0),
        (0.5 - y**2/2, y <= 1), (1 - y, True))

    g = Piecewise((1 - x, Interval(0, 1).contains(x)),
        (1 + x, Interval(-1, 0).contains(x)), (0, True))
    assert integrate(g, (x, -5, 1)) == 1
    assert integrate(g, (x, -5, y)).subs(y, 1) == 1
    assert integrate(g, (x, y, 1)).subs(y, -5) == 1
    assert integrate(g, (x, 1, -5)) == -1
    assert integrate(g, (x, 1, y)).subs(y, -5) == -1
    assert integrate(g, (x, y, -5)).subs(y, 1) == -1
    assert integrate(g, (x, -5, y)) == Piecewise(
        (-y**2/2 + y + 0.5, Interval(0, 1).contains(y)),
        (y**2/2 + y + 0.5, Interval(-1, 0).contains(y)),
        (0, y <= -1), (1, True))
    assert integrate(g, (x, y, 1)) == Piecewise(
        (y**2/2 - y + 0.5, Interval(0, 1).contains(y)),
        (-y**2/2 - y + 0.5, Interval(-1, 0).contains(y)),
        (1, y <= -1), (0, True))

    g = Piecewise((0, Or(x <= -1, x >= 1)), (1 - x, x > 0), (1 + x, True))
    assert integrate(g, (x, -5, 1)) == 1
    assert integrate(g, (x, -5, y)).subs(y, 1) == 1
    assert integrate(g, (x, y, 1)).subs(y, -5) == 1
    assert integrate(g, (x, 1, -5)) == -1
    assert integrate(g, (x, 1, y)).subs(y, -5) == -1
    assert integrate(g, (x, y, -5)).subs(y, 1) == -1
    assert integrate(g, (x, -5, y)) == Piecewise((0, y <= -1), (1, y >= 1),
        (-y**2/2 + y + 0.5, y > 0), (y**2/2 + y + 0.5, True))
    assert integrate(g, (x, y, 1)) == Piecewise((1, y <= -1), (0, y >= 1),
        (y**2/2 - y + 0.5, y > 0), (-y**2/2 - y + 0.5, True))


def test_piecewise_integrate_inequality_conditions():
    c1, c2 = symbols("c1 c2", positive=True)
    g = Piecewise((0, c1*x > 1), (1, c1*x > 0), (0, True))
    assert integrate(g, (x, -oo, 0)) == 0
    assert integrate(g, (x, -5, 0)) == 0
    assert integrate(g, (x, 0, 5)) == Min(5, 1/c1)
    assert integrate(g, (x, 0, oo)) == 1/c1

    g = Piecewise((0, c1*x + c2*y > 1), (1, c1*x + c2*y > 0), (0, True))
    assert integrate(g, (x, -oo, 0)).subs(y, 0) == 0
    assert integrate(g, (x, -5, 0)).subs(y, 0) == 0
    assert integrate(g, (x, 0, 5)).subs(y, 0) == Min(5, 1/c1)
    assert integrate(g, (x, 0, oo)).subs(y, 0) == 1/c1


def test_piecewise_integrate_symbolic_conditions():
    from sympy.abc import a, b, x, y
    p0 = Piecewise((0, Or(x < a, x > b)), (1, True))
    p1 = Piecewise((0, x < a), (0, x > b), (1, True))
    p2 = Piecewise((0, x > b), (0, x < a), (1, True))
    p3 = Piecewise((0, x < a), (1, x < b), (0, True))
    p4 = Piecewise((0, x > b), (1, x > a), (0, True))
    p5 = Piecewise((1, And(a < x, x < b)), (0, True))
    assert integrate(p0, (x, -oo, y)) == Min(b, y) - Min(a, b, y)
    assert integrate(p1, (x, -oo, y)) == Min(b, y) - Min(a, b, y)
    assert integrate(p2, (x, -oo, y)) == Min(b, y) - Min(a, b, y)
    assert integrate(p3, (x, -oo, y)) == Min(b, y) - Min(a, b, y)
    assert integrate(p4, (x, -oo, y)) == Min(b, y) - Min(a, b, y)
    assert integrate(p5, (x, -oo, y)) == Min(b, y) - Min(a, b, y)
    assert integrate(p0, (x, y, oo)) == Max(a, b, y) - Max(a, y)
    assert integrate(p1, (x, y, oo)) == Max(a, b, y) - Max(a, y)
    assert integrate(p2, (x, y, oo)) == Max(a, b, y) - Max(a, y)
    assert integrate(p3, (x, y, oo)) == Max(a, b, y) - Max(a, y)
    assert integrate(p4, (x, y, oo)) == Max(a, b, y) - Max(a, y)
    assert integrate(p5, (x, y, oo)) == Max(a, b, y) - Max(a, y)

    assert integrate(p0, x) == Piecewise((0, Or(x < a, x > b)), (x, True))
    assert integrate(p1, x) == Piecewise((0, Or(x < a, x > b)), (x, True))
    assert integrate(p2, x) == Piecewise((0, Or(x < a, x > b)), (x, True))

    p1 = Piecewise((0, x < a), (0.5, x > b), (1, True))
    p2 = Piecewise((0.5, x > b), (0, x < a), (1, True))
    p3 = Piecewise((0, x < a), (1, x < b), (0.5, True))
    p4 = Piecewise((0.5, x > b), (1, x > a), (0, True))
    p5 = Piecewise((1, And(a < x, x < b)), (0.5, x > b), (0, True))
    assert integrate(p1, (x, -oo, y)) == 0.5*y + 0.5*Min(b, y) - Min(a, b, y)
    assert integrate(p2, (x, -oo, y)) == 0.5*y + 0.5*Min(b, y) - Min(a, b, y)
    assert integrate(p3, (x, -oo, y)) == 0.5*y + 0.5*Min(b, y) - Min(a, b, y)
    assert integrate(p4, (x, -oo, y)) == 0.5*y + 0.5*Min(b, y) - Min(a, b, y)
    assert integrate(p5, (x, -oo, y)) == 0.5*y + 0.5*Min(b, y) - Min(a, b, y)


def test_piecewise_integrate_independent_conditions():
    p = Piecewise((0, Eq(y, 0)), (x*y, True))
    assert integrate(p, (x, 1, 3)) == \
        Piecewise((0, Eq(y, 0)), (4*y, True))


def test_piecewise_simplify():
    p = Piecewise(((x**2 + 1)/x**2, Eq(x*(1 + x) - x**2, 0)),
                  ((-1)**x*(-1), True))
    assert p.simplify() == \
        Piecewise((1 + 1/x**2, Eq(x, 0)), ((-1)**(x + 1), True))


def test_piecewise_solve():
    abs2 = Piecewise((-x, x <= 0), (x, x > 0))
    f = abs2.subs(x, x - 2)
    assert solve(f, x) == [2]
    assert solve(f - 1, x) == [1, 3]

    f = Piecewise(((x - 2)**2, x >= 0), (1, True))
    assert solve(f, x) == [2]

    g = Piecewise(((x - 5)**5, x >= 4), (f, True))
    assert solve(g, x) == [2, 5]

    g = Piecewise(((x - 5)**5, x >= 4), (f, x < 4))
    assert solve(g, x) == [2, 5]

    g = Piecewise(((x - 5)**5, x >= 2), (f, x < 2))
    assert solve(g, x) == [5]

    g = Piecewise(((x - 5)**5, x >= 2), (f, True))
    assert solve(g, x) == [5]

    g = Piecewise(((x - 5)**5, x >= 2), (f, True), (10, False))
    assert solve(g, x) == [5]

    g = Piecewise(((x - 5)**5, x >= 2),
                  (-x + 2, x - 2 <= 0), (x - 2, x - 2 > 0))
    assert solve(g, x) == [5]

# See issue 1253 (enhance the solver to handle inequalities).


@XFAIL
def test_piecewise_solve2():
    f = Piecewise(((x - 2)**2, x >= 0), (0, True))
    assert solve(f, x) == [2, Interval(0, oo, True, True)]


def test_piecewise_fold():

    p = Piecewise((x, x < 1), (1, 1 <= x))

    assert piecewise_fold(x*p) == Piecewise((x**2, x < 1), (x, 1 <= x))
    assert piecewise_fold(p + p) == Piecewise((2*x, x < 1), (2, 1 <= x))
    assert piecewise_fold(Piecewise((1, x < 0), (2, True))
                          + Piecewise((10, x < 0), (-10, True))) == \
        Piecewise((11, x < 0), (-8, True))

    p1 = Piecewise((0, x < 0), (x, x <= 1), (0, True))
    p2 = Piecewise((0, x < 0), (1 - x, x <= 1), (0, True))

    p = 4*p1 + 2*p2
    assert integrate(
        piecewise_fold(p), (x, -oo, oo)) == integrate(2*x + 2, (x, 0, 1))


def test_piecewise_fold_piecewise_in_cond():
    p1 = Piecewise((cos(x), x < 0), (0, True))
    p2 = Piecewise((0, Eq(p1, 0)), (p1 / Abs(p1), True))
    p3 = piecewise_fold(p2)
    assert(p2.subs(x, -pi/2) == 0.0)
    assert(p2.subs(x, 1) == 0.0)
    assert(p2.subs(x, -pi/4) == 1.0)
    p4 = Piecewise((0, Eq(p1, 0)), (1,True))
    assert(piecewise_fold(p4) == Piecewise(
        (0, Or(And(Eq(cos(x), 0), x < 0), Not(x < 0))), (1, True)))

    r1 = 1 < Piecewise((1, x < 1), (3, True))
    assert(piecewise_fold(r1) == Not(x < 1))

    p5 = Piecewise((1, x < 0), (3, True))
    p6 = Piecewise((1, x < 1), (3, True))
    p7 = piecewise_fold(Piecewise((1, p5 < p6), (0, True)))
    assert(Piecewise((1, And(Not(x < 1), x < 0)), (0, True)))


@XFAIL
def test_piecewise_fold_piecewise_in_cond_2():
    p1 = Piecewise((cos(x), x < 0), (0, True))
    p2 = Piecewise((0, Eq(p1, 0)), (1 / p1, True))
    p3 = Piecewise((0, Or(And(Eq(cos(x), 0), x < 0), Not(x < 0))),
        (1 / cos(x), True))
    assert(piecewise_fold(p2) == p3)


def test_piecewise_fold_expand():
    p1 = Piecewise((1, Interval(0, 1, False, True).contains(x)), (0, True))

    p2 = piecewise_fold(expand((1 - x)*p1))
    assert p2 == Piecewise((1 - x, Interval(0, 1, False, True).contains(x)),
        (Piecewise((-x, Interval(0, 1, False, True).contains(x)), (0, True)), True))

    p2 = expand(piecewise_fold((1 - x)*p1))
    assert p2 == Piecewise(
        (1 - x, Interval(0, 1, False, True).contains(x)), (0, True))


def test_piecewise_duplicate():
    p = Piecewise((x, x < -10), (x**2, x <= -1), (x, 1 < x))
    assert p == Piecewise(*p.args)


def test_doit():
    p1 = Piecewise((x, x < 1), (x**2, -1 <= x), (x, 3 < x))
    p2 = Piecewise((x, x < 1), (Integral(2 * x), -1 <= x), (x, 3 < x))
    assert p2.doit() == p1
    assert p2.doit(deep=False) == p2


def test_piecewise_interval():
    p1 = Piecewise((x, Interval(0, 1).contains(x)), (0, True))
    assert p1.subs(x, -0.5) == 0
    assert p1.subs(x, 0.5) == 0.5
    assert p1.diff(x) == Piecewise((1, Interval(0, 1).contains(x)), (0, True))
    assert integrate(
        p1, x) == Piecewise((x**2/2, Interval(0, 1).contains(x)), (0, True))


def test_piecewise_collapse():
    p1 = Piecewise((x, x < 0), (x**2, x > 1))
    p2 = Piecewise((p1, x < 0), (p1, x > 1))
    assert p2 == Piecewise((x, x < 0), (x**2, 1 < x))

    p1 = Piecewise((Piecewise((x, x < 0), (1, True)), True))
    assert p1 == Piecewise((Piecewise((x, x < 0), (1, True)), True))


def test_piecewise_lambdify():
    p = Piecewise(
        (x**2, x < 0),
        (x, Interval(0, 1, False, True).contains(x)),
        (2 - x, x >= 1),
        (0, True)
    )

    f = lambdify(x, p)
    assert f(-2.0) == 4.0
    assert f(0.0) == 0.0
    assert f(0.5) == 0.5
    assert f(2.0) == 0.0


def test_piecewise_series():
    from sympy import sin, cos, O
    p1 = Piecewise((sin(x), x < 0), (cos(x), x > 0))
    p2 = Piecewise((x + O(x**2), x < 0), (1 + O(x**2), x > 0))
    assert p1.nseries(x, n=2) == p2


def test_piecewise_as_leading_term():
    p1 = Piecewise((1/x, x > 1), (0, True))
    p2 = Piecewise((x, x > 1), (0, True))
    p3 = Piecewise((1/x, x > 1), (x, True))
    p4 = Piecewise((x, x > 1), (1/x, True))
    p5 = Piecewise((1/x, x > 1), (x, True))
    p6 = Piecewise((1/x, x < 1), (x, True))
    p7 = Piecewise((x, x < 1), (1/x, True))
    p8 = Piecewise((x, x > 1), (1/x, True))
    assert p1.as_leading_term(x) == 0
    assert p2.as_leading_term(x) == 0
    assert p3.as_leading_term(x) == x
    assert p4.as_leading_term(x) == 1/x
    assert p5.as_leading_term(x) == x
    assert p6.as_leading_term(x) == 1/x
    assert p7.as_leading_term(x) == x
    assert p8.as_leading_term(x) == 1/x


def test_piecewise_complex():
    p1 = Piecewise((2, x < 0), (1, 0 <= x))
    p2 = Piecewise((2*I, x < 0), (I, 0 <= x))
    p3 = Piecewise((I*x, x > 1), (1 + I, True))
    p4 = Piecewise((-I*conjugate(x), x > 1), (1 - I, True))

    assert conjugate(p1) == p1
    assert conjugate(p2) == piecewise_fold(-p2)
    assert conjugate(p3) == p4

    assert p1.is_imaginary is False
    assert p1.is_real is True
    assert p2.is_imaginary is True
    assert p2.is_real is False
    assert p3.is_imaginary is None
    assert p3.is_real is None

    assert p1.as_real_imag() == (p1, 0)
    assert p2.as_real_imag() == (0, -I*p2)


def test_conjugate_transpose():
    A, B = symbols("A B", commutative=False)
    p = Piecewise((A*B**2, x > 0), (A**2*B, True))
    assert p.adjoint() == \
        Piecewise((adjoint(A*B**2), x > 0), (adjoint(A**2*B), True))
    assert p.conjugate() == \
        Piecewise((conjugate(A*B**2), x > 0), (conjugate(A**2*B), True))
    assert p.transpose() == \
        Piecewise((transpose(A*B**2), x > 0), (transpose(A**2*B), True))


def test_piecewise_evaluate():
    assert Piecewise((x, True)) == x
    assert Piecewise((x, True), evaluate=True) == x
    p = Piecewise((x, True), evaluate=False)
    assert p != x
    assert p.is_Piecewise
    assert all(isinstance(i, Basic) for i in p.args)
