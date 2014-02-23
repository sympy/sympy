"""Tests for tools for solving inequalities and systems of inequalities. """

from sympy import (And, Eq, FiniteSet, Ge, Gt, im, Interval, Le, Lt, Ne, oo,
        Or, Q, re, S, sin, sqrt, Symbol, Union)
from sympy.assumptions import assuming
from sympy.abc import x, y
from sympy.solvers.inequalities import (reduce_inequalities,
                                        reduce_rational_inequalities,
                                        solve_univariate_inequality)
from sympy.utilities.pytest import raises

inf = oo.evalf()


def test_reduce_poly_inequalities_real_interval():
    with assuming(Q.real(x), Q.real(y)):
        assert reduce_rational_inequalities(
            [[Eq(x**2, 0)]], x, relational=False) == FiniteSet(0)
        assert reduce_rational_inequalities(
            [[Le(x**2, 0)]], x, relational=False) == FiniteSet(0)
        assert reduce_rational_inequalities(
            [[Lt(x**2, 0)]], x, relational=False) == S.EmptySet
        assert reduce_rational_inequalities(
            [[Ge(x**2, 0)]], x, relational=False) == Interval(-oo, oo)
        assert reduce_rational_inequalities(
            [[Gt(x**2, 0)]], x, relational=False) == FiniteSet(0).complement
        assert reduce_rational_inequalities(
            [[Ne(x**2, 0)]], x, relational=False) == FiniteSet(0).complement

        assert reduce_rational_inequalities(
            [[Eq(x**2, 1)]], x, relational=False) == FiniteSet(-1, 1)
        assert reduce_rational_inequalities(
            [[Le(x**2, 1)]], x, relational=False) == Interval(-1, 1)
        assert reduce_rational_inequalities(
            [[Lt(x**2, 1)]], x, relational=False) == Interval(-1, 1, True, True)
        assert reduce_rational_inequalities([[Ge(x**2, 1)]], x, relational=False) == Union(Interval(-oo, -1), Interval(1, oo))
        assert reduce_rational_inequalities(
            [[Gt(x**2, 1)]], x, relational=False) == Interval(-1, 1).complement
        assert reduce_rational_inequalities(
            [[Ne(x**2, 1)]], x, relational=False) == FiniteSet(-1, 1).complement
        assert reduce_rational_inequalities([[Eq(
            x**2, 1.0)]], x, relational=False) == FiniteSet(-1.0, 1.0).evalf()
        assert reduce_rational_inequalities(
            [[Le(x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0)
        assert reduce_rational_inequalities([[Lt(
            x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0, True, True)
        assert reduce_rational_inequalities([[Ge(x**2, 1.0)]], x, relational=False) == Union(Interval(-inf, -1.0), Interval(1.0, inf))
        assert reduce_rational_inequalities([[Gt(x**2, 1.0)]], x, relational=False) == Union(Interval(-inf, -1.0, right_open=True), Interval(1.0, inf, left_open=True))
        assert reduce_rational_inequalities([[Ne(
            x**2, 1.0)]], x, relational=False) == FiniteSet(-1.0, 1.0).complement

        s = sqrt(2)

        assert reduce_rational_inequalities([[Lt(
            x**2 - 1, 0), Gt(x**2 - 1, 0)]], x, relational=False) == S.EmptySet
        assert reduce_rational_inequalities([[Le(x**2 - 1, 0), Ge(
            x**2 - 1, 0)]], x, relational=False) == FiniteSet(-1, 1)
        assert reduce_rational_inequalities([[Le(x**2 - 2, 0), Ge(x**2 - 1, 0)]], x, relational=False) == Union(Interval(-s, -1, False, False), Interval(1, s, False, False))
        assert reduce_rational_inequalities([[Le(x**2 - 2, 0), Gt(x**2 - 1, 0)]], x, relational=False) == Union(Interval(-s, -1, False, True), Interval(1, s, True, False))
        assert reduce_rational_inequalities([[Lt(x**2 - 2, 0), Ge(x**2 - 1, 0)]], x, relational=False) == Union(Interval(-s, -1, True, False), Interval(1, s, False, True))
        assert reduce_rational_inequalities([[Lt(x**2 - 2, 0), Gt(x**2 - 1, 0)]], x, relational=False) == Union(Interval(-s, -1, True, True), Interval(1, s, True, True))
        assert reduce_rational_inequalities([[Lt(x**2 - 2, 0), Ne(x**2 - 1, 0)]], x, relational=False) == Union(Interval(-s, -1, True, True), Interval(-1, 1, True, True), Interval(1, s, True, True))



def test_reduce_poly_inequalities_real_relational():
    with assuming(Q.real(x), Q.real(y)):
        assert reduce_rational_inequalities(
            [[Eq(x**2, 0)]], x, relational=True) == Eq(x, 0)
        assert reduce_rational_inequalities(
            [[Le(x**2, 0)]], x, relational=True) == Eq(x, 0)
        assert reduce_rational_inequalities(
            [[Lt(x**2, 0)]], x, relational=True) is False
        assert reduce_rational_inequalities(
            [[Ge(x**2, 0)]], x, relational=True) is True
        assert reduce_rational_inequalities(
            [[Gt(x**2, 0)]], x, relational=True) == Or(Lt(x, 0), Gt(x, 0))
        assert reduce_rational_inequalities(
            [[Ne(x**2, 0)]], x, relational=True) == Or(Lt(x, 0), Gt(x, 0))

        assert reduce_rational_inequalities(
            [[Eq(x**2, 1)]], x, relational=True) == Or(Eq(x, -1), Eq(x, 1))
        assert reduce_rational_inequalities(
            [[Le(x**2, 1)]], x, relational=True) == And(Le(-1, x), Le(x, 1))
        assert reduce_rational_inequalities(
            [[Lt(x**2, 1)]], x, relational=True) == And(Lt(-1, x), Lt(x, 1))
        assert reduce_rational_inequalities(
            [[Ge(x**2, 1)]], x, relational=True) == Or(Le(x, -1), Ge(x, 1))
        assert reduce_rational_inequalities(
            [[Gt(x**2, 1)]], x, relational=True) == Or(Lt(x, -1), Gt(x, 1))
        assert reduce_rational_inequalities([[Ne(x**2, 1)]], x, relational=True) == Or(
            Lt(x, -1), And(Lt(-1, x), Lt(x, 1)), Gt(x, 1))

        assert reduce_rational_inequalities(
            [[Le(x**2, 1.0)]], x, relational=True) == And(Le(-1.0, x), Le(x, 1.0))
        assert reduce_rational_inequalities(
            [[Lt(x**2, 1.0)]], x, relational=True) == And(Lt(-1.0, x), Lt(x, 1.0))
        assert reduce_rational_inequalities(
            [[Ge(x**2, 1.0)]], x, relational=True) == Or(Le(x, -1.0), Ge(x, 1.0))
        assert reduce_rational_inequalities(
            [[Gt(x**2, 1.0)]], x, relational=True) == Or(Lt(x, -1.0), Gt(x, 1.0))
        assert reduce_rational_inequalities([[Ne(x**2, 1.0)]], x, relational=True) == \
            Or(Lt(x, -1.0), And(Lt(-1.0, x), Lt(x, 1.0)), Gt(x, 1.0))


def test_reduce_poly_inequalities_complex_relational():
    cond = Eq(im(x), 0)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=True) == And(Eq(re(x), 0), cond)
    assert reduce_rational_inequalities(
        [[Le(x**2, 0)]], x, relational=True) == And(Eq(re(x), 0), cond)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=True) == False
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=True) == cond
    assert reduce_rational_inequalities([[Gt(x**2, 0)]], x, relational=True) == \
        And(Or(Lt(re(x), 0), Gt(re(x), 0)), cond)
    assert reduce_rational_inequalities([[Ne(x**2, 0)]], x, relational=True) == \
        And(Or(Lt(re(x), 0), Gt(re(x), 0)), cond)

    assert reduce_rational_inequalities([[Eq(x**2, 1)]], x, relational=True) == \
        And(Or(Eq(re(x), -1), Eq(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Le(x**2, 1)]], x, relational=True) == \
        And(And(Le(-1, re(x)), Le(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Lt(x**2, 1)]], x, relational=True) == \
        And(And(Lt(-1, re(x)), Lt(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Ge(x**2, 1)]], x, relational=True) == \
        And(Or(Le(re(x), -1), Ge(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Gt(x**2, 1)]], x, relational=True) == \
        And(Or(Lt(re(x), -1), Gt(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Ne(x**2, 1)]], x, relational=True) == \
        And(Or(Lt(re(x), -1), And(Lt(-1, re(x)), Lt(re(x), 1)), Gt(re(x), 1)), cond)

    assert reduce_rational_inequalities([[Le(x**2, 1.0)]], x, relational=True) == \
        And(And(Le(-1.0, re(x)), Le(re(x), 1.0)), cond)
    assert reduce_rational_inequalities([[Lt(x**2, 1.0)]], x, relational=True) == \
        And(And(Lt(-1.0, re(x)), Lt(re(x), 1.0)), cond)
    assert reduce_rational_inequalities([[Ge(x**2, 1.0)]], x, relational=True) == \
        And(Or(Le(re(x), -1.0), Ge(re(x), 1.0)), cond)
    assert reduce_rational_inequalities([[Gt(x**2, 1.0)]], x, relational=True) == \
        And(Or(Lt(re(x), -1.0), Gt(re(x), 1.0)), cond)
    assert reduce_rational_inequalities([[Ne(x**2, 1.0)]], x, relational=True) == \
        And(Or(Lt(re(x), -1.0), And(Lt(-1.0, re(x)), Lt(re(x), 1.0)), Gt(re(x), 1.0)), cond)


def test_reduce_rational_inequalities_real_relational():
    def OpenInterval(a, b):
        return Interval(a, b, True, True)
    def LeftOpenInterval(a, b):
        return Interval(a, b, True, False)
    def RightOpenInterval(a, b):
        return Interval(a, b, False, True)

    with assuming(Q.real(x)):
        assert reduce_rational_inequalities([[(x**2 + 3*x + 2)/(x**2 - 16) >= 0]], x, relational=False) == \
            Union(OpenInterval(-oo, -4), Interval(-2, -1), OpenInterval(4, oo))

        assert reduce_rational_inequalities([[((-2*x - 10)*(3 - x))/((x**2 + 5)*(x - 2)**2) < 0]], x, relational=False) == \
            Union(OpenInterval(-5, 2), OpenInterval(2, 3))

        assert reduce_rational_inequalities([[(x + 1)/(x - 5) <= 0]], x, assume=Q.real(x), relational=False) == \
            RightOpenInterval(-1, 5)

        assert reduce_rational_inequalities([[(x**2 + 4*x + 3)/(x - 1) > 0]], x, assume=Q.real(x), relational=False) == \
            Union(OpenInterval(-3, -1), OpenInterval(1, oo))

        assert reduce_rational_inequalities([[(x**2 - 16)/(x - 1)**2 < 0]], x, assume=Q.real(x), relational=False) == \
            Union(OpenInterval(-4, 1), OpenInterval(1, 4))

        assert reduce_rational_inequalities([[(3*x + 1)/(x + 4) >= 1]], x, assume=Q.real(x), relational=False) == \
            Union(OpenInterval(-oo, -4), RightOpenInterval(S(3)/2, oo))

        assert reduce_rational_inequalities([[(x - 8)/x <= 3 - x]], x, assume=Q.real(x), relational=False) == \
            Union(LeftOpenInterval(-oo, -2), LeftOpenInterval(0, 4))


def test_reduce_abs_inequalities():
    real = Q.real(x)

    assert reduce_inequalities(
        abs(x - 5) < 3, assume=real) == And(Lt(2, x), Lt(x, 8))
    assert reduce_inequalities(
        abs(2*x + 3) >= 8, assume=real) == Or(Le(x, -S(11)/2), Ge(x, S(5)/2))
    assert reduce_inequalities(abs(x - 4) + abs(
        3*x - 5) < 7, assume=real) == And(Lt(S(1)/2, x), Lt(x, 4))
    assert reduce_inequalities(abs(x - 4) + abs(3*abs(x) - 5) < 7, assume=real) == Or(And(S(-2) < x, x < -1), And(S(1)/2 < x, x < 4))

    raises(NotImplementedError, lambda: reduce_inequalities(abs(x - 5) < 3))


def test_reduce_inequalities_boolean():
    assert reduce_inequalities(
        [Eq(x**2, 0), True]) == And(Eq(re(x), 0), Eq(im(x), 0))
    assert reduce_inequalities([Eq(x**2, 0), False]) is False


def test_reduce_inequalities_assume():
    assert reduce_inequalities(
        [Le(x**2, 1), Q.real(x)]) == And(Le(-1, x), Le(x, 1))
    assert reduce_inequalities(
        [Le(x**2, 1)], Q.real(x)) == And(Le(-1, x), Le(x, 1))


def test_reduce_inequalities_multivariate():
    assert reduce_inequalities([Ge(x**2, 1), Ge(y**2, 1)]) == \
        And(And(Or(Le(re(x), -1), Ge(re(x), 1)), Eq(im(x), 0)),
            And(Or(Le(re(y), -1), Ge(re(y), 1)), Eq(im(y), 0)))


def test_reduce_inequalities_errors():
    raises(NotImplementedError, lambda: reduce_inequalities(Ge(sin(x) + x, 1)))
    raises(NotImplementedError, lambda: reduce_inequalities(Ge(x**2*y + y, 1)))
    raises(NotImplementedError, lambda: reduce_inequalities(Ge(sqrt(2)*x, 1)))


def test_hacky_inequalities():
    assert reduce_inequalities(x + y < 1, symbols=[x]) == (x < 1 - y)
    assert reduce_inequalities(x + y >= 1, symbols=[x]) == (x >= 1 - y)


def test_issue_3244():
    eq = -3*x**2/2 - 45*x/4 + S(33)/2 > 0
    assert reduce_inequalities(eq, Q.real(x)) == \
        And(x < -S(15)/4 + sqrt(401)/4, -sqrt(401)/4 - S(15)/4 < x)


def test_solve_univariate_inequality():
    x = Symbol('x', real=True)
    isolve = solve_univariate_inequality
    assert isolve(x**2 >= 4, x, relational=False) == Union(Interval(-oo, -2), Interval(2, oo))
    assert isolve(x**2 >= 4, x) == Or(x <= -2, x >= 2)
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x, relational=False) == \
        Union(Interval(1, 2), Interval(3, oo))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x) == \
        Or(And(S.One <= x, x <= 2), x >= 3)
    # github's issue #2785:
    assert isolve(x**3 - 2*x - 1 > 0, x, relational=False) == \
        Union(Interval(-1, -sqrt(5)/2 + S(1)/2, True, True),
              Interval(S(1)/2 + sqrt(5)/2, oo, True, True))
    # github's issue #2794:
    assert isolve(x**3 - x**2 + x - 1 > 0, x, relational=False) == \
        Interval(1, oo, True)
