"""Tests for tools for solving inequalities and systems of inequalities. """

from sympy import (And, Eq, FiniteSet, Ge, Gt, im, Interval, Le, Lt, Ne, oo,
                   Or, Q, re, S, sin, sqrt, Symbol, Union, Integral, Sum,
                   Function, Float)
from sympy.solvers.inequalities import (reduce_inequalities,
                                        reduce_rational_inequalities,
                                        solve_univariate_inequality as isolve)
from sympy.utilities.pytest import raises
from sympy.polys.rootoftools import RootOf
from sympy.solvers.solvers import solve

inf = oo.evalf()


def test_reduce_poly_inequalities_real_interval():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[Le(x**2, 0)]], x, relational=False) == FiniteSet(0)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=False) == S.EmptySet
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=False) == Interval(-oo, oo)
    assert reduce_rational_inequalities(
        [[Gt(x**2, 0)]], x, relational=False) == FiniteSet(0).complement(S.Reals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=False) == FiniteSet(0).complement(S.Reals)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 1)]], x, relational=False) == FiniteSet(-1, 1)
    assert reduce_rational_inequalities(
        [[Le(x**2, 1)]], x, relational=False) == Interval(-1, 1)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 1)]], x, relational=False) == Interval(-1, 1, True, True)
    assert reduce_rational_inequalities([[Ge(x**2, 1)]], x, relational=False) == Union(Interval(-oo, -1), Interval(1, oo))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 1)]], x, relational=False) == Interval(-1, 1).complement(S.Reals)
    assert reduce_rational_inequalities(
        [[Ne(x**2, 1)]], x, relational=False) == FiniteSet(-1, 1).complement(S.Reals)
    assert reduce_rational_inequalities([[Eq(
        x**2, 1.0)]], x, relational=False) == FiniteSet(-1.0, 1.0).evalf()
    assert reduce_rational_inequalities(
        [[Le(x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0)
    assert reduce_rational_inequalities([[Lt(
        x**2, 1.0)]], x, relational=False) == Interval(-1.0, 1.0, True, True)
    assert reduce_rational_inequalities([[Ge(x**2, 1.0)]], x, relational=False) == Union(Interval(-inf, -1.0), Interval(1.0, inf))
    assert reduce_rational_inequalities([[Gt(x**2, 1.0)]], x, relational=False) == Union(Interval(-inf, -1.0, right_open=True), Interval(1.0, inf, left_open=True))
    assert reduce_rational_inequalities([[Ne(
        x**2, 1.0)]], x, relational=False) == FiniteSet(-1.0, 1.0).complement(S.Reals)

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
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=True) == Eq(x, 0)
    assert reduce_rational_inequalities(
        [[Le(x**2, 0)]], x, relational=True) == Eq(x, 0)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=True) == False
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=True) == And(Lt(-oo, x), Lt(x, oo))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 0)]], x, relational=True) == Or(And(Lt(-oo, x), Lt(x, 0)), And(Lt(0, x), Lt(x, oo)))
    assert reduce_rational_inequalities(
        [[Ne(x**2, 0)]], x, relational=True) == Or(And(Lt(-oo, x), Lt(x, 0)), And(Lt(0, x), Lt(x, oo)))

    assert reduce_rational_inequalities(
        [[Eq(x**2, 1)]], x, relational=True) == Or(Eq(x, -1), Eq(x, 1))
    assert reduce_rational_inequalities(
        [[Le(x**2, 1)]], x, relational=True) == And(Le(-1, x), Le(x, 1))
    assert reduce_rational_inequalities(
        [[Lt(x**2, 1)]], x, relational=True) == And(Lt(-1, x), Lt(x, 1))
    assert reduce_rational_inequalities(
        [[Ge(x**2, 1)]], x, relational=True) == Or(And(Le(1, x), Lt(x, oo)), And(Le(x, -1), Lt(-oo, x)))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 1)]], x, relational=True) == Or(And(Lt(1, x), Lt(x, oo)), And(Lt(x, -1), Lt(-oo, x)))
    assert reduce_rational_inequalities([[Ne(x**2, 1)]], x, relational=True) == Or(
            And(Lt(-oo, x), Lt(x, -1)), And(Lt(-1, x), Lt(x, 1)), And(Lt(1, x), Lt(x, oo)))

    assert reduce_rational_inequalities(
        [[Le(x**2, 1.0)]], x, relational=True) == And(Le(-1.0, x), Le(x, 1.0))
    assert reduce_rational_inequalities(
        [[Lt(x**2, 1.0)]], x, relational=True) == And(Lt(-1.0, x), Lt(x, 1.0))
    assert reduce_rational_inequalities(
        [[Ge(x**2, 1.0)]], x, relational=True) == Or(And(Lt(Float('-inf'), x), Le(x, -1.0)),
                                                         And(Le(1.0, x), Lt(x, Float('+inf'))))
    assert reduce_rational_inequalities(
        [[Gt(x**2, 1.0)]], x, relational=True) == Or(And(Lt(Float('-inf'), x), Lt(x, -1.0)),
                                                         And(Lt(1.0, x), Lt(x, Float('+inf'))))
    assert reduce_rational_inequalities([[Ne(x**2, 1.0)]], x, relational=True) == \
            Or(And(Lt(-1.0, x), Lt(x, 1.0)), And(Lt(Float('-inf'), x), Lt(x, -1.0)),
               And(Lt(1.0, x), Lt(x, Float('+inf'))))


def test_reduce_poly_inequalities_complex_relational():
    x = Symbol('x')
    cond = Eq(im(x), 0)

    assert reduce_rational_inequalities(
        [[Eq(x**2, 0)]], x, relational=True) == And(Eq(re(x), 0), cond)
    assert reduce_rational_inequalities(
        [[Le(x**2, 0)]], x, relational=True) == And(Eq(re(x), 0), cond)
    assert reduce_rational_inequalities(
        [[Lt(x**2, 0)]], x, relational=True) == False
    assert reduce_rational_inequalities(
        [[Ge(x**2, 0)]], x, relational=True) == And(cond, Lt(-oo, re(x)), Lt(re(x), oo))
    assert reduce_rational_inequalities([[Gt(x**2, 0)]], x, relational=True) == \
        And(cond, Or(And(Lt(-oo, re(x)), Lt(re(x), 0)), And(Lt(0, re(x)), Lt(re(x), oo))))
    assert reduce_rational_inequalities([[Ne(x**2, 0)]], x, relational=True) == \
        And(cond, Or(And(Lt(-oo, re(x)), Lt(re(x), 0)), And(Lt(0, re(x)), Lt(re(x), oo))))

    assert reduce_rational_inequalities([[Eq(x**2, 1)]], x, relational=True) == \
        And(Or(Eq(re(x), -1), Eq(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Le(x**2, 1)]], x, relational=True) == \
        And(And(Le(-1, re(x)), Le(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Lt(x**2, 1)]], x, relational=True) == \
        And(And(Lt(-1, re(x)), Lt(re(x), 1)), cond)
    assert reduce_rational_inequalities([[Ge(x**2, 1)]], x, relational=True) == \
        And(cond, Or(And(Le(1, re(x)), Lt(re(x), oo)), And(Le(re(x), -1), Lt(-oo, re(x)))))
    assert reduce_rational_inequalities([[Gt(x**2, 1)]], x, relational=True) == \
        And(cond, Or(And(Lt(-oo, re(x)), Lt(re(x), -1)), And(Lt(1, re(x)), Lt(re(x), oo))))
    assert reduce_rational_inequalities([[Ne(x**2, 1)]], x, relational=True) == \
        And(cond, Or(And(Lt(-oo, re(x)), Lt(re(x), -1)),
                     And(Lt(-1, re(x)), Lt(re(x), 1)), And(Lt(1, re(x)), Lt(re(x), oo))))

    assert reduce_rational_inequalities([[Le(x**2, 1.0)]], x, relational=True) == \
        And(And(Le(-1.0, re(x)), Le(re(x), 1.0)), cond)
    assert reduce_rational_inequalities([[Lt(x**2, 1.0)]], x, relational=True) == \
        And(And(Lt(-1.0, re(x)), Lt(re(x), 1.0)), cond)
    assert reduce_rational_inequalities([[Ge(x**2, 1.0)]], x, relational=True) == \
        And(cond, Or(And(Le(1.0, re(x)), re(x) < Float('+inf')),
                     And(Le(re(x), -1.0), Float('-inf') < re(x))))
    assert reduce_rational_inequalities([[Gt(x**2, 1.0)]], x, relational=True) == \
        And(cond, Or(And(Float('-inf') < re(x), re(x) < -1.0),
                     And(Lt(1.0, re(x)), re(x) < Float('+inf'))))
    assert reduce_rational_inequalities([[Ne(x**2, 1.0)]], x, relational=True) == \
        And(cond, Or(And(Float('-inf') < re(x), Lt(re(x), -1.0)),
                     And(Lt(-1.0, re(x)), re(x) < 1.0), And(Lt(1.0, re(x)), re(x) < Float('+inf'))))


def test_reduce_rational_inequalities_real_relational():
    def OpenInterval(a, b):
        return Interval(a, b, True, True)
    def LeftOpenInterval(a, b):
        return Interval(a, b, True, False)
    def RightOpenInterval(a, b):
        return Interval(a, b, False, True)

    x = Symbol('x', real=True)

    assert reduce_rational_inequalities([[(x**2 + 3*x + 2)/(x**2 - 16) >= 0]], x, relational=False) == \
        Union(OpenInterval(-oo, -4), Interval(-2, -1), OpenInterval(4, oo))

    assert reduce_rational_inequalities([[((-2*x - 10)*(3 - x))/((x**2 + 5)*(x - 2)**2) < 0]], x, relational=False) == \
        Union(OpenInterval(-5, 2), OpenInterval(2, 3))

    assert reduce_rational_inequalities([[(x + 1)/(x - 5) <= 0]], x, relational=False) == \
        RightOpenInterval(-1, 5)

    assert reduce_rational_inequalities([[(x**2 + 4*x + 3)/(x - 1) > 0]], x, relational=False) == \
        Union(OpenInterval(-3, -1), OpenInterval(1, oo))

    assert reduce_rational_inequalities([[(x**2 - 16)/(x - 1)**2 < 0]], x, relational=False) == \
        Union(OpenInterval(-4, 1), OpenInterval(1, 4))

    assert reduce_rational_inequalities([[(3*x + 1)/(x + 4) >= 1]], x, relational=False) == \
        Union(OpenInterval(-oo, -4), RightOpenInterval(S(3)/2, oo))

    assert reduce_rational_inequalities([[(x - 8)/x <= 3 - x]], x, relational=False) == \
        Union(LeftOpenInterval(-oo, -2), LeftOpenInterval(0, 4))


def test_reduce_abs_inequalities():
    x = Symbol('x', real=True)

    assert reduce_inequalities(abs(x - 5) < 3) == And(Lt(2, x), Lt(x, 8))
    assert reduce_inequalities(
        abs(2*x + 3) >= 8) == Or(And(Le(S(5)/2, x), Lt(x, oo)), And(Le(x, -S(11)/2), Lt(-oo, x)))
    assert reduce_inequalities(abs(x - 4) + abs(
        3*x - 5) < 7) == And(Lt(S(1)/2, x), Lt(x, 4))
    assert reduce_inequalities(abs(x - 4) + abs(3*abs(x) - 5) < 7) == \
        Or(And(S(-2) < x, x < -1), And(S(1)/2 < x, x < 4))

    x = Symbol('x')
    raises(NotImplementedError, lambda: reduce_inequalities(abs(x - 5) < 3))


def test_reduce_inequalities_boolean():
    x = Symbol('x')

    assert reduce_inequalities(
        [Eq(x**2, 0), True]) == And(Eq(re(x), 0), Eq(im(x), 0))
    assert reduce_inequalities([Eq(x**2, 0), False]) is False


def test_reduce_inequalities_multivariate():
    x = Symbol('x')
    y = Symbol('y')

    assert reduce_inequalities([Ge(x**2, 1), Ge(y**2, 1)]) == \
        And(Eq(im(x), 0), Eq(im(y), 0), Or(And(Le(1, re(x)), Lt(re(x), oo)),
                                           And(Le(re(x), -1), Lt(-oo, re(x)))),
            Or(And(Le(1, re(y)), Lt(re(y), oo)), And(Le(re(y), -1), Lt(-oo, re(y)))))


def test_reduce_inequalities_errors():
    x = Symbol('x')
    y = Symbol('y')

    raises(NotImplementedError, lambda: reduce_inequalities(Ge(sin(x) + x, 1)))
    raises(NotImplementedError, lambda: reduce_inequalities(Ge(x**2*y + y, 1)))
    raises(NotImplementedError, lambda: reduce_inequalities(Ge(sqrt(2)*x, 1)))


def test_hacky_inequalities():
    x = Symbol('x')
    y = Symbol('y')

    assert reduce_inequalities(x + y < 1, symbols=[x]) == (x < 1 - y)
    assert reduce_inequalities(x + y >= 1, symbols=[x]) == (x >= 1 - y)


def test_issue_6343():
    x = Symbol('x', real=True)

    eq = -3*x**2/2 - 45*x/4 + S(33)/2 > 0
    assert reduce_inequalities(eq) == \
        And(x < -S(15)/4 + sqrt(401)/4, -sqrt(401)/4 - S(15)/4 < x)


def test_issue_8235():
    x = Symbol('x', real=True)
    assert reduce_inequalities(x**2 - 1 < 0) == \
        And(S(-1) < x, x < S(1))
    assert reduce_inequalities(x**2 - 1 <= 0) == \
        And(S(-1) <= x, x <= 1)
    assert reduce_inequalities(x**2 - 1 > 0) == \
        Or(And(-oo < x, x < -1), And(x < oo, S(1) < x))
    assert reduce_inequalities(x**2 - 1 >= 0) == \
        Or(And(-oo < x, x <= S(-1)), And(S(1) <= x, x < oo))

    eq = x**8 + x**2 - 9
    sol = solve(eq >= 0)
    known_sol = Or(And(-oo < x, RootOf(x**8 + x**2 - 9, 1) <= x, x < oo), \
            And(-oo < x, x < oo, x <= RootOf(x**8 + x**2 - 9, 0)))
    assert sol == known_sol


def test_issue_5526():
    x = Symbol('x')
    y = Symbol('y')

    assert reduce_inequalities(S(0) <= x + Integral(y**2, (y, 1, 3)) - 1, [x]) == \
        (-Integral(y**2, (y, 1, 3)) + 1 <= x)
    f = Function('f')
    e = Sum(f(x), (x, 1, 3))
    assert reduce_inequalities(S(0) <= x + e + y**2, [x]) == \
        (-y**2 - Sum(f(x), (x, 1, 3)) <= x)


def test_solve_univariate_inequality():
    x = Symbol('x', real=True)
    assert isolve(x**2 >= 4, x, relational=False) == Union(Interval(-oo, -2), Interval(2, oo))
    assert isolve(x**2 >= 4, x) == Or(And(Le(2, x), Lt(x, oo)), And(Le(x, -2), Lt(-oo, x)))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x, relational=False) == \
        Union(Interval(1, 2), Interval(3, oo))
    assert isolve((x - 1)*(x - 2)*(x - 3) >= 0, x) == \
        Or(And(Le(1, x), Le(x, 2)), And(Le(3, x), Lt(x, oo)))
    # issue 2785:
    assert isolve(x**3 - 2*x - 1 > 0, x, relational=False) == \
        Union(Interval(-1, -sqrt(5)/2 + S(1)/2, True, True),
              Interval(S(1)/2 + sqrt(5)/2, oo, True, True))
    # issue 2794:
    assert isolve(x**3 - x**2 + x - 1 > 0, x, relational=False) == \
        Interval(1, oo, True)
