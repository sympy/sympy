"""Tests for tools for solving inequalities and systems of inequalities. """

from sympy.solvers.inequalities import (
    solve_poly_inequality, solve_poly_inequalities,
)

from sympy import (
    Symbol, Interval, Eq, Ne, Lt, Le, Gt, Ge, Or,
    And, pi, oo, sqrt,
)

from sympy.utilities.pytest import raises
from sympy.abc import x

inf = oo.evalf()

def test_solve_poly_inequality():
    a = Symbol('a', real=True)

    assert solve_poly_inequality(Eq(a**2, 0)) == ([Interval(0, 0)], True, a)
    assert solve_poly_inequality(Le(a**2, 0)) == ([Interval(0, 0)], True, a)
    assert solve_poly_inequality(Lt(a**2, 0)) == ([], True, a)
    assert solve_poly_inequality(Ge(a**2, 0)) == ([Interval(-oo, oo)], True, a)
    assert solve_poly_inequality(Gt(a**2, 0)) == ([Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)], True, a)
    assert solve_poly_inequality(Ne(a**2, 0)) == ([Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)], True, a)

    assert solve_poly_inequality(Eq(a**2, 1)) == ([Interval(-1,-1), Interval(1, 1)], True, a)
    assert solve_poly_inequality(Le(a**2, 1)) == ([Interval(-1, 1)], True, a)
    assert solve_poly_inequality(Lt(a**2, 1)) == ([Interval(-1, 1, True, True)], True, a)
    assert solve_poly_inequality(Ge(a**2, 1)) == ([Interval(-oo, -1), Interval(1, oo)], True, a)
    assert solve_poly_inequality(Gt(a**2, 1)) == ([Interval(-oo, -1, right_open=True), Interval(1, oo, left_open=True)], True, a)
    assert solve_poly_inequality(Ne(a**2, 1)) == ([Interval(-oo, -1, right_open=True), Interval(-1, 1, True, True), Interval(1, oo, left_open=True)], True, a)

    assert solve_poly_inequality(Eq(a**2, 1.0)) == ([Interval(-1,-1), Interval(1, 1)], False, a)
    assert solve_poly_inequality(Le(a**2, 1.0)) == ([Interval(-1, 1)], False, a)
    assert solve_poly_inequality(Lt(a**2, 1.0)) == ([Interval(-1, 1, True, True)], False, a)
    assert solve_poly_inequality(Ge(a**2, 1.0)) == ([Interval(-oo, -1), Interval(1, oo)], False, a)
    assert solve_poly_inequality(Gt(a**2, 1.0)) == ([Interval(-oo, -1, right_open=True), Interval(1, oo, left_open=True)], False, a)
    assert solve_poly_inequality(Ne(a**2, 1.0)) == ([Interval(-oo, -1, right_open=True), Interval(-1, 1, True, True), Interval(1, oo, left_open=True)], False, a)

    raises(ValueError, "solve_poly_inequality(Gt(x**2, 0))")
    raises(ValueError, "solve_poly_inequality(Gt(a**2, pi))")

def test_solve_poly_inequalities():
    a = Symbol('a', real=True)
    b = Symbol('b', real=True)

    assert solve_poly_inequalities(Eq(a**2, 0)) == [Interval(0, 0)]
    assert solve_poly_inequalities(Le(a**2, 0)) == [Interval(0, 0)]
    assert solve_poly_inequalities(Lt(a**2, 0)) == []
    assert solve_poly_inequalities(Ge(a**2, 0)) == [Interval(-oo, oo)]
    assert solve_poly_inequalities(Gt(a**2, 0)) == [Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)]
    assert solve_poly_inequalities(Ne(a**2, 0)) == [Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)]

    assert solve_poly_inequalities(Eq(a**2, 1)) == [Interval(-1,-1), Interval(1, 1)]
    assert solve_poly_inequalities(Le(a**2, 1)) == [Interval(-1, 1)]
    assert solve_poly_inequalities(Lt(a**2, 1)) == [Interval(-1, 1, True, True)]
    assert solve_poly_inequalities(Ge(a**2, 1)) == [Interval(-oo, -1), Interval(1, oo)]
    assert solve_poly_inequalities(Gt(a**2, 1)) == [Interval(-oo, -1, right_open=True), Interval(1, oo, left_open=True)]
    assert solve_poly_inequalities(Ne(a**2, 1)) == [Interval(-oo, -1, right_open=True), Interval(-1, 1, True, True), Interval(1, oo, left_open=True)]

    assert solve_poly_inequalities(Eq(a**2, 1.0)) == [Interval(-1.0,-1.0), Interval(1.0, 1.0)]
    assert solve_poly_inequalities(Le(a**2, 1.0)) == [Interval(-1.0, 1.0)]
    assert solve_poly_inequalities(Lt(a**2, 1.0)) == [Interval(-1.0, 1.0, True, True)]
    assert solve_poly_inequalities(Ge(a**2, 1.0)) == [Interval(-inf, -1.0), Interval(1.0, inf)]
    assert solve_poly_inequalities(Gt(a**2, 1.0)) == [Interval(-inf, -1.0, right_open=True), Interval(1.0, inf, left_open=True)]
    assert solve_poly_inequalities(Ne(a**2, 1.0)) == [Interval(-inf, -1.0, right_open=True), Interval(-1.0, 1.0, True, True), Interval(1.0, inf, left_open=True)]

    assert solve_poly_inequalities(Eq(a**2, 0), relational=True) == Eq(a, 0)
    assert solve_poly_inequalities(Le(a**2, 0), relational=True) == Eq(a, 0)
    assert solve_poly_inequalities(Lt(a**2, 0), relational=True) == False
    assert solve_poly_inequalities(Ge(a**2, 0), relational=True) == True
    assert solve_poly_inequalities(Gt(a**2, 0), relational=True) == Or(Lt(a, 0), Lt(0, a))
    assert solve_poly_inequalities(Ne(a**2, 0), relational=True) == Or(Lt(a, 0), Lt(0, a))

    assert solve_poly_inequalities(Eq(a**2, 1), relational=True) == Or(Eq(a, -1), Eq(a, 1))
    assert solve_poly_inequalities(Le(a**2, 1), relational=True) == And(Le(-1, a), Le(a, 1))
    assert solve_poly_inequalities(Lt(a**2, 1), relational=True) == And(Lt(-1, a), Lt(a, 1))
    assert solve_poly_inequalities(Ge(a**2, 1), relational=True) == Or(Le(a, -1), Le(1, a))
    assert solve_poly_inequalities(Gt(a**2, 1), relational=True) == Or(Lt(a, -1), Lt(1, a))
    assert solve_poly_inequalities(Ne(a**2, 1), relational=True) == Or(Lt(a, -1), And(Lt(-1, a), Lt(a, 1)), Lt(1, a))

    assert solve_poly_inequalities(Eq(a**2, 1.0), relational=True) == Or(Eq(a, -1.0), Eq(a, 1.0))
    assert solve_poly_inequalities(Le(a**2, 1.0), relational=True) == And(Le(-1.0, a), Le(a, 1.0))
    assert solve_poly_inequalities(Lt(a**2, 1.0), relational=True) == And(Lt(-1.0, a), Lt(a, 1.0))
    assert solve_poly_inequalities(Ge(a**2, 1.0), relational=True) == Or(Le(a, -1.0), Le(1.0, a))
    assert solve_poly_inequalities(Gt(a**2, 1.0), relational=True) == Or(Lt(a, -1.0), Lt(1.0, a))
    assert solve_poly_inequalities(Ne(a**2, 1.0), relational=True) == Or(Lt(a, -1.0), And(Lt(-1.0, a), Lt(a, 1.0)), Lt(1.0, a))

    s = sqrt(2)

    assert solve_poly_inequalities([Lt(a**2 - 1, 0), Gt(a**2 - 1, 0)]) == []
    assert solve_poly_inequalities([Le(a**2 - 1, 0), Ge(a**2 - 1, 0)]) == [Interval(-1,-1), Interval(1, 1)]
    assert solve_poly_inequalities([Le(a**2 - 2, 0), Ge(a**2 - 1, 0)]) == [Interval(-s, -1, False, False), Interval(1, s, False, False)]
    assert solve_poly_inequalities([Le(a**2 - 2, 0), Gt(a**2 - 1, 0)]) == [Interval(-s, -1, False, True), Interval(1, s, True, False)]
    assert solve_poly_inequalities([Lt(a**2 - 2, 0), Ge(a**2 - 1, 0)]) == [Interval(-s, -1, True, False), Interval(1, s, False, True)]
    assert solve_poly_inequalities([Lt(a**2 - 2, 0), Gt(a**2 - 1, 0)]) == [Interval(-s, -1, True, True), Interval(1, s, True, True)]
    assert solve_poly_inequalities([Lt(a**2 - 2, 0), Ne(a**2 - 1, 0)]) == [Interval(-s, -1, True, True), Interval(-1, 1, True, True), Interval(1, s, True, True)]

    raises(ValueError, "solve_poly_inequalities([Ge(a**2, 0), Ge(b**2, 0)])")
