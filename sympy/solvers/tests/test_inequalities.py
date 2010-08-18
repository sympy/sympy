"""Tests for tools for solving inequalities and systems of inequalities. """

from sympy.solvers.inequalities import solve_poly_inequalities

from sympy import (
    Symbol, Interval, Eq, Ne, Lt, Le, Gt, Ge, Or, And, pi, oo, sqrt,
    raises, Q, Assume, global_assumptions, re, im, sin, DomainError)

from sympy.utilities.pytest import raises
from sympy.abc import x, y

inf = oo.evalf()

x_assume = Assume(x, Q.real)
y_assume = Assume(y, Q.real)

def test_solve_poly_inequalities_real_interval():
    global_assumptions.add(x_assume)
    global_assumptions.add(y_assume)

    assert solve_poly_inequalities(Eq(x**2, 0), relational=False) == [Interval(0, 0)]
    assert solve_poly_inequalities(Le(x**2, 0), relational=False) == [Interval(0, 0)]
    assert solve_poly_inequalities(Lt(x**2, 0), relational=False) == []
    assert solve_poly_inequalities(Ge(x**2, 0), relational=False) == [Interval(-oo, oo)]
    assert solve_poly_inequalities(Gt(x**2, 0), relational=False) == [Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)]
    assert solve_poly_inequalities(Ne(x**2, 0), relational=False) == [Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)]

    assert solve_poly_inequalities(Eq(x**2, 1), relational=False) == [Interval(-1,-1), Interval(1, 1)]
    assert solve_poly_inequalities(Le(x**2, 1), relational=False) == [Interval(-1, 1)]
    assert solve_poly_inequalities(Lt(x**2, 1), relational=False) == [Interval(-1, 1, True, True)]
    assert solve_poly_inequalities(Ge(x**2, 1), relational=False) == [Interval(-oo, -1), Interval(1, oo)]
    assert solve_poly_inequalities(Gt(x**2, 1), relational=False) == [Interval(-oo, -1, right_open=True), Interval(1, oo, left_open=True)]
    assert solve_poly_inequalities(Ne(x**2, 1), relational=False) == [Interval(-oo, -1, right_open=True), Interval(-1, 1, True, True), Interval(1, oo, left_open=True)]

    assert solve_poly_inequalities(Eq(x**2, 1.0), relational=False) == [Interval(-1.0,-1.0), Interval(1.0, 1.0)]
    assert solve_poly_inequalities(Le(x**2, 1.0), relational=False) == [Interval(-1.0, 1.0)]
    assert solve_poly_inequalities(Lt(x**2, 1.0), relational=False) == [Interval(-1.0, 1.0, True, True)]
    assert solve_poly_inequalities(Ge(x**2, 1.0), relational=False) == [Interval(-inf, -1.0), Interval(1.0, inf)]
    assert solve_poly_inequalities(Gt(x**2, 1.0), relational=False) == [Interval(-inf, -1.0, right_open=True), Interval(1.0, inf, left_open=True)]
    assert solve_poly_inequalities(Ne(x**2, 1.0), relational=False) == [Interval(-inf, -1.0, right_open=True), Interval(-1.0, 1.0, True, True), Interval(1.0, inf, left_open=True)]

    s = sqrt(2)

    assert solve_poly_inequalities([Lt(x**2 - 1, 0), Gt(x**2 - 1, 0)], relational=False) == []
    assert solve_poly_inequalities([Le(x**2 - 1, 0), Ge(x**2 - 1, 0)], relational=False) == [Interval(-1,-1), Interval(1, 1)]
    assert solve_poly_inequalities([Le(x**2 - 2, 0), Ge(x**2 - 1, 0)], relational=False) == [Interval(-s, -1, False, False), Interval(1, s, False, False)]
    assert solve_poly_inequalities([Le(x**2 - 2, 0), Gt(x**2 - 1, 0)], relational=False) == [Interval(-s, -1, False, True), Interval(1, s, True, False)]
    assert solve_poly_inequalities([Lt(x**2 - 2, 0), Ge(x**2 - 1, 0)], relational=False) == [Interval(-s, -1, True, False), Interval(1, s, False, True)]
    assert solve_poly_inequalities([Lt(x**2 - 2, 0), Gt(x**2 - 1, 0)], relational=False) == [Interval(-s, -1, True, True), Interval(1, s, True, True)]
    assert solve_poly_inequalities([Lt(x**2 - 2, 0), Ne(x**2 - 1, 0)], relational=False) == [Interval(-s, -1, True, True), Interval(-1, 1, True, True), Interval(1, s, True, True)]

    global_assumptions.remove(x_assume)
    global_assumptions.remove(y_assume)

def test_solve_poly_inequalities_real_relational():
    global_assumptions.add(x_assume)
    global_assumptions.add(y_assume)

    assert solve_poly_inequalities(Eq(x**2, 0), relational=True) == Eq(x, 0)
    assert solve_poly_inequalities(Le(x**2, 0), relational=True) == Eq(x, 0)
    assert solve_poly_inequalities(Lt(x**2, 0), relational=True) == False
    assert solve_poly_inequalities(Ge(x**2, 0), relational=True) == True
    assert solve_poly_inequalities(Gt(x**2, 0), relational=True) == Or(Lt(x, 0), Lt(0, x))
    assert solve_poly_inequalities(Ne(x**2, 0), relational=True) == Or(Lt(x, 0), Lt(0, x))

    assert solve_poly_inequalities(Eq(x**2, 1), relational=True) == Or(Eq(x, -1), Eq(x, 1))
    assert solve_poly_inequalities(Le(x**2, 1), relational=True) == And(Le(-1, x), Le(x, 1))
    assert solve_poly_inequalities(Lt(x**2, 1), relational=True) == And(Lt(-1, x), Lt(x, 1))
    assert solve_poly_inequalities(Ge(x**2, 1), relational=True) == Or(Le(x, -1), Le(1, x))
    assert solve_poly_inequalities(Gt(x**2, 1), relational=True) == Or(Lt(x, -1), Lt(1, x))
    assert solve_poly_inequalities(Ne(x**2, 1), relational=True) == Or(Lt(x, -1), And(Lt(-1, x), Lt(x, 1)), Lt(1, x))

    assert solve_poly_inequalities(Eq(x**2, 1.0), relational=True) == Or(Eq(x, -1.0), Eq(x, 1.0))
    assert solve_poly_inequalities(Le(x**2, 1.0), relational=True) == And(Le(-1.0, x), Le(x, 1.0))
    assert solve_poly_inequalities(Lt(x**2, 1.0), relational=True) == And(Lt(-1.0, x), Lt(x, 1.0))
    assert solve_poly_inequalities(Ge(x**2, 1.0), relational=True) == Or(Le(x, -1.0), Le(1.0, x))
    assert solve_poly_inequalities(Gt(x**2, 1.0), relational=True) == Or(Lt(x, -1.0), Lt(1.0, x))
    assert solve_poly_inequalities(Ne(x**2, 1.0), relational=True) == Or(Lt(x, -1.0), And(Lt(-1.0, x), Lt(x, 1.0)), Lt(1.0, x))

    global_assumptions.remove(x_assume)
    global_assumptions.remove(y_assume)

def test_solve_poly_inequalities_complex():
    cond = Eq(im(x), 0)

    assert solve_poly_inequalities(Eq(x**2, 0), relational=True) == And(Eq(re(x), 0), cond)
    assert solve_poly_inequalities(Le(x**2, 0), relational=True) == And(Eq(re(x), 0), cond)
    assert solve_poly_inequalities(Lt(x**2, 0), relational=True) == False
    assert solve_poly_inequalities(Ge(x**2, 0), relational=True) == cond
    assert solve_poly_inequalities(Gt(x**2, 0), relational=True) == And(Or(Lt(re(x), 0), Lt(0, re(x))), cond)
    assert solve_poly_inequalities(Ne(x**2, 0), relational=True) == And(Or(Lt(re(x), 0), Lt(0, re(x))), cond)

    assert solve_poly_inequalities(Eq(x**2, 1), relational=True) == And(Or(Eq(re(x), -1), Eq(re(x), 1)), cond)
    assert solve_poly_inequalities(Le(x**2, 1), relational=True) == And(And(Le(-1, re(x)), Le(re(x), 1)), cond)
    assert solve_poly_inequalities(Lt(x**2, 1), relational=True) == And(And(Lt(-1, re(x)), Lt(re(x), 1)), cond)
    assert solve_poly_inequalities(Ge(x**2, 1), relational=True) == And(Or(Le(re(x), -1), Le(1, re(x))), cond)
    assert solve_poly_inequalities(Gt(x**2, 1), relational=True) == And(Or(Lt(re(x), -1), Lt(1, re(x))), cond)
    assert solve_poly_inequalities(Ne(x**2, 1), relational=True) == And(Or(Lt(re(x), -1), And(Lt(-1, re(x)), Lt(re(x), 1)), Lt(1, re(x))), cond)

    assert solve_poly_inequalities(Eq(x**2, 1.0), relational=True) == And(Or(Eq(re(x), -1.0), Eq(re(x), 1.0)), cond)
    assert solve_poly_inequalities(Le(x**2, 1.0), relational=True) == And(And(Le(-1.0, re(x)), Le(re(x), 1.0)), cond)
    assert solve_poly_inequalities(Lt(x**2, 1.0), relational=True) == And(And(Lt(-1.0, re(x)), Lt(re(x), 1.0)), cond)
    assert solve_poly_inequalities(Ge(x**2, 1.0), relational=True) == And(Or(Le(re(x), -1.0), Le(1.0, re(x))), cond)
    assert solve_poly_inequalities(Gt(x**2, 1.0), relational=True) == And(Or(Lt(re(x), -1.0), Lt(1.0, re(x))), cond)
    assert solve_poly_inequalities(Ne(x**2, 1.0), relational=True) == And(Or(Lt(re(x), -1.0), And(Lt(-1.0, re(x)), Lt(re(x), 1.0)), Lt(1.0, re(x))), cond)

def test_solve_poly_inequalities_multivariate():
    assert solve_poly_inequalities([Ge(x**2, 1), Ge(y**2, 1)]) == \
        And(And(Or(Le(re(x), -1), Le(1, re(x))), Eq(im(x), 0)),
            And(Or(Le(re(y), -1), Le(1, re(y))), Eq(im(y), 0)))

def test_solve_poly_inequalities_errors():
    raises(NotImplementedError, "solve_poly_inequalities(Ge(sin(x) + x, 1))")
    raises(NotImplementedError, "solve_poly_inequalities(Ge(x**2*y + y, 1))")

    raises(DomainError, "solve_poly_inequalities(Ge(sqrt(2)*x, 1))")
