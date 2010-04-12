"""Tests for tools for solving inequalities and systems of inequalities. """

from sympy.solvers.inequalities import solve_poly_inequality

from sympy import Symbol, Interval, Eq, Ne, Lt, Le, Gt, Ge, pi, oo

from sympy.utilities.pytest import raises
from sympy.abc import x

def test_solve_poly_inequality():
    a = Symbol('a', real=True)

    assert solve_poly_inequality(Eq(a**2, 0)) == [Interval(0, 0)]
    assert solve_poly_inequality(Le(a**2, 0)) == [Interval(0, 0)]
    assert solve_poly_inequality(Lt(a**2, 0)) == []
    assert solve_poly_inequality(Ge(a**2, 0)) == [Interval(-oo, oo)]
    assert solve_poly_inequality(Gt(a**2, 0)) == [Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)]
    assert solve_poly_inequality(Ne(a**2, 0)) == [Interval(-oo, 0, right_open=True), Interval(0, oo, left_open=True)]

    raises(ValueError, "solve_poly_inequality(Gt(x**2, 0))")
    raises(ValueError, "solve_poly_inequality(Gt(a**2, pi))")
