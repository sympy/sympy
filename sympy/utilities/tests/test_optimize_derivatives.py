"""
Tests for sympy.utilities.optimize_derivatives
"""

from sympy import symbols, exp, sin, cse
from sympy.utilities.optimize_derivatives import (
    generate_optimized_derivatives,
    count_operations,
)

x, y, z = symbols("x y z")


def _build_subs(replacements):
    """Helper to rebuild substitution dictionary from CSE output."""
    subs = {}
    for sym, val in replacements:
        subs[sym] = val.subs(subs)
    return subs


def test_returns_replacements_and_reduced():
    expr = (x + y) ** 2
    repl, reduced = generate_optimized_derivatives(expr, (x, y))

    assert isinstance(repl, list)
    assert isinstance(reduced, list)


def test_reduced_length_matches_layout():
    n = 2
    expr = (x + y) ** 3
    _, reduced = generate_optimized_derivatives(expr, (x, y))

    assert len(reduced) == 1 + n + n * n


def test_reduced_length_three_variables():
    n = 3
    expr = (x + y + z) ** 2
    _, reduced = generate_optimized_derivatives(expr, (x, y, z))

    assert len(reduced) == 1 + n + n * n


def test_gradient_values_are_correct():
    expr = x ** 2 + x * y + y ** 2
    repl, reduced = generate_optimized_derivatives(expr, (x, y))

    subs = _build_subs(repl)

    f_red = reduced[0].subs(subs)
    gx_red = reduced[1].subs(subs)
    gy_red = reduced[2].subs(subs)

    from sympy import simplify
    assert simplify(f_red - expr) == 0
    assert simplify(gx_red - expr.diff(x)) == 0
    assert simplify(gy_red - expr.diff(y)) == 0


def test_hessian_is_symmetric():
    expr = exp(x * y) + sin(x + y)
    repl, reduced = generate_optimized_derivatives(expr, (x, y))

    subs = _build_subs(repl)

    n = 2
    hess_reduced = reduced[1 + n:]

    from sympy import simplify
    h01 = hess_reduced[1].subs(subs)
    h10 = hess_reduced[2].subs(subs)

    assert simplify(h01 - h10) == 0



def test_joint_cse_finds_shared_subexpressions():
    expr = (x + y) ** 3

    repl, _ = generate_optimized_derivatives(expr, (x, y))

    assert len(repl) > 0

def test_joint_cse_runs_successfully():
    expr = (x + y) ** 4 + (x + y) ** 2
    repl, reduced = generate_optimized_derivatives(expr, (x, y))

    assert isinstance(repl, list)
    assert isinstance(reduced, list)

def test_joint_saves_ops_vs_separate():
    expr = exp(x + y) * (x ** 2 + y ** 2)
    vars_ = (x, y)

    repl_j, red_j = generate_optimized_derivatives(expr, vars_)
    ops_joint = count_operations(repl_j, red_j)

    f = expr
    grad = [expr.diff(v) for v in vars_]
    hess = [expr.diff(vi).diff(vj) for vi in vars_ for vj in vars_]

    def _ops(exprs):
        r, red = cse(exprs)
        return count_operations(r, red)

    ops_separate = _ops([f]) + _ops(grad) + _ops(hess)

    assert ops_joint <= ops_separate


def test_single_variable():
    expr = x ** 4 + x ** 2 + 1
    _, reduced = generate_optimized_derivatives(expr, (x,))

    assert len(reduced) == 3


def test_constant_expression():
    from sympy import Integer

    expr = Integer(5)
    _, reduced = generate_optimized_derivatives(expr, (x, y))

    assert len(reduced) == 1 + 2 + 4


def test_linear_expression():
    expr = 3 * x + 2 * y
    repl, reduced = generate_optimized_derivatives(expr, (x, y))

    subs = _build_subs(repl)
    hess_entries = [r.subs(subs) for r in reduced[3:]]

    assert all(e == 0 for e in hess_entries)