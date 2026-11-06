"""Tests for the ``sympy.simplify._cse_diff.py`` module."""
from __future__ import annotations

import pytest

from sympy.core.symbol import (Symbol, symbols)
from sympy.core.numbers import Integer
from sympy.core.function import Function
from sympy.core import Derivative
from sympy.functions.elementary.exponential import exp, log
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.physics.mechanics import dynamicsymbols
from sympy.simplify._cse_diff import (SparseJacobianIR,
                                      _forward_jacobian,
                                      _remove_cse_from_derivative,
                                      _forward_jacobian_cse,
                                      _forward_jacobian_norm_in_cse_out,
                                      _forward_jacobian_sparse_cse,
                                      _get_dependencies,
                                      _propagate_support)
from sympy.simplify.simplify import simplify
from sympy.matrices import Matrix, eye

from sympy.testing.pytest import raises
from sympy.functions.elementary.trigonometric import (cos, sin, tan)
from sympy.simplify.trigsimp import trigsimp

from sympy import cse


w = Symbol('w')
x = Symbol('x')
y = Symbol('y')
z = Symbol('z')

q1, q2, q3 = dynamicsymbols('q1 q2 q3')

# Define the custom functions
k = Function('k')(x, y)
f = Function('f')(k, z)

zero = Integer(0)
one = Integer(1)
two = Integer(2)
neg_one = Integer(-1)


def _backsubstitute_matrix(replacements, matrix):
    sub_rep = dict(replacements)
    for sym, expr in replacements:
        sub_rep[sym] = expr.xreplace(sub_rep)
    return matrix.xreplace(sub_rep)


@pytest.mark.parametrize(
    'expr, wrt',
    [
        ([zero], [x]),
        ([one], [x]),
        ([two], [x]),
        ([neg_one], [x]),
        ([x], [x]),
        ([y], [x]),
        ([x + y], [x]),
        ([x*y], [x]),
        ([x**2], [x]),
        ([x**y], [x]),
        ([exp(x)], [x]),
        ([sin(x)], [x]),
        ([tan(x)], [x]),
        ([zero, one, x, y, x*y, x + y], [x, y]),
        ([((x/y) + sin(x/y) - exp(y))*((x/y) - exp(y))], [x, y]),
        ([w*tan(y*z)/(x - tan(y*z)), w*x*tan(y*z)/(x - tan(y*z))], [w, x, y, z]),
        ([q1**2 + q2, q2**2 + q3, q3**2 + q1], [q1, q2, q3]),
        ([f + Derivative(f, x) + k + 2*x], [x])
    ]
)


def test_forward_jacobian(expr, wrt):
    expr = ImmutableDenseMatrix([expr]).T
    wrt = ImmutableDenseMatrix([wrt]).T
    jacobian = _forward_jacobian(expr, wrt)
    zeros = ImmutableDenseMatrix.zeros(*jacobian.shape)
    assert simplify(jacobian - expr.jacobian(wrt)) == zeros


def test_process_cse():
    x, y, z = symbols('x y z')
    f = Function('f')
    k = Function('k')
    expr = Matrix([f(k(x,y), z) + Derivative(f(k(x,y), z), x) + k(x,y) + 2*x])
    repl, reduced = cse(expr)
    p_repl, p_reduced = _remove_cse_from_derivative(repl, reduced)

    x0 = symbols('x0')
    x1 = symbols('x1')

    expected_output = (
        [(x0, k(x, y)), (x1, f(x0, z))],
        [Matrix([2 * x + x0 + x1 + Derivative(f(k(x, y), z), x)])]
    )

    assert p_repl == expected_output[0], f"Expected {expected_output[0]}, but got {p_repl}"
    assert p_reduced == expected_output[1], f"Expected {expected_output[1]}, but got {p_reduced}"


def test_get_dependencies_mixed_wrt_types():
    t = symbols('t')
    expr = q1 + q2.diff(t) + x
    wrt = [x, q1, q2.diff(t)]

    assert _get_dependencies(expr, wrt) == {0, 1, 2}


def test_io_matrix_type():
    x, y, z = symbols('x y z')
    expr = ImmutableDenseMatrix([
        x * y + y * z + x * y * z,
        x ** 2 + y ** 2 + z ** 2,
        x * y + x * z + y * z
    ])
    wrt = ImmutableDenseMatrix([x, y, z])

    replacements, reduced_expr = cse(expr)

    # Test _forward_jacobian_core
    replacements_core, jacobian_core, precomputed_fs_core = _forward_jacobian_cse(replacements, reduced_expr, wrt)
    assert isinstance(jacobian_core[0], type(reduced_expr[0])), "Jacobian should be a Matrix of the same type as the input"

    # Test _forward_jacobian_norm_in_dag_out
    replacements_norm, jacobian_norm, precomputed_fs_norm = _forward_jacobian_norm_in_cse_out(
        expr, wrt)
    assert isinstance(jacobian_norm[0], type(reduced_expr[0])), "Jacobian should be a Matrix of the same type as the input"

    # Test _forward_jacobian
    jacobian = _forward_jacobian(expr, wrt)
    assert isinstance(jacobian, type(expr)), "Jacobian should be a Matrix of the same type as the input"


def test_forward_jacobian_sparse_cse_matches_dense_reduced_output():
    expr = Matrix([x*y + y*z + x*y*z, x**2 + y**2 + z**2, x*y + x*z + y*z])
    wrt = Matrix([x, y, z])

    replacements, reduced_expr = cse(expr)

    dense_replacements, dense_jacobian, _ = _forward_jacobian_cse(replacements, reduced_expr, wrt)
    sparse_replacements, sparse_jacobian, _ = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt)

    assert sparse_replacements == dense_replacements
    assert sparse_jacobian == dense_jacobian


def test_forward_jacobian_sparse_cse_handles_dynamicsymbols():
    expr = Matrix([q1**2 + q2, q2**2 + q3, q3**2 + q1])
    wrt = Matrix([q1, q2, q3])

    replacements, reduced_expr = cse(expr)
    _, sparse_jacobian, _ = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt)

    expanded = _backsubstitute_matrix(replacements, sparse_jacobian[0])
    assert simplify(expanded - expr.jacobian(wrt)) == Matrix.zeros(3, 3)


def test_sparse_cse_matches_public_forward_jacobian_for_basic_case():
    expr = Matrix([sin(x + y), exp(x*y), log(x + z)])
    wrt = Matrix([x, y, z])

    replacements, reduced_expr = cse(expr)
    _, sparse_jacobian, _ = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt)

    expanded = _backsubstitute_matrix(replacements, sparse_jacobian[0])
    assert simplify(expanded - _forward_jacobian(expr, wrt)) == Matrix.zeros(3, 3)


def test_sparse_cse_derivative_in_wrt_is_treated_as_leaf():
    t = symbols('t')
    q1_local, q2_local = dynamicsymbols('q1 q2')
    expr = Matrix([Derivative(q1_local, t) + q1_local,
                   q2_local + Derivative(q1_local, t)**2])
    wrt = Matrix([q1_local, Derivative(q1_local, t)])

    replacements, reduced_expr = cse(expr)
    _, sparse_jacobian, _ = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt)

    expanded = _backsubstitute_matrix(replacements, sparse_jacobian[0])
    assert simplify(expanded - expr.jacobian(wrt)) == Matrix.zeros(2, 2)


def test_sparse_cse_derivative_not_in_wrt_is_semantic_boundary():
    t = symbols('t')
    q1_local, q2_local = dynamicsymbols('q1 q2')
    expr = Matrix([Derivative(q1_local, t) + q1_local, q2_local + Derivative(q1_local, t)])
    wrt = Matrix([q1_local, q2_local])

    replacements, reduced_expr = cse(expr)
    _, sparse_jacobian, _ = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt)

    expanded = _backsubstitute_matrix(replacements, sparse_jacobian[0])
    assert simplify(expanded - expr.jacobian(wrt)) == Matrix.zeros(2, 2)


def test_sparse_cse_unsupported_function_preserves_jacobian():
    f_local = Function('f')
    expr = Matrix([f_local(x, y) + x, f_local(x, y)*z])
    wrt = Matrix([x, y, z])

    replacements, reduced_expr = cse(expr)
    ir = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt, return_mode="ir")

    expanded = _backsubstitute_matrix(replacements, ir.to_matrix(reduced_expr[0].__class__))
    assert simplify(expanded - expr.jacobian(wrt)) == Matrix.zeros(2, 3)


def test_sparse_cse_edge_case_zero_jacobian():
    expr = Matrix([1, z])
    wrt = Matrix([x, y])

    replacements, reduced_expr = cse(expr)
    ir = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt, return_mode="ir")

    assert ir.to_matrix(reduced_expr[0].__class__) == Matrix.zeros(2, 2)


def test_sparse_jacobian_ir_roundtrip_helpers():
    ir = SparseJacobianIR(
        shape=(2, 3),
        rows=[0, 0, 1],
        cols=[0, 2, 1],
        vals=[x, y, z],
        row_slices=[(0, 2), (2, 3)],
        wrt=[x, y, z],
        intermediates=[],
    )

    assert ir.to_matrix(Matrix) == Matrix([[x, 0, y], [0, z, 0]])
    assert ir.to_coo_tuple() == ([], (2, 3), [0, 0, 1], [0, 2, 1], [x, y, z])
    assert ir.row_vals(0) == [(0, x), (2, y)]
    assert ir.row_vals(1) == [(1, z)]


def test_forward_jacobian_sparse_cse_return_mode_ir_matches_matrix_mode():
    expr = Matrix([x*y + z, x + y])
    wrt = Matrix([x, y, z])
    replacements, reduced_expr = cse(expr)

    ir = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt, return_mode="ir")
    matrix_result = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt)

    assert ir.to_matrix(reduced_expr[0].__class__) == matrix_result[1][0]
    assert ir.intermediates == matrix_result[0]


def test_forward_jacobian_sparse_cse_rejects_unknown_return_mode():
    expr = Matrix([x + y])
    wrt = Matrix([x, y])
    replacements, reduced_expr = cse(expr)

    try:
        _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt, return_mode="bad")
    except ValueError as exc:
        assert "return_mode" in str(exc)
    else:
        raise AssertionError("Expected ValueError for unsupported return mode")


def test_propagate_support_treats_derivative_as_boundary():
    t = symbols('t')
    expr = Derivative(q1, t) + x

    assert _propagate_support(expr, {x: 0, q1: 1}, {}) == frozenset((0,))
    assert _propagate_support(expr, {x: 0, Derivative(q1, t): 1}, {}) == frozenset((0, 1))


def test_sparse_cse_structure_prepass_preserves_result():
    f_local = Function('f')
    expr = Matrix([f_local(x, y) + z, f_local(x, y)*z + x])
    wrt = Matrix([x, y, z])
    replacements, reduced_expr = cse(expr)

    default_result = _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt)
    prepass_result = _forward_jacobian_sparse_cse(
        replacements, reduced_expr, wrt, structure_prepass=True)

    default_expanded = _backsubstitute_matrix(replacements, default_result[1][0])
    prepass_expanded = _backsubstitute_matrix(replacements, prepass_result[1][0])
    assert simplify(default_expanded - prepass_expanded) == Matrix.zeros(2, 3)


def test_sparse_cse_structure_prepass_cache_miss_keeps_fallback_derivatives():
    expr = Function('f')(x, y) + z
    wrt = Matrix([x, y, z])
    replacements, reduced_expr = cse(Matrix([expr]))

    _, sparse_jacobian, _ = _forward_jacobian_sparse_cse(
        replacements, reduced_expr, wrt, structure_prepass=True)

    assert simplify(
        _backsubstitute_matrix(replacements, sparse_jacobian[0]) - Matrix([expr]).jacobian(wrt)
    ) == Matrix.zeros(1, 3)


def test_sparse_cse_row_value_cse_preserves_jacobian():
    expr = Matrix([sin(x + y) + exp(x + y)])
    wrt = Matrix([x, y])
    replacements, reduced_expr = cse(expr)

    row_result = _forward_jacobian_sparse_cse(
        replacements, reduced_expr, wrt, value_cse="row")
    expanded = _backsubstitute_matrix(row_result[0], row_result[1][0])

    assert simplify(expanded - expr.jacobian(wrt)) == Matrix.zeros(1, 2)
    assert len(row_result[0]) >= len(replacements)


def test_sparse_cse_global_value_cse_preserves_jacobian():
    expr = Matrix([
        sin(x + y) + exp(x + y),
        (sin(x + y) + exp(x + y))**2,
    ])
    wrt = Matrix([x, y])
    replacements, reduced_expr = cse(expr)

    global_result = _forward_jacobian_sparse_cse(
        replacements, reduced_expr, wrt, value_cse="global")
    expanded = _backsubstitute_matrix(global_result[0], global_result[1][0])

    assert simplify(expanded - expr.jacobian(wrt)) == Matrix.zeros(2, 2)
    assert len(global_result[0]) >= len(replacements)


def test_sparse_cse_value_cse_return_mode_ir_keeps_intermediates():
    expr = Matrix([sin(x + y) + exp(x + y)])
    wrt = Matrix([x, y])
    replacements, reduced_expr = cse(expr)

    ir = _forward_jacobian_sparse_cse(
        replacements, reduced_expr, wrt, return_mode="ir", value_cse="row")

    assert ir.intermediates
    assert simplify(
        _backsubstitute_matrix(ir.intermediates, ir.to_matrix(Matrix)) - expr.jacobian(wrt)
    ) == Matrix.zeros(1, 2)


def test_sparse_cse_rejects_unknown_value_cse_mode():
    expr = Matrix([x + y])
    wrt = Matrix([x, y])
    replacements, reduced_expr = cse(expr)

    try:
        _forward_jacobian_sparse_cse(replacements, reduced_expr, wrt, value_cse="bad")
    except ValueError as exc:
        assert "value_cse" in str(exc)
    else:
        raise AssertionError("Expected ValueError for unsupported value_cse")


def test_jacobian_hessian():
    L = Matrix(1, 2, [x**2*y, 2*y**2 + x*y])
    syms = [x, y]
    assert _forward_jacobian(L, syms) == Matrix([[2*x*y, x**2], [y, 4*y + x]])

    L = Matrix(1, 2, [x, x**2*y**3])
    assert _forward_jacobian(L, syms) == Matrix([[1, 0], [2*x*y**3, x**2*3*y**2]])


def test_jacobian_metrics():
    rho, phi = symbols("rho,phi")
    X = Matrix([rho * cos(phi), rho * sin(phi)])
    Y = Matrix([rho, phi])
    J = _forward_jacobian(X, Y)
    assert J == X.jacobian(Y.T)
    assert J == (X.T).jacobian(Y)
    assert J == (X.T).jacobian(Y.T)
    g = J.T * eye(J.shape[0]) * J
    g = g.applyfunc(trigsimp)
    assert g == Matrix([[1, 0], [0, rho ** 2]])


def test_jacobian2():
    rho, phi = symbols("rho,phi")
    X = Matrix([rho * cos(phi), rho * sin(phi), rho ** 2])
    Y = Matrix([rho, phi])
    J = Matrix([
        [cos(phi), -rho * sin(phi)],
        [sin(phi), rho * cos(phi)],
        [2 * rho, 0],
    ])
    assert _forward_jacobian(X, Y) == J


def test_issue_4564():
    X = Matrix([exp(x + y + z), exp(x + y + z), exp(x + y + z)])
    Y = Matrix([x, y, z])
    for i in range(1, 3):
        for j in range(1, 3):
            X_slice = X[:i, :]
            Y_slice = Y[:j, :]
            J = _forward_jacobian(X_slice, Y_slice)
            assert J.rows == i
            assert J.cols == j
            for k in range(j):
                assert J[:, k] == X_slice


def test_nonvectorJacobian():
    X = Matrix([[exp(x + y + z), exp(x + y + z)],
                [exp(x + y + z), exp(x + y + z)]])
    raises(TypeError, lambda: _forward_jacobian(X, Matrix([x, y, z])))
    X = X[0, :]
    Y = Matrix([[x, y], [x, z]])
    raises(TypeError, lambda: _forward_jacobian(X, Y))
    raises(TypeError, lambda: _forward_jacobian(X, Matrix([[x, y], [x, z]])))
