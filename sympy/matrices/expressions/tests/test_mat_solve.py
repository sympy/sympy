from sympy import Eq, Matrix, MatrixSymbol, Symbol, solve
from sympy.matrices.expressions.mat_solve import matrix_solve
from sympy.matrices.expressions.special import Identity
from sympy.testing.pytest import raises


def test_matrix_solve_left():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(A*X, B), X)

    assert result == {X: A.I*B}
    assert A*result[X] == B


def test_matrix_solve_right():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(X*A, B), X)

    assert result == {X: B*A.I}
    assert result[X]*A == B
    assert B*A.I != A.I*B


def test_matrix_solve_affine_left():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(A*X + A, B), X)

    assert result == {X: A.I*(B - A)}


def test_matrix_solve_affine_right():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(X*A + A, B), X)

    assert result == {X: (B - A)*A.I}


def test_matrix_solve_combined_left_coefficients():
    A = MatrixSymbol('A', 2, 2)
    C = MatrixSymbol('C', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(A*X + C*X, B), X)

    assert result == {X: (A + C).I*B}


def test_matrix_solve_combined_right_coefficients():
    A = MatrixSymbol('A', 2, 2)
    C = MatrixSymbol('C', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(X*A + X*C, B), X)

    assert result == {X: B*(A + C).I}


def test_matrix_solve_scalar_coefficients():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    left = matrix_solve(Eq(2*A*X, B), X)
    right = matrix_solve(Eq(2*X*A, B), X)

    assert left == {X: (2*A).I*B}
    assert right == {X: B*(2*A).I}


def test_matrix_solve_scalar_times_bare_unknown():
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(2*X, B), X)

    assert result == {X: Identity(2).I*B/2}


def test_matrix_solve_bare_unknown_with_left_term():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(X + A*X, B), X)

    assert result == {X: (Identity(2) + A).I*B}


def test_matrix_solve_bare_unknown_with_right_term():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(X + X*A, B), X)

    assert result == {X: B*(Identity(2) + A).I}


def test_matrix_solve_multiple_matrix_factors():
    A = MatrixSymbol('A', 2, 2)
    C = MatrixSymbol('C', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    left = matrix_solve(Eq(A*C*X, B), X)
    right = matrix_solve(Eq(X*A*C, B), X)

    assert left == {X: (A*C).I*B}
    assert right == {X: B*(A*C).I}


def test_matrix_solve_other_matrix_symbol_is_parameter():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    C = MatrixSymbol('C', 2, 2)
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)

    result = matrix_solve(Eq(A*X + B*Y, C), X)

    assert result == {X: A.I*(C - B*Y)}


def test_matrix_solve_rejects_mixed_sides():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    C = MatrixSymbol('C', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    raises(
        NotImplementedError,
        lambda: matrix_solve(Eq(A*X + X*C, B), X)
    )


def test_matrix_solve_rejects_two_sided_product():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    C = MatrixSymbol('C', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    raises(
        NotImplementedError,
        lambda: matrix_solve(Eq(A*X*B, C), X)
    )


def test_matrix_solve_rejects_nonlinear_equation():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    raises(
        NotImplementedError,
        lambda: matrix_solve(Eq(A*(X**2), B), X)
    )


def test_matrix_solve_rejects_nested_unknown():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    C = MatrixSymbol('C', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    raises(
        NotImplementedError,
        lambda: matrix_solve(Eq(A*(X + B)*C, C), X)
    )


def test_matrix_solve_rejects_missing_unknown():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    raises(ValueError, lambda: matrix_solve(Eq(A, B), X))


def test_matrix_solve_rejects_rectangular_coefficient():
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 3, 2)

    raises(
        NotImplementedError,
        lambda: matrix_solve(Eq(A*X, B), X)
    )


def test_matrix_solve_symbolic_dimensions():
    n = Symbol('n', integer=True)
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, n)
    X = MatrixSymbol('X', n, n)

    result = matrix_solve(Eq(A*X, B), X)

    assert result == {X: A.I*B}


def test_matrix_solve_concrete_right_multiplication():
    A = Matrix([[1, 1], [0, 1]])
    B = Matrix([[2, 0], [0, 3]])
    X = MatrixSymbol('X', 2, 2)

    result = matrix_solve(Eq(X*A, B), X)

    assert result[X].doit() == B*A.inv()
    assert result[X].doit()*A == B


def test_solve_matrixsymbol_left():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = solve(Eq(A*X, B), X)

    assert result == {X: A.I*B}


def test_solve_matrixsymbol_right():
    A = MatrixSymbol('A', 2, 2)
    B = MatrixSymbol('B', 2, 2)
    X = MatrixSymbol('X', 2, 2)

    result = solve(Eq(X*A, B), X)

    assert result == {X: B*A.I}
