"""Static type checks for MatrixBase and MatrixExpr arithmetic."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing_extensions import assert_type
    from sympy import Matrix, MatrixSymbol
    from sympy.matrices.dense import MutableDenseMatrix
    from sympy.matrices.expressions.matexpr import MatrixExpr

    M = Matrix([[1, 2], [3, 4]])
    X = MatrixSymbol("X", 2, 2)

    assert_type(M + M, MutableDenseMatrix)
    assert_type(M - M, MutableDenseMatrix)
    assert_type(M * M, MutableDenseMatrix)
    assert_type(M @ M, MutableDenseMatrix)
    assert_type(M * 2, MutableDenseMatrix)
    assert_type(2 * M, MutableDenseMatrix)

    assert_type(M + X, MatrixExpr)
    assert_type(X + M, MatrixExpr)

    assert_type(M - X, MatrixExpr)
    assert_type(X - M, MatrixExpr)

    assert_type(M * X, MatrixExpr)
    assert_type(X * M, MatrixExpr)

    assert_type(M @ X, MatrixExpr)
    assert_type(X @ M, MatrixExpr)

    assert_type(X + X, MatrixExpr)
    assert_type(X - X, MatrixExpr)
    assert_type(X * X, MatrixExpr)
    assert_type(X @ X, MatrixExpr)
    assert_type(-X, MatrixExpr)
    assert_type(2 * X, MatrixExpr)
    assert_type(X * 2, MatrixExpr)
    assert_type(X / 2, MatrixExpr)
