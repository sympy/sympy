from sympy.core import Lambda, S, symbols
from sympy.functions import adjoint, conjugate, transpose
from sympy.matrices import eye, Matrix, ShapeError
from sympy.matrices.expressions import (
    Adjoint, Identity, FunctionMatrix, MatrixExpr, MatrixSymbol, Determinant,
    det, ZeroMatrix
)
from sympy.utilities.pytest import raises

n = symbols('n', integer=True)
A = MatrixSymbol('A', n, n)
B = MatrixSymbol('B', n, n)
C = MatrixSymbol('C', 3, 4)


def test_det():
    assert isinstance(Determinant(A), Determinant)
    assert not isinstance(Determinant(A), MatrixExpr)
    raises(ShapeError, lambda: Determinant(C))
    assert det(eye(3)) == 1
    assert det(Matrix(3, 3, [1, 3, 2, 4, 1, 3, 2, 5, 2])) == 17
    A / det(A)  # Make sure this is possible

    raises(TypeError, lambda: Determinant(S.One))

    assert Determinant(A).arg is A
