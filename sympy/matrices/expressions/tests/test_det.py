from sympy.core import Lambda, S, symbols
from sympy.functions import adjoint, conjugate, transpose
from sympy.matrices import eye, Matrix, ShapeError
from sympy.matrices.expressions import (
    Adjoint, Identity, FunctionMatrix, MatrixExpr, MatrixSymbol, Det,
    ZeroMatrix
)
from sympy.utilities.pytest import raises

n = symbols('n', integer=True)
A = MatrixSymbol('A', n, n)
B = MatrixSymbol('B', n, n)
C = MatrixSymbol('C', 3, 4)


def test_det():
    assert isinstance(Det(A), Det)
    assert not isinstance(Det(A), MatrixExpr)
    raises(ShapeError, lambda: Det(C))
    assert Det(eye(3)) == 1
    assert Det(Matrix(3, 3, [1, 3, 2, 4, 1, 3, 2, 5, 2])) == 17
    A / Det(A)  # Make sure this is possible

    raises(TypeError, lambda: Det(S.One))

    assert Det(A).arg is A
