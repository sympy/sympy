from sympy.core import Lambda, S, symbols
from sympy.functions import adjoint, conjugate, transpose
from sympy.matrices import eye, Matrix, ShapeError
from sympy.matrices.expressions import (
    Adjoint, Identity, FunctionMatrix, MatrixExpr, MatrixSymbol, Trace,
    ZeroMatrix
)
from sympy.utilities.pytest import raises

n = symbols('n', integer=True)
A = MatrixSymbol('A', n, n)
B = MatrixSymbol('B', n, n)
C = MatrixSymbol('C', 3, 4)


def test_trace():
    assert isinstance(Trace(A), Trace)
    assert not isinstance(Trace(A), MatrixExpr)
    raises(ShapeError, lambda: Trace(C))
    assert Trace(eye(3)) == 3
    assert Trace(Matrix(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9])) == 15

    assert adjoint(Trace(A)) == Trace(Adjoint(A))
    assert conjugate(Trace(A)) == Trace(Adjoint(A))
    assert transpose(Trace(A)) == Trace(A)

    A / Trace(A)  # Make sure this is possible

    # Some easy simplifications
    assert Trace(Identity(5)) == 5
    assert Trace(ZeroMatrix(5, 5)) == 0
    assert Trace(2*A*B) == 2 * Trace(A*B)
    assert Trace(A.T) == Trace(A)

    i, j = symbols('i j')
    F = FunctionMatrix(3, 3, Lambda((i, j), i + j))
    assert Trace(F).doit() == (0 + 0) + (1 + 1) + (2 + 2)

    raises(TypeError, lambda: Trace(S.One))

    assert Trace(A).arg is A
