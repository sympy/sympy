from sympy.functions import adjoint, conjugate, transpose
from sympy.matrices.expressions import MatrixSymbol, Adjoint, Trace, Transpose
from sympy.matrices import eye, Matrix
from sympy import symbols, S

n, m, l, k, p = symbols('n m l k p', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)


def test_transpose():
    Sq = MatrixSymbol('Sq', n, n)

    assert Transpose(A).shape == (m, n)
    assert Transpose(A*B).shape == (l, n)
    assert transpose(Transpose(A)) == A

    assert transpose(eye(3)) == eye(3)

    assert transpose(S(5)) == S(5)

    assert transpose(Matrix([[1, 2], [3, 4]])) == Matrix([[1, 3], [2, 4]])

    assert transpose(Trace(Sq)) == Trace(Sq)
