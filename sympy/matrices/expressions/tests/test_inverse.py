from sympy.core import symbols, S
from sympy.matrices.expressions import MatrixSymbol, Inverse, MatPow, ZeroMatrix, OneMatrix
from sympy.matrices.common import NonSquareMatrixError, NonInvertibleMatrixError
from sympy.matrices import eye, Identity
from sympy.testing.pytest import raises
from sympy import refine, Q

n, m, l = symbols('n m l', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)
D = MatrixSymbol('D', n, n)
E = MatrixSymbol('E', m, n)


def test_inverse():
    assert Inverse(C).args == (C, S.NegativeOne)
    assert Inverse(C).shape == (n, n)
    assert Inverse(A*E).shape == (n, n)
    assert Inverse(E*A).shape == (m, m)
    assert Inverse(C).inverse() == C
    assert Inverse(Inverse(C)).doit() == C
    assert isinstance(Inverse(Inverse(C)), Inverse)

    assert Inverse(*Inverse(E*A).args) == Inverse(E*A)

    assert C.inverse().inverse() == C

    assert C.inverse()*C == Identity(C.rows)

    assert Identity(n).inverse() == Identity(n)
    assert (3*Identity(n)).inverse() == Identity(n)/3

    # Simplifies Muls if possible (i.e. submatrices are square)
    assert (C*D).inverse() == D.I*C.I
    # But still works when not possible
    assert isinstance((A*E).inverse(), Inverse)
    assert Inverse(C*D).doit(inv_expand=False) == Inverse(C*D)

    assert Inverse(eye(3)).doit() == eye(3)
    assert Inverse(eye(3)).doit(deep=False) == eye(3)

    assert OneMatrix(1, 1).I == Identity(1)
    assert isinstance(OneMatrix(n, n).I, Inverse)

def test_inverse_non_invertible():
    raises(NonSquareMatrixError, lambda: Inverse(A))
    raises(NonSquareMatrixError, lambda: Inverse(A*B))
    raises(NonSquareMatrixError, lambda: ZeroMatrix(n, m).I)
    raises(NonInvertibleMatrixError, lambda: ZeroMatrix(n, n).I)
    raises(NonSquareMatrixError, lambda: OneMatrix(n, m).I)
    raises(NonInvertibleMatrixError, lambda: OneMatrix(2, 2).I)

def test_refine():
    assert refine(C.I, Q.orthogonal(C)) == C.T


def test_inverse_matpow_canonicalization():
    A = MatrixSymbol('A', 3, 3)
    assert Inverse(MatPow(A, 3)).doit() == MatPow(Inverse(A), 3).doit()
