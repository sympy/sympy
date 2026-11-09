from __future__ import annotations
from sympy.core import symbols, S
from sympy.matrices.expressions import MatrixSymbol, Inverse, MatPow, ZeroMatrix, OneMatrix
from sympy.matrices.exceptions import NonInvertibleMatrixError, NonSquareMatrixError
from sympy.matrices import eye, Identity
from sympy.testing.pytest import raises
from sympy.assumptions.ask import Q
from sympy.assumptions.refine import refine

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
    raises(NonInvertibleMatrixError, lambda: ZeroMatrix(n, n).I)
    raises(NonInvertibleMatrixError, lambda: OneMatrix(2, 2).I)

def test_refine():
    assert refine(C.I, Q.orthogonal(C)) == C.T


def test_inverse_matpow_canonicalization():
    A = MatrixSymbol('A', 3, 3)
    assert Inverse(MatPow(A, 3)).doit() == MatPow(Inverse(A), 3).doit()


def test_nonsquare_error():
    A = MatrixSymbol('A', 3, 4)
    raises(NonSquareMatrixError, lambda: Inverse(A))


def test_adjoint_trnaspose_conjugate():
    A = MatrixSymbol('A', n, n)
    assert A.transpose().inverse() == A.inverse().transpose()
    assert A.conjugate().inverse() == A.inverse().conjugate()
    assert A.adjoint().inverse() == A.inverse().adjoint()


def test_inverse_derivative():
    from sympy.core.symbol import Symbol
    from sympy.matrices.expressions.special import MatrixUnit
    from sympy.tensor.array.expressions.array_expressions import PermuteDims
    x = Symbol('x')
    X = MatrixSymbol('X', 3, 3)
    Z = MatrixSymbol('Z', 3, 3)

    # Scalar differentiation of matrix inverse: d(M^-1)/dx = - M^-1 (dM/dx) M^-1
    inv_expr = Inverse(x*X + Z)
    assert inv_expr.diff(x) == -(x*X + Z)**(-1) * X * (x*X + Z)**(-1)

    # Constant matrix derivative is ZeroMatrix
    assert Inverse(X).diff(x) == ZeroMatrix(3, 3)

    # Non-linear scalar argument with chain rule
    assert Inverse(x**2 * X + Z).diff(x) == (-2*x)*(x**2*X + Z)**(-1) * X * (x**2*X + Z)**(-1)

    # Matrix element differentiation
    assert Inverse(X).diff(X[0, 1]) == -X**(-1) * MatrixUnit(3, 3, 0, 1) * X**(-1)

    # Matrix-by-matrix differentiation preserves tensor array engine
    assert isinstance(Inverse(X).diff(X), PermuteDims)
