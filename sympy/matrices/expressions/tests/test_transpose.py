from __future__ import annotations
from sympy.functions import adjoint, conjugate, transpose
from sympy.matrices.expressions import MatrixSymbol, Adjoint, trace, Transpose, Conjugate
from sympy.matrices import eye, Matrix
from sympy.assumptions.ask import Q
from sympy.assumptions.refine import refine
from sympy.core.singleton import S
from sympy.core.symbol import symbols

n, m, l, k, p = symbols('n m l k p', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)


def test_transpose():
    Sq = MatrixSymbol('Sq', n, n)

    assert transpose(A) == Transpose(A)
    assert Transpose(A).shape == (m, n)
    assert Transpose(A*B).shape == (l, n)
    assert transpose(Transpose(A)) == A
    assert isinstance(Transpose(Transpose(A)), Transpose)

    assert adjoint(Transpose(A)) == Conjugate(A)
    assert Adjoint(Transpose(A)).doit() == Conjugate(A)
    assert conjugate(Transpose(A)) == Adjoint(A)

    assert Transpose(eye(3)).doit() == eye(3)

    assert Transpose(S(5)).doit() == S(5)

    assert Transpose(Matrix([[1, 2], [3, 4]])).doit() == Matrix([[1, 3], [2, 4]])

    assert transpose(trace(Sq)) == trace(Sq)
    assert trace(Transpose(Sq)) == trace(Sq)

    assert Transpose(Sq)[0, 1] == Sq[1, 0]

    assert Transpose(A*B).doit() == Transpose(B) * Transpose(A)


def test_transpose_MatAdd_MatMul():
    # Issue 16807
    from sympy.functions.elementary.trigonometric import cos

    x = symbols('x')
    M = MatrixSymbol('M', 3, 3)
    N = MatrixSymbol('N', 3, 3)

    assert (N + (cos(x) * M)).T == cos(x)*M.T + N.T


def test_refine():
    assert refine(C.T, Q.symmetric(C)) == C


def test_transpose1x1():
    m = MatrixSymbol('m', 1, 1)
    assert m == refine(m.T)
    assert m == refine(m.T.T)

def test_issue_9817():
    from sympy.matrices.expressions import Identity
    v = MatrixSymbol('v', 3, 1)
    A = MatrixSymbol('A', 3, 3)
    x = Matrix([i + 1 for i in range(3)])
    X = Identity(3)
    quadratic = v.T * A * v
    subbed = quadratic.xreplace({v:x, A:X})
    assert subbed.as_explicit() == Matrix([[14]])


def test_transpose_derivative():
    from sympy.functions.elementary.trigonometric import sin, cos
    from sympy.matrices.expressions.special import MatrixUnit, ZeroMatrix

    x, y = symbols('x y')
    X = MatrixSymbol('X', 2, 3)
    Y = MatrixSymbol('Y', 3, 2)
    M = MatrixSymbol('M', 3, 3)

    # Differentiating Transpose(x * X) wrt scalar x
    expr1 = Transpose(x * X)
    assert expr1.shape == (3, 2)
    d1 = expr1.diff(x)
    assert d1 == X.T
    assert d1.shape == (3, 2)

    # Numeric/explicit check
    X_num = Matrix([[1, 2, 3], [4, 5, 6]])
    assert expr1.subs(X, X_num).as_explicit().diff(x) == d1.subs(X, X_num).as_explicit()

    # Differentiating Transpose(M) wrt matrix element M[i, j]
    d2 = (M.T).diff(M[0, 1])
    assert d2 == MatrixUnit(3, 3, 1, 0)
    assert (M.T).diff(M[1, 2]) == MatrixUnit(3, 3, 2, 1)

    M_num = Matrix(M)
    assert (M.T).diff(M[0, 1]).as_explicit() == (M_num.T).diff(M[0, 1])

    # Linear combination with trigonometric scalars
    expr3 = Transpose(sin(x) * X + cos(x) * Y.T)
    assert expr3.shape == (3, 2)
    d3 = expr3.diff(x)
    assert d3 == cos(x) * X.T - sin(x) * Y
    assert d3.shape == (3, 2)

    # Independent scalar yields ZeroMatrix with transposed shape
    d4 = Transpose(X).diff(y)
    assert d4 == ZeroMatrix(3, 2)
    assert d4.shape == (3, 2)
