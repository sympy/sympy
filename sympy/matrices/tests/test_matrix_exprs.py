from sympy.utilities.pytest import raises
from sympy import S, symbols, Symbol
from sympy.matrices import (eye, MatrixSymbol, Transpose, Inverse, ShapeError,
        MatMul)

def test_transpose():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)

    assert Transpose(A).shape == (m,n)
    assert Transpose(A*B).shape == (l,n)
    assert Transpose(Transpose(A)) == A

    assert Transpose(eye(3)) == eye(3)

def test_inverse():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    C = MatrixSymbol('C', n, n)

    raises(ShapeError, "Inverse(A)")
    assert Inverse(Inverse(C)) == C

    assert Inverse(eye(3)) == eye(3)

def test_shape():
    n, m, l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)
    assert A.shape == (n, m)
    assert (A*B).shape == (n, l)
    raises(ShapeError, 'B*A')

def test_matexpr():
    n, m, l = symbols('n m l', integer=True)
    x = Symbol('x')
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', m, l)

    assert (x*A).shape == A.shape
    assert (x*A).__class__ == MatMul
    assert 2*A - A - A == S.Zero

def test_BlockMatrix():
    n,m,l,k = symbols('n m l k', integer=True)
    A = MatrixSymbol('A', n, m)
    B = MatrixSymbol('B', n, k)
    C = MatrixSymbol('C', l, m)
    D = MatrixSymbol('D', l, k)
    X = BlockMatrix(Matrix([[A,B],[C,D]]))

    assert X.shape == (l+n, k+m)
    assert Transpose(X) == BlockMatrix(Matrix([[A.T, C.T], [B.T, D.T]]))
    assert Transpose(X).shape == X.shape[::-1]
    assert X.blockshape == (2,2)

    E = MatrixSymbol('E', m, 1)
    F = MatrixSymbol('F', k, 1)

    Y = BlockMatrix(Matrix([[E], [F]]))

    assert (X*Y).shape = (l+n, 1)
    assert block_collapse(X*Y)[0,0] == A*E + B*F
    assert block_collapse(X*Y)[1,0] == C*E + D*F
    assert Transpose(block_collapse(Transpose(X*Y))) == block_collapse(X*Y)

def test_BlockDiagMatrix():
    n,m,l = symbols('n m l', integer=True)
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', m, m)
    C = MatrixSymbol('C', l, l)

    X = BlockDiagMatrix(A,B,C)

    assert X[1,1] == B
    assert X.shape == (n+m+l, n+m+l)
    assert all(X[i,j].is_ZeroMatrix if i!=j else X[i,j] in [A,B,C]
            for i in range(3) for j in range(3))

    assert block_collapse(X.I * X).is_Identity

    assert block_collapse(X*X) == BlockDiagMatrix(A**2, B**2, C**2)



