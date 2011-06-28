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




    pass
