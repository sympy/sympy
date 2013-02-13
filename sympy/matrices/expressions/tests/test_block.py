from sympy.matrices.expressions.block import Block
from sympy.matrices.expressions import MatrixSymbol
from sympy.abc import a, b, c, d, k, l, m, n
from sympy.utilities.pytest import raises

X = MatrixSymbol('X', n, m)
Y = MatrixSymbol('Y', m, k)

def test_shape():
    B = Block(X, (a, b), (c, d))
    assert B.shape == (b - a, d - c)

def test_entry():
    B = Block(X, (a, b), (c, d))
    assert B[0,0] == X[a, c]
    assert B[k,l] == X[a+k, c+l]
    raises(IndexError, lambda : Block(X, 1, (2, 5))[1, 0])

def test_on_diag():
    assert not Block(X, (a, b), (c, d)).on_diag
    assert Block(X, (a, b), (a, b)).on_diag

def test_inputs():
    assert Block(X, 1, (2, 5)) == Block(X, (1, 2), (2, 5))
    assert Block(X, 1, (2, 5)).shape == (1, 3)
