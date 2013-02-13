from sympy.matrices.expressions.block import Block
from sympy.matrices.expressions import MatrixSymbol
from sympy.abc import a, b, c, d, k, l, m, n

X = MatrixSymbol('X', n, m)
Y = MatrixSymbol('Y', m, k)

def test_shape():
    B = Block(X, (a, b), (c, d))
    assert B.shape == (b - a, d - c)

def test_entry():
    B = Block(X, (a, b), (c, d))
    assert B[0,0] == X[a, c]
    assert B[k,l] == X[a+k, c+l]


