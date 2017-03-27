from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices.expressions.diagonal import DiagonalMatrix, DiagonalOf
from sympy import Symbol, ask, Q

n = Symbol('n')
x = MatrixSymbol('x', n, 1)
X = MatrixSymbol('X', n, n)
D = DiagonalMatrix(x)
d = DiagonalOf(X)

def test_DiagonalMatrix():
    x = MatrixSymbol('x', n, 1)
    assert D.shape == (n, n)
    assert D[1, 2] == 0
    assert D[1, 1] == x[1, 0]
    
    i = Symbol('i')
    j = Symbol('j')
    x = MatrixSymbol('x', 3, 3)
    ij = DiagonalMatrix(x)[i, j]
    assert ij != 0
    assert ij.subs({i:0, j:0}) == x[0, 0]
    assert ij.subs({i:0, j:1}) == 0
    assert ij.subs({i:1, j:1}) == x[1, 1]

def test_DiagonalMatrix_Assumptions():
    assert ask(Q.diagonal(D))

def test_DiagonalOf():
    assert d.shape == (n, 1)
    assert d[2, 0] == X[2, 2]
