from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices.expressions.diagonal import DiagonalMatrix
from sympy import Symbol, ask, Q

n = Symbol('n')
x = MatrixSymbol('x', n, 1)
D = DiagonalMatrix(x)

def test_DiagonalMatrix():
    assert D.shape == (n, n)
    assert D[1, 2] == 0
    assert D[1, 1] == x[1, 0]

def test_DiagonalMatrix_Assumptions():
    assert ask(Q.diagonal(D))
