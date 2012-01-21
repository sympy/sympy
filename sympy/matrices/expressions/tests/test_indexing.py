from sympy import *
from sympy.utilities.pytest import raises, XFAIL

k,l,m,n = symbols('k l m n', integer=True)
i,j = symbols('i j', integer=True)

W = MatrixSymbol('W', k, l)
X = MatrixSymbol('X', l, m)
Y = MatrixSymbol('Y', l, m)
Z = MatrixSymbol('Z', m, n)

A = MatrixSymbol('A', 2,2)
B = MatrixSymbol('B', 2,2)
x = MatrixSymbol('x', 1,2)
y = MatrixSymbol('x', 2,1)

def test_symbolic_indexing():
    x12 = X[1,2]
    assert x12 == Symbol(X.name)(1,2)
    assert X[i,j] == Symbol(X.name)(i,j)

def test_add_index():
    assert (X+Y)[i,j] == X[i,j] + Y[i,j]

def test_mul_index():
    assert (A*y)[0,0] == A[0,0]*y[0,0] + A[0,1]*y[1,0]
    assert (A*B).to_explicit() == A.to_explicit() * B.to_explicit()

def test_Identity_index():
    I = Identity(3)
    assert I[0,0] == I[1,1] == I[2,2] == 1
    assert I[1,0] == I[0,1] == I[2,1] == 0
    raises(ValueError, "I[3,3]")

def test_block_index():
    I = Identity(3)
    Z = ZeroMatrix(3,3)
    B = BlockMatrix([[I,I],[I,I]])
    BB = BlockMatrix([[eye(3), eye(3)], [eye(3), eye(3)]])
    assert B[0,0] == B[3,0] == B[0,3] == B[3,3] == 1
    assert B[4,3] == B[5,1] == 0

    BB = BlockMatrix([[eye(3), eye(3)], [eye(3), eye(3)]])
    assert B.to_explicit() == BB.to_explicit()

    BI = BlockMatrix([[I, Z], [Z, I]])

    assert BI.to_explicit() == eye(6)

