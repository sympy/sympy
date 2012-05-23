from sympy import (symbols, MatrixSymbol, Symbol, MatPow, BlockMatrix,
        Identity, ZeroMatrix, ImmutableMatrix, eye)
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
    assert x12 == Symbol(X.name+"_12")
    assert X[i,j] == Symbol(X.name+"_ij")

def test_add_index():
    assert (X+Y)[i,j] == X[i,j] + Y[i,j]

def test_mul_index():
    assert (A*y)[0,0] == A[0,0]*y[0,0] + A[0,1]*y[1,0]
    assert (A*B).as_mutable() == (A.as_mutable() * B.as_mutable())
    X = MatrixSymbol('X', n, m)
    Y = MatrixSymbol('Y', m, k)
    # Using str to avoid dealing with a Dummy variable
    assert str((X*Y)[4,2]) == "Sum(X(4, _k)*Y(_k, 2), (_k, 0, m - 1))"

def test_pow_index():
    Q = MatPow(A, 2)
    assert Q[0,0] == A[0,0]**2 + A[0,1]*A[1,0]

def test_transpose_index():
    assert X.T[i,j] == X[j,i]

def test_Identity_index():
    I = Identity(3)
    assert I[0,0] == I[1,1] == I[2,2] == 1
    assert I[1,0] == I[0,1] == I[2,1] == 0
    raises(IndexError, lambda: I[3,3])

def test_block_index():
    I = Identity(3)
    Z = ZeroMatrix(3,3)
    B = BlockMatrix([[I,I],[I,I]])
    e3 = ImmutableMatrix(eye(3))
    BB = BlockMatrix([[e3, e3], [e3, e3]])
    assert B[0,0] == B[3,0] == B[0,3] == B[3,3] == 1
    assert B[4,3] == B[5,1] == 0

    BB = BlockMatrix([[e3, e3], [e3, e3]])
    assert B.as_explicit() == BB.as_explicit()

    BI = BlockMatrix([[I, Z], [Z, I]])

    assert BI.as_explicit().equals(eye(6))

def test_slicing():
    raises(NotImplementedError, lambda: W[3,:])
    A.as_explicit()[0,:] # does not raise an error

def test_errors():
    raises(IndexError, lambda: Identity(2)[1,2,3,4,5])
    raises(IndexError, lambda: Identity(2)[[1,2,3,4,5]])
