from sympy.core import symbols
from sympy.utilities.pytest import raises

from sympy.matrices import (ShapeError, Matrix, MatrixSymbol, MatrixExpr,
    randMatrix)
from sympy.matrices.expressions import HadamardProduct, hadamard_product

n, m, k = symbols('n,m,k')
Z = MatrixSymbol('Z', n, n)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', n, m)
C = MatrixSymbol('C', m, k)
M = randMatrix(3)
N = randMatrix(3)
P = MatrixSymbol('P', 3, 3)
Q = MatrixSymbol('Q', 3, 3)

def test_HadamardProduct():
    assert HadamardProduct(A, B, A).shape == A.shape
    assert HadamardProduct(P, M, Q).shape == Q.shape

    raises(ShapeError, lambda: HadamardProduct(A, B.T))
    raises(TypeError,  lambda: HadamardProduct(A, n))
    raises(TypeError,  lambda: HadamardProduct(A, 1))

    assert HadamardProduct(A, 2*B, -A)[1, 1] == \
            -2 * A[1, 1] * B[1, 1] * A[1, 1]
    assert HadamardProduct(M, 3*P, N)[0, 1] == \
            3 * M[0, 1] * N[0, 1] * P[0, 1]

    mix = HadamardProduct(Z*A, B)*C
    assert mix.shape == (n, k)
    mix = HadamardProduct(M*N, P)*Q
    assert mix.shape == P.shape

    assert set(HadamardProduct(A, B, A).T.args) == set((A.T, A.T, B.T))

def test_HadamardProduct_isnt_commutative():
    assert HadamardProduct(A, B) != HadamardProduct(B, A)

def test_mixed_indexing():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    Z = MatrixSymbol('Z', 2, 2)

    assert (X*HadamardProduct(Y, Z))[0, 0] == \
            X[0, 0]*Y[0, 0]*Z[0, 0] + X[0, 1]*Y[1, 0]*Z[1, 0]
    assert (P*HadamardProduct(M, Q))[1, 1] == \
            P[1, 0]*M[0, 1]*Q[0, 1] + P[1, 1]*M[1, 1]*Q[1, 1] + \
            P[1, 2]*M[2, 1]*Q[2, 1]

def test_canonicalize():
    X = MatrixSymbol('X', 2, 2)
    expr = HadamardProduct(X, check=False)
    assert isinstance(expr, HadamardProduct)
    expr2 = expr.doit() # unpack is called
    assert isinstance(expr2, MatrixSymbol)

    expr = HadamardProduct(M, check=False)
    assert isinstance(expr, HadamardProduct)
    expr2 = expr.doit()
    assert isinstance(expr2, MatrixExpr)

    expr = HadamardProduct(M, P)
    assert isinstance(expr, HadamardProduct)
    expr2 = expr.doit()
    assert isinstance(expr2, MatrixExpr)

def test_hadamard():
    m, n, p = symbols('m, n, p', integer=True)
    A = MatrixSymbol('A', m, n)
    B = MatrixSymbol('B', m, n)
    C = MatrixSymbol('C', m, p)
    with raises(TypeError):
        hadamard_product()
    assert hadamard_product(A) == A
    assert isinstance(hadamard_product(A, B), HadamardProduct)
    assert hadamard_product(A, B).doit() == hadamard_product(A, B)
    with raises(ShapeError):
        hadamard_product(A, C)

    assert isinstance(hadamard_product(M, N, P), HadamardProduct)
    assert hadamard_product(M, N, P) == hadamard_product(P, M, N).doit()
