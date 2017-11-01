from sympy import (symbols, MatrixSymbol, MatPow, BlockMatrix, KroneckerDelta,
        Identity, ZeroMatrix, ImmutableMatrix, eye, Sum, MatMul)
from sympy.utilities.pytest import raises
from sympy.matrices.expressions.matexpr import MatrixExpr

k, l, m, n = symbols('k l m n', integer=True)
i, j = symbols('i j', integer=True)

W = MatrixSymbol('W', k, l)
X = MatrixSymbol('X', l, m)
Y = MatrixSymbol('Y', l, m)
Z = MatrixSymbol('Z', m, n)

A = MatrixSymbol('A', 2, 2)
B = MatrixSymbol('B', 2, 2)
x = MatrixSymbol('x', 1, 2)
y = MatrixSymbol('x', 2, 1)


def test_symbolic_indexing():
    x12 = X[1, 2]
    assert all(s in str(x12) for s in ['1', '2', X.name])
    # We don't care about the exact form of this.  We do want to make sure
    # that all of these features are present


def test_add_index():
    assert (X + Y)[i, j] == X[i, j] + Y[i, j]


def test_mul_index():
    assert (A*y)[0, 0] == A[0, 0]*y[0, 0] + A[0, 1]*y[1, 0]
    assert (A*B).as_mutable() == (A.as_mutable() * B.as_mutable())
    X = MatrixSymbol('X', n, m)
    Y = MatrixSymbol('Y', m, k)

    result = (X*Y)[4,2]
    expected = Sum(X[4, i]*Y[i, 2], (i, 0, m - 1))
    assert result.args[0].dummy_eq(expected.args[0], i)
    assert result.args[1][1:] == expected.args[1][1:]


def test_pow_index():
    Q = MatPow(A, 2)
    assert Q[0, 0] == A[0, 0]**2 + A[0, 1]*A[1, 0]


def test_transpose_index():
    assert X.T[i, j] == X[j, i]


def test_Identity_index():
    I = Identity(3)
    assert I[0, 0] == I[1, 1] == I[2, 2] == 1
    assert I[1, 0] == I[0, 1] == I[2, 1] == 0
    raises(IndexError, lambda: I[3, 3])


def test_block_index():
    I = Identity(3)
    Z = ZeroMatrix(3, 3)
    B = BlockMatrix([[I, I], [I, I]])
    e3 = ImmutableMatrix(eye(3))
    BB = BlockMatrix([[e3, e3], [e3, e3]])
    assert B[0, 0] == B[3, 0] == B[0, 3] == B[3, 3] == 1
    assert B[4, 3] == B[5, 1] == 0

    BB = BlockMatrix([[e3, e3], [e3, e3]])
    assert B.as_explicit() == BB.as_explicit()

    BI = BlockMatrix([[I, Z], [Z, I]])

    assert BI.as_explicit().equals(eye(6))


def test_slicing():
    A.as_explicit()[0, :]  # does not raise an error


def test_errors():
    raises(IndexError, lambda: Identity(2)[1, 2, 3, 4, 5])
    raises(IndexError, lambda: Identity(2)[[1, 2, 3, 4, 5]])


def test_matrix_expression_from_index_summation():
    from sympy.abc import a,b,c,d
    A = MatrixSymbol("A", k, k)
    B = MatrixSymbol("B", k, k)
    C = MatrixSymbol("C", k, k)

    expr = Sum(W[a,b]*X[b,c]*Z[c,d], (b, 0, l), (c, 0, m))
    assert MatrixExpr.from_index_summation(expr, a) == W*X*Z
    expr = Sum(W.T[b,a]*X[b,c]*Z[c,d], (b, 0, l), (c, 0, m))
    assert MatrixExpr.from_index_summation(expr, a) == W*X*Z
    expr = Sum(A[b, a]*B[b, c]*C[c, d], (b, 0, k), (c, 0, k))
    assert MatrixSymbol.from_index_summation(expr, a) == A.T*B*C
    expr = Sum(A[b, a]*B[c, b]*C[c, d], (b, 0, k), (c, 0, k))
    assert MatrixSymbol.from_index_summation(expr, a) == A.T*B.T*C
    expr = Sum(C[c, d]*A[b, a]*B[c, b], (b, 0, k), (c, 0, k))
    assert MatrixSymbol.from_index_summation(expr, a) == A.T*B.T*C
    expr = Sum(A[a, b] + B[a, b], (a, 0, k-1), (b, 0, k-1))
    assert MatrixExpr.from_index_summation(expr, a) == A + B
    expr = Sum((A[a, b] + B[a, b])*C[b, c], (b, 0, k-1))
    assert MatrixExpr.from_index_summation(expr, a) == (A+B)*C
    expr = Sum((A[a, b] + B[b, a])*C[b, c], (b, 0, k-1))
    assert MatrixExpr.from_index_summation(expr, a) == (A+B.T)*C
    expr = Sum(A[a, b]*A[b, c]*A[c, d], (b, 0, k-1), (c, 0, k-1))
    assert MatrixExpr.from_index_summation(expr, a) == MatMul(A, A, A)
    expr = Sum(A[a, b]*A[b, c]*B[c, d], (b, 0, k-1), (c, 0, k-1))
    assert MatrixExpr.from_index_summation(expr, a) == MatMul(A, A, B)

    # Parse nested sums:
    expr = Sum(A[a, b]*Sum(B[b, c]*C[c, d], (c, 0, k-1)), (b, 0, k-1))
    assert MatrixExpr.from_index_summation(expr, a) == A*B*C

    # Test Kronecker delta:
    expr = Sum(A[a, b]*KroneckerDelta(b, c)*B[c, d], (b, 0, k-1), (c, 0, k-1))
    assert MatrixExpr.from_index_summation(expr, a) == A*B
