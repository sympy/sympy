from sympy import (symbols, MatrixSymbol, MatPow, BlockMatrix, Symbol,
        Identity, ZeroMatrix, ImmutableMatrix, eye, Sum, Dummy)
from sympy.utilities.pytest import raises
from sympy.matrices.expressions.matexpr import MatrixElement

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
    n = symbols("n")
    Q2 = A**n
    assert Q2[0, 0] == MatrixElement(Q2, 0, 0)


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


def test_matrix_expression_to_indices():
    i, j = symbols("i, j")
    i1, i2, i3 = symbols("i_1:4")

    def replace_dummies(expr):
        repl = {i: Symbol(i.name) for i in expr.atoms(Dummy)}
        return expr.xreplace(repl)

    expr = W*X*Z
    assert replace_dummies(expr._entry(i, j)) == \
        Sum(W[i, i1]*X[i1, i2]*Z[i2, j], (i1, 0, l-1), (i2, 0, m-1))
    expr = Z.T*X.T*W.T
    assert replace_dummies(expr._entry(i, j)) == \
        Sum(W[j, i2]*X[i2, i1]*Z[i1, i], (i1, 0, m-1), (i2, 0, l-1))

    expr = W*X*Z + W*Y*Z
    assert replace_dummies(expr._entry(i, j)) == \
        Sum(W[i, i1]*X[i1, i2]*Z[i2, j], (i1, 0, l-1), (i2, 0, m-1)) +\
        Sum(W[i, i1]*Y[i1, i2]*Z[i2, j], (i1, 0, l-1), (i2, 0, m-1))

    expr = W*(X + Y)*Z
    assert replace_dummies(expr._entry(i, j)) == \
            Sum(W[i, i1]*(X[i1, i2] + Y[i1, i2])*Z[i2, j], (i1, 0, l-1), (i2, 0, m-1))

    expr = A*B**2*A
    #assert replace_dummies(expr._entry(i, j)) == \
    #        Sum(A[i, i1]*B[i1, i2]*B[i2, i3]*A[i3, j], (i1, 0, 1), (i2, 0, 1), (i3, 0, 1))
