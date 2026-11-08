from __future__ import annotations
from sympy.matrices.expressions.slice import MatrixSlice
from sympy.matrices.expressions import MatrixSymbol
from sympy.abc import a, b, c, d, k, l, m, n
from sympy.testing.pytest import raises, XFAIL
from sympy.functions.elementary.integers import ceiling
from sympy.assumptions import assuming, Q


X = MatrixSymbol('X', n, m)
Y = MatrixSymbol('Y', m, k)

def test_shape():
    B = MatrixSlice(X, (a, b), (c, d))
    assert B.shape == (b - a, d - c)

def test_entry():
    B = MatrixSlice(X, (a, b), (c, d))
    assert B[0,0] == X[a, c]
    assert B[k,l] == X[a+k, c+l]
    raises(IndexError, lambda : MatrixSlice(X, 1, (2, 5))[1, 0])

    assert X[1::2, :][1, 3] == X[1+2, 3]
    assert X[:, 1::2][3, 1] == X[3, 1+2]

def test_on_diag():
    assert not MatrixSlice(X, (a, b), (c, d)).on_diag
    assert MatrixSlice(X, (a, b), (a, b)).on_diag

def test_inputs():
    assert MatrixSlice(X, 1, (2, 5)) == MatrixSlice(X, (1, 2), (2, 5))
    assert MatrixSlice(X, 1, (2, 5)).shape == (1, 3)

def test_slicing():
    assert X[1:5, 2:4] == MatrixSlice(X, (1, 5), (2, 4))
    assert X[1, 2:4] == MatrixSlice(X, 1, (2, 4))
    assert X[1:5, :].shape == (4, X.shape[1])
    assert X[:, 1:5].shape == (X.shape[0], 4)

    assert X[::2, ::2].shape == (ceiling(n/2), ceiling(m/2))
    assert X[2, :] == MatrixSlice(X, 2, (0, m))
    assert X[k, :] == MatrixSlice(X, k, (0, m))

    # Concrete slices with step > 1 must not truncate boundary elements
    from sympy.matrices import Matrix
    X5 = MatrixSymbol('X5', 5, 5)
    assert X5[::2, ::2].shape == (3, 3)
    M5 = Matrix(5, 5, lambda i, j: 10*i + j)
    assert X5[::2, ::2].subs(X5, M5).as_explicit() == M5[::2, ::2]

    X6 = MatrixSymbol('X6', 6, 6)
    assert X6[::2, ::2].shape == (3, 3)
    M6 = Matrix(6, 6, lambda i, j: 10*i + j)
    assert X6[::2, ::2].subs(X6, M6).as_explicit() == M6[::2, ::2]

    X4 = MatrixSymbol('X4', 4, 4)
    assert X4[::2, ::2].shape == (2, 2)
    M4 = Matrix(4, 4, lambda i, j: 10*i + j)
    assert X4[::2, ::2].subs(X4, M4).as_explicit() == M4[::2, ::2]

    X7 = MatrixSymbol('X7', 7, 7)
    assert X7[::3, ::3].shape == (3, 3)
    M7 = Matrix(7, 7, lambda i, j: 10*i + j)
    assert X7[::3, ::3].subs(X7, M7).as_explicit() == M7[::3, ::3]
    assert X7[1:6:2, :].shape == (3, 7)
    assert X7[2:7:3, :].shape == (2, 7)


def test_exceptions():
    X = MatrixSymbol('x', 10, 20)
    raises(IndexError, lambda: X[0:12, 2])
    raises(IndexError, lambda: X[0:9, 22])
    raises(IndexError, lambda: X[-1:5, 2])

@XFAIL
def test_symmetry():
    X = MatrixSymbol('x', 10, 10)
    Y = X[:5, 5:]
    with assuming(Q.symmetric(X)):
        assert Y.T == X[5:, :5]

def test_slice_of_slice():
    X = MatrixSymbol('x', 10, 10)
    assert X[2, :][:, 3][0, 0] == X[2, 3]
    assert X[:5, :5][:4, :4] == X[:4, :4]
    assert X[1:5, 2:6][1:3, 2] == X[2:4, 4]
    assert X[1:9:2, 2:6][1:3, 2] == X[3:7:2, 4]

def test_negative_index():
    X = MatrixSymbol('x', 10, 10)
    assert X[-1, :] == X[9, :]
