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

def test_exceptions():
    X = MatrixSymbol('x', 10, 20)
    # Out-of-range single indices raise, like indexing a Python sequence.
    raises(IndexError, lambda: X[12, 2])
    raises(IndexError, lambda: X[0:9, 22])
    # Out-of-range slice bounds are clamped to match Python slicing (#18411).
    assert X[0:12, 2] == X[0:10, 2]
    assert X[-1:5, 2].shape == (0, 1)


def test_slice_matches_explicit():
    # A symbolic slice must agree with slicing the explicit matrix, including
    # stepped, negative-step and out-of-range slices. Non-regression test for
    # https://github.com/sympy/sympy/issues/18411
    A = MatrixSymbol('A', 4, 5)
    explicit = A.as_explicit()
    slices = [
        slice(1, 4, 2), slice(5, 0, -1), slice(6, None, -1), slice(0, 5, -1),
        slice(None, None, -1), slice(0, 10), slice(0, 4, 3), slice(3, None, -2),
    ]
    for rs in slices:
        for cs in slices:
            expected = explicit[rs, cs]
            if 0 in expected.shape:
                # empty result: just check the shape agrees
                assert A[rs, cs].shape == expected.shape
            else:
                assert A[rs, cs].as_explicit() == expected

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
