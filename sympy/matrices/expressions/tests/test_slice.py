from __future__ import annotations
from sympy.matrices.expressions.slice import MatrixSlice
from sympy.matrices.expressions import MatrixSymbol
from sympy.abc import a, b, c, d, k, l, m, n
from sympy.testing.pytest import raises, XFAIL
from sympy.functions.elementary.integers import ceiling
from sympy.assumptions import assuming, Q


X = MatrixSymbol('X', n, m)
Y = MatrixSymbol('Y', m, k)
X10 = MatrixSymbol('X', 10, 20)
X_sq = MatrixSymbol('X_sq', 10, 10)


def test_shape():
    B = MatrixSlice(X10, (1, 5), (2, 4))
    assert B.shape == (4, 2)

    B_sym = MatrixSlice(X, (a, b), (c, d))
    from sympy.core.containers import Tuple
    assert Tuple(*B_sym.shape).subs(
        {n: 10, m: 20, a: 1, b: 5, c: 2, d: 4}
    ) == (4, 2)


def test_entry():
    B = MatrixSlice(X10, (1, 5), (2, 4))
    assert B[0, 0] == X10[1, 2]
    assert B[k, l] == X10[1+k, 2+l]
    raises(IndexError, lambda: MatrixSlice(X10, 1, (2, 5))[1, 0])

    assert X10[1::2, :][1, 3] == X10[1+2, 3]
    assert X10[:, 1::2][3, 1] == X10[3, 1+2]

    B_sym = MatrixSlice(X, (a, b), (c, d))
    assert B_sym[0, 0].subs({n: 10, m: 20, a: 1, c: 2}) == X10[1, 2]


def test_on_diag():
    assert not MatrixSlice(X_sq, (1, 5), (2, 6)).on_diag
    assert MatrixSlice(X_sq, (1, 5), (1, 5)).on_diag


def test_inputs():
    assert MatrixSlice(X10, 1, (2, 5)) == MatrixSlice(X10, (1, 2), (2, 5))
    assert MatrixSlice(X10, 1, (2, 5)).shape == (1, 3)


def test_slicing():
    assert X10[1:5, 2:4] == MatrixSlice(X10, (1, 5), (2, 4))
    assert X10[1, 2:4] == MatrixSlice(X10, 1, (2, 4))
    assert X10[1:5, :].shape == (4, 20)
    assert X10[:, 1:5].shape == (10, 4)
    assert X10[2, :] == MatrixSlice(X10, 2, (0, 20))

    from sympy.core.symbol import Symbol
    n_pos = Symbol('n', integer=True, positive=True)
    m_pos = Symbol('m', integer=True, positive=True)
    X_pos = MatrixSymbol('X_pos', n_pos, m_pos)
    assert X_pos[::2, ::2].shape == (ceiling(n_pos/2), ceiling(m_pos/2))

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
    # Integer index out of bounds raises IndexError
    raises(IndexError, lambda: X[0:9, 22])
    raises(IndexError, lambda: X[12, 2])
    raises(IndexError, lambda: X[-11, 2])
    # Zero step raises ValueError
    raises(
        ValueError,
        lambda: MatrixSlice(X, slice(None, None, 0), slice(None))
    )


def test_slice_clipping():
    X = MatrixSymbol('x', 10, 20)
    # Slices beyond bounds clip to matrix shape (Python/NumPy array slicing)
    assert X[0:12, :].shape == (10, 20)
    assert X[0:100, 0:100].shape == (10, 20)
    assert X[-1:5, :].shape == (0, 20)
    assert X[-10:-1, :].shape == (9, 20)


def test_negative_step_slicing():
    from sympy.matrices import Matrix
    M = Matrix(4, 4, range(16))
    X = MatrixSymbol('X', 4, 4)

    # Reverse rows and cols
    rev = X[::-1, ::-1]
    assert rev.shape == (4, 4)
    M_rev = rev.subs(X, M).as_explicit()
    assert M_rev == M[::-1, ::-1]
    assert M_rev[0, 0] == M[3, 3]

    # Negative step with non-unit stride
    s = X[3:0:-2, :]
    assert s.shape == (2, 4)
    assert s.subs(X, M).as_explicit() == M[3:0:-2, :]


def test_symbolic_slicing_substitution_invariance():
    # Issue #30999: shape and entries must be consistent under substitutions
    from sympy.matrices import Matrix
    from sympy.core.symbol import Symbol
    A = MatrixSymbol('A', 5, 4)
    M = Matrix(5, 4, range(20))
    # Slicing from negative index with step 2
    sl = MatrixSlice(A, slice(-4, None, 2), slice(None))
    assert sl.shape == (2, 4)
    assert sl.subs(A, M).as_explicit() == M[-4::2, :]

    # Slicing beyond bounds
    sl_out = MatrixSlice(A, slice(0, 100), slice(None))
    assert sl_out.shape == (5, 4)
    assert sl_out.subs(A, M).as_explicit() == M

    # Negative step on symbolic matrix
    sl_rev = MatrixSlice(A, slice(None, None, -1), slice(None))
    assert sl_rev.shape == (5, 4)
    assert sl_rev.subs(A, M).as_explicit() == M[::-1, :]

    # S.Y. Lee Test 1: Slicing negative index on symbolic dimension matrix
    k = Symbol('k', integer=True, positive=True)
    X_k = MatrixSymbol('X', k, 4)
    s_neg = MatrixSlice(X_k, slice(-6, None), slice(None)).subs(k, 5)
    assert s_neg.shape == (5, 4)

    # S.Y. Lee Test 2: Slicing with symbolic start parameter
    sym_a = Symbol('a')
    s_sym = MatrixSlice(A, slice(sym_a, None, 2), slice(None))
    assert s_sym.subs(sym_a, -4).shape == (2, 4)


def test_issue_18411():
    # Issue #18411: Matrix slicing inconsistency with python for negative step
    from sympy.matrices import Matrix
    A = MatrixSymbol('A', 4, 4)
    M = Matrix(4, 4, range(16))

    # S.Y. Lee test case from #18411:
    sl1 = A[5:0:-1, 5:0:-1]
    assert sl1.shape == (3, 3)
    assert sl1.subs(A, M).as_explicit() == M.as_explicit()[5:0:-1, 5:0:-1]

    # Oscar Benjamin test case from #18411:
    sl2 = A[6::-1, 6::-1]
    assert sl2.shape == (4, 4)
    assert sl2.subs(A, M).as_explicit() == M[::-1, ::-1]


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
