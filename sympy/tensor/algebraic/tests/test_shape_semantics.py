from __future__ import annotations

from contextlib import contextmanager

@contextmanager
def raises(expected_exc):
    """Minimal pytest.raises replacement for environments without pytest."""
    try:
        yield
    except expected_exc:
        pass
    except Exception as e:
        raise AssertionError(f"Expected {expected_exc.__name__} but got {type(e).__name__}: {e}")
    else:
        raise AssertionError(f"Expected {expected_exc.__name__} but no exception was raised")

from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.algebraic_tensor import (
    AlgebraicTensor,
    ShapeMismatchError,
)
from sympy.tensor.algebraic.pure_tensor import PureTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)


def test_tensor_shape_vector_matrix():
    v = MatrixSymbol("v", 3, 1)
    M = MatrixSymbol("M", 1, 5)
    pt = PureTensor(v, M)
    assert pt.tensor_shape == ((3, 1), (1, 5))


def test_tensor_shape_transposed_vector():
    x = MatrixSymbol("x", 1, 3)
    y = MatrixSymbol("y", 3, 1)
    pt = PureTensor(x, y)
    assert pt.tensor_shape == ((1, 3), (3, 1))


def test_tensor_shape_mixed_spaces():
    M1 = MatrixSymbol("M1", 2, 3)
    M2 = MatrixSymbol("M2", 5, 7)
    pt = PureTensor(M1, M2)
    assert pt.tensor_shape == ((2, 3), (5, 7))


def test_tensor_shape_four_factors():
    M1 = MatrixSymbol("M1", 2, 3)
    M2 = MatrixSymbol("M2", 3, 4)
    M3 = MatrixSymbol("M3", 4, 5)
    M4 = MatrixSymbol("M4", 5, 6)
    pt = PureTensor(M1, M2, M3, M4)
    assert pt.tensor_shape == ((2, 3), (3, 4), (4, 5), (5, 6))
    assert pt.num_factors == 4


def test_algebraic_tensor_add_mixed_factor_spaces():
    M1 = MatrixSymbol("M1", 2, 3)
    M2 = MatrixSymbol("M2", 5, 7)
    N1 = MatrixSymbol("N1", 2, 3)
    N2 = MatrixSymbol("N2", 5, 7)
    pt1 = PureTensor(M1, M2)
    pt2 = PureTensor(N1, N2)
    result = pt1 + pt2
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((2, 3), (5, 7))


def test_algebraic_tensor_different_factor_count_raises():
    pt_two = PureTensor(A, C)
    pt_one = PureTensor(D)
    with raises(ShapeMismatchError):
        pt_two + pt_one
