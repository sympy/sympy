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
from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor, algebraic_tensor_product


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)
I4 = MatrixSymbol("I4", 4, 4)


def test_pure_tensor_single_factor():
    pt = AlgebraicPureTensor(A)
    assert pt is A


def test_pure_tensor_construction():
    pt = AlgebraicPureTensor(A, C)
    assert isinstance(pt, AlgebraicPureTensor)
    assert pt.is_AlgebraicPureTensor
    assert not pt.is_commutative
    assert pt.num_factors == 2
    assert pt.factors == (A, C)
    assert pt.tensor_shape == ((3, 4), (4, 5))


def test_pure_tensor_three_factors():
    pt = AlgebraicPureTensor(A, I4, C)
    assert pt.num_factors == 3
    assert pt.factors == (A, I4, C)
    assert pt.tensor_shape == ((3, 4), (4, 4), (4, 5))


def test_pure_tensor_shape_chain():
    pt = AlgebraicPureTensor(A, C)
    assert pt.tensor_shape == ((3, 4), (4, 5))
    pt2 = AlgebraicPureTensor(C, D)
    assert pt2.tensor_shape == ((4, 5), (3, 5))


def test_pure_tensor_no_args():
    with raises(ValueError):
        AlgebraicPureTensor()


def test_pure_tensor_no_shape():
    from sympy import Symbol
    x = Symbol("x")
    with raises(TypeError):
        AlgebraicPureTensor(x)


def test_pure_tensor_number_rejected():
    with raises(TypeError):
        AlgebraicPureTensor(5)


def test_pure_tensor_repr():
    pt = AlgebraicPureTensor(A, C)
    r = repr(pt)
    assert "AlgebraicPureTensor" in r
    assert "A" in r
    assert "C" in r


def test_pure_tensor_str():
    pt = AlgebraicPureTensor(A, C)
    assert "A" in str(pt)
    assert "C" in str(pt)


def test_tensor_product_convenience():
    pt = algebraic_tensor_product(A, C)
    assert isinstance(pt, AlgebraicPureTensor)
    assert pt.tensor_shape == ((3, 4), (4, 5))
