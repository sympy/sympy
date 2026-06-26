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
from sympy.tensor.algebraic.zero_tensor import ZeroTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)


def test_algebraic_tensor_add_same_shape():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    result = pt_a + pt_b
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))
    assert result.is_Add


def test_algebraic_tensor_add_different_shape():
    pt_ac = PureTensor(A, C)
    pt_cd = PureTensor(C, D)
    with raises(ShapeMismatchError):
        pt_ac + pt_cd


def test_algebraic_tensor_three_terms_same_shape():
    ab = PureTensor(A, I3)
    bc = PureTensor(B, I3)
    result = ab + bc + PureTensor(A, I3)
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))


def test_algebraic_tensor_flatten_nested():
    t1 = PureTensor(A, I3) + PureTensor(B, I3)
    t2 = PureTensor(A, I3) + PureTensor(B, I3)
    result = AlgebraicTensor(t1, t2)
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_single_arg_unwrap():
    pt = PureTensor(A, C)
    result = AlgebraicTensor(pt)
    assert result is pt


def test_algebraic_tensor_single_algebraic_unwrap():
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = AlgebraicTensor(at)
    assert result is at


def test_algebraic_tensor_single_zerotensor_unwrap():
    z = ZeroTensor(((3, 4),))
    result = AlgebraicTensor(z)
    assert result is z


def test_algebraic_tensor_no_args():
    with raises(ValueError):
        AlgebraicTensor()


def test_algebraic_tensor_with_coeff():
    pt = PureTensor(A, C)
    result = AlgebraicTensor(2 * pt, 3 * pt)
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (4, 5))


def test_algebraic_tensor_with_zero_tensor_anchor():
    z = ZeroTensor(((3, 4), (4, 5)))
    pt = PureTensor(A, C)
    result = AlgebraicTensor(pt, z)
    assert isinstance(result, AlgebraicTensor)
    assert result.has_zero_term()


def test_algebraic_tensor_neg():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    at = pt_a + pt_b
    neg = -at
    assert isinstance(neg, AlgebraicTensor)


def test_algebraic_tensor_terms():
    at = PureTensor(A, I3) + PureTensor(B, I3)
    terms = at.terms
    assert len(terms) >= 2


def test_algebraic_tensor_repr_str():
    at = PureTensor(A, I3) + PureTensor(B, I3)
    assert "AlgebraicTensor" in repr(at)
    assert "+" in str(at)


def test_algebraic_tensor_has_zero_term_false():
    at = PureTensor(A, I3) + PureTensor(B, I3)
    assert not at.has_zero_term()


def test_algebraic_tensor_has_zero_term_true():
    z = ZeroTensor(((3, 4), (3, 3)))
    at = AlgebraicTensor(PureTensor(A, I3), z)
    assert at.has_zero_term()


def test_algebraic_tensor_is_flags():
    at = PureTensor(A, I3) + PureTensor(B, I3)
    assert at.is_AlgebraicTensor
    assert at.is_Add
