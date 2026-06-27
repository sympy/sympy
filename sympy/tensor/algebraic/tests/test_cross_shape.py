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
from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor, algebraic_tensor_product
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor, algebraic_zero_tensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
E = MatrixSymbol("E", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)


def test_pure_tensors_same_shape_addition():
    pt_a = AlgebraicPureTensor(A, I3)
    pt_b = AlgebraicPureTensor(B, I3)
    assert isinstance(pt_a, AlgebraicPureTensor)
    assert isinstance(pt_b, AlgebraicPureTensor)
    assert pt_a.tensor_shape == pt_b.tensor_shape == ((3, 4), (3, 3))

    result = pt_a + pt_b
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))
    assert result.is_Add
    assert result.is_AlgebraicTensor


def test_pure_tensors_different_shape_raises():
    pt_ac = AlgebraicPureTensor(A, C)
    pt_cd = AlgebraicPureTensor(C, D)
    assert pt_ac.tensor_shape == ((3, 4), (4, 5))
    assert pt_cd.tensor_shape == ((4, 5), (3, 5))
    assert pt_ac.tensor_shape != pt_cd.tensor_shape

    with raises(ShapeMismatchError):
        pt_ac + pt_cd


def test_chain_addition_same_shape():
    pt_d = AlgebraicPureTensor(D, I3)
    pt_e = AlgebraicPureTensor(E, I3)
    result = pt_d + pt_e
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 5), (3, 3))


def test_mul_then_add():
    pt_a = AlgebraicPureTensor(A, C)
    result = AlgebraicTensor(2 * pt_a, 3 * pt_a)
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (4, 5))


def test_zero_tensor_mixed_shape_error():
    z = AlgebraicZeroTensor(((3, 4), (3, 3)))
    pt_ac = AlgebraicPureTensor(A, C)
    with raises(ShapeMismatchError):
        z + pt_ac


def test_display_types():
    pt_a = AlgebraicPureTensor(A, I3)
    pt_b = AlgebraicPureTensor(B, I3)
    pt_c = AlgebraicPureTensor(C)
    sum_ab = pt_a + pt_b

    assert type(pt_a).__name__ == "AlgebraicPureTensor"
    assert type(pt_b).__name__ == "AlgebraicPureTensor"
    assert type(pt_c).__name__ == "MatrixSymbol"
    assert pt_c is C
    assert type(sum_ab).__name__ == "AlgebraicTensor"

    assert pt_a.tensor_shape == ((3, 4), (3, 3))
    assert pt_b.tensor_shape == ((3, 4), (3, 3))
    assert sum_ab.tensor_shape == ((3, 4), (3, 3))


def test_add_dispatcher_pure_tensor():
    from sympy.core.add import Add
    pt_a = AlgebraicPureTensor(A, I3)
    pt_b = AlgebraicPureTensor(B, I3)
    result = pt_a + pt_b
    assert isinstance(result, AlgebraicTensor)
    assert not isinstance(result, Add)


def test_add_dispatcher_algebraic_tensor():
    at1 = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    at2 = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    result = at1 + at2
    assert isinstance(result, AlgebraicTensor)


def test_init_reexports():
    from sympy.tensor.algebraic import (
        AlgebraicTensor as AT,
        AlgebraicPureTensor as PT,
        ShapeMismatchError as SME,
        AlgebraicZeroTensor as ZT,
        algebraic_tensor_product as tp,
        algebraic_zero_tensor as zt,
    )
    assert AT is AlgebraicTensor
    assert PT is AlgebraicPureTensor
    assert SME is ShapeMismatchError
    assert ZT is AlgebraicZeroTensor
    assert tp is algebraic_tensor_product
    assert zt is algebraic_zero_tensor
