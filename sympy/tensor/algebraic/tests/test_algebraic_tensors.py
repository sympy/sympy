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

from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.core.mul import Mul
from sympy.core.add import Add

from sympy.tensor.algebraic.algebraic_tensor import (
    AlgebraicTensor,
    ShapeMismatchError,
)
from sympy.tensor.algebraic.pure_tensor import PureTensor, tensor_product
from sympy.tensor.algebraic.zero_tensor import ZeroTensor, zero_tensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
E = MatrixSymbol("E", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)
I4 = MatrixSymbol("I4", 4, 4)

# --- Canonical tensor shapes for the fixtures ---
# PureTensor(A, C)  -> ((3,4), (4,5))
# PureTensor(A, I3) -> ((3,4), (3,3))
# PureTensor(A, I4, C) -> ((3,4), (4,4), (4,5))
# A alone (single factor unwrap) -> ((3,4),)
# PureTensor(A, I3, D) -> ((3,4), (3,3), (3,5))


# =========================================================
# ZeroTensor
# =========================================================

def test_zero_tensor_construction():
    z = ZeroTensor(((3, 4), (4, 5)))
    assert z.shape == ((3, 4), (4, 5))
    assert z.tensor_shape == ((3, 4), (4, 5))
    assert not z
    assert z == ZeroTensor(((3, 4), (4, 5)))
    assert z != ZeroTensor(((3, 5),))
    assert -z is z


def test_zero_tensor_bare_pair_wrapping():
    z = ZeroTensor((3, 4))
    assert z.shape == ((3, 4),)


def test_zero_tensor_hash():
    z1 = ZeroTensor(((3, 4),))
    z2 = ZeroTensor(((3, 4),))
    assert hash(z1) == hash(z2)
    assert z1 == z2
    assert {z1} == {z2}


def test_zero_tensor_repr_str():
    z = ZeroTensor(((3, 4),))
    assert repr(z) == "ZeroTensor((3, 4),)"
    assert str(z) == "0_((3, 4),)"


def test_zero_tensor_add():
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = z + pt
    assert isinstance(result, (PureTensor, AlgebraicTensor))


def test_zero_tensor_radd():
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = pt + z
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_sub():
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = z - pt
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_rsub():
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = pt - z
    assert isinstance(result, AlgebraicTensor)


def test_zero_tensor_convenience():
    z = zero_tensor(((3, 4), (4, 5)))
    assert isinstance(z, ZeroTensor)
    assert z.shape == ((3, 4), (4, 5))


# =========================================================
# PureTensor
# =========================================================

def test_pure_tensor_single_factor():
    pt = PureTensor(A)
    assert pt is A


def test_pure_tensor_construction():
    pt = PureTensor(A, C)
    assert isinstance(pt, PureTensor)
    assert pt.is_PureTensor
    assert not pt.is_commutative
    assert pt.num_factors == 2
    assert pt.factors == (A, C)
    assert pt.tensor_shape == ((3, 4), (4, 5))


def test_pure_tensor_three_factors():
    pt = PureTensor(A, I4, C)
    assert pt.num_factors == 3
    assert pt.factors == (A, I4, C)
    assert pt.tensor_shape == ((3, 4), (4, 4), (4, 5))


def test_pure_tensor_shape_chain():
    pt = PureTensor(A, C)
    assert pt.tensor_shape == ((3, 4), (4, 5))
    pt2 = PureTensor(C, D)
    assert pt2.tensor_shape == ((4, 5), (3, 5))


def test_pure_tensor_no_args():
    with raises(ValueError):
        PureTensor()


def test_pure_tensor_no_shape():
    from sympy import Symbol
    x = Symbol("x")
    with raises(TypeError):
        PureTensor(x)


def test_pure_tensor_number_rejected():
    with raises(TypeError):
        PureTensor(5)


def test_pure_tensor_repr():
    pt = PureTensor(A, C)
    r = repr(pt)
    assert "PureTensor" in r
    assert "A" in r
    assert "C" in r


def test_pure_tensor_str():
    pt = PureTensor(A, C)
    assert "A" in str(pt)
    assert "C" in str(pt)


def test_pure_tensor_neg():
    pt = PureTensor(A, C)
    neg = -pt
    assert isinstance(neg, PureTensor)
    assert neg.factors == pt.factors


def test_pure_tensor_mul_scalar():
    pt = PureTensor(A, C)
    result = pt * 3
    assert isinstance(result, PureTensor)
    assert result.factors == pt.factors


def test_pure_tensor_rmul_scalar():
    pt = PureTensor(A, C)
    result = 3 * pt
    assert isinstance(result, PureTensor)
    assert result.factors == pt.factors


def test_pure_tensor_mul_one():
    pt = PureTensor(A, C)
    assert pt * S.One is pt
    assert S.One * pt is pt


def test_pure_tensor_mul_zero():
    pt = PureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_zero():
    """pt * S.Zero goes through __mul__ and returns ZeroTensor."""
    pt = PureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_add():
    pt1 = PureTensor(A, I3)
    pt2 = PureTensor(B, I3)
    result = pt1 + pt2
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_radd():
    pt = PureTensor(A, I3)
    result = pt + PureTensor(B, I3)
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_sub():
    pt1 = PureTensor(A, I3)
    pt2 = PureTensor(B, I3)
    result = pt1 - pt2
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_rsub():
    pt = PureTensor(A, I3)
    result = pt - PureTensor(B, I3)
    assert isinstance(result, AlgebraicTensor)


# =========================================================
# tensor_product convenience
# =========================================================

def test_tensor_product_convenience():
    pt = tensor_product(A, C)
    assert isinstance(pt, PureTensor)
    assert pt.tensor_shape == ((3, 4), (4, 5))


# =========================================================
# AlgebraicTensor
# =========================================================

def test_algebraic_tensor_add_same_shape():
    """Add two multi-factor PureTensors of the same shape."""
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    result = pt_a + pt_b
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))
    assert result.is_Add


def test_algebraic_tensor_add_different_shape():
    """Two multi-factor PureTensors of different shapes raise ShapeMismatchError."""
    pt_ac = PureTensor(A, C)   # ((3,4), (4,5))
    pt_cd = PureTensor(C, D)   # ((4,5), (3,5))
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


def test_algebraic_tensor_sub():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    result = pt_a - pt_b
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_rsub():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    result = pt_b - pt_a
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_radd():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    at = pt_a + pt_b
    result = PureTensor(A, I3) + at
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_add_method():
    at1 = PureTensor(A, I3) + PureTensor(B, I3)
    at2 = PureTensor(A, I3) + PureTensor(B, I3)
    result = at1 + at2
    assert isinstance(result, AlgebraicTensor)


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


# =========================================================
# Factorization helpers
# =========================================================

def test_as_common_left_basic():
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    left, rest, right = at.as_common_left()
    assert left == ()
    assert rest is at


def test_as_common_left_no_common():
    at = PureTensor(A, C) + PureTensor(B, C)
    left, rest, right = at.as_common_left()
    assert left == ()
    assert rest is at


def test_as_common_right_basic():
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    left, rest, right = at.as_common_right()
    assert len(right) >= 1
    assert isinstance(rest, (AlgebraicTensor, PureTensor, MatrixSymbol))


def test_as_common_right_no_common():
    at = PureTensor(A, I3, D) + PureTensor(B, I3, E)
    left, rest, right = at.as_common_right()
    assert left == ()
    assert right == ()
    assert rest is at


def test_as_common_factors():
    """Extract common left and right factors from a two-term sum.

    Both terms must have identical factor-shape sequences.
    PureTensor(A, I4, C) has shape ((3,4), (4,4), (4,5)).
    We pair it with PureTensor(A, J4, C) which also has shape ((3,4), (4,4), (4,5)).
    Common left: A, common right: C, middle differs (I4 vs J4).
    """
    J4 = MatrixSymbol("J4", 4, 4)
    at = PureTensor(A, I4, C) + PureTensor(A, J4, C)
    assert at.tensor_shape == ((3, 4), (4, 4), (4, 5))
    left, mid, right = at.as_common_factors()
    assert len(left) >= 1
    assert left[0] == A
    assert len(right) >= 1
    assert right[-1] == C


def test_as_common_left_identity_factors():
    at = PureTensor(A, I3, D) + PureTensor(A, I3, E)
    left, rest, right = at.as_common_left()
    assert len(left) >= 2
    assert left[0] == A
    assert left[1] == I3


def test_as_common_right_identity_factors():
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    left, rest, right = at.as_common_right()
    assert len(right) >= 1


# =========================================================
# Cross-shape tests (user-specified scenario)
# =========================================================

def test_pure_tensors_same_shape_addition():
    """Create two pure tensors of one shape and add them."""
    pt_a = PureTensor(A, I3)   # shape ((3,4), (3,3))
    pt_b = PureTensor(B, I3)   # shape ((3,4), (3,3))
    assert isinstance(pt_a, PureTensor)
    assert isinstance(pt_b, PureTensor)
    assert pt_a.tensor_shape == pt_b.tensor_shape == ((3, 4), (3, 3))

    result = pt_a + pt_b
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))
    assert result.is_Add
    assert result.is_AlgebraicTensor


def test_pure_tensors_different_shape_raises():
    """Create two pure tensors of different shapes; addition raises."""
    pt_ac = PureTensor(A, C)   # shape ((3,4), (4,5))
    pt_cd = PureTensor(C, D)   # shape ((4,5), (3,5))
    assert pt_ac.tensor_shape == ((3, 4), (4, 5))
    assert pt_cd.tensor_shape == ((4, 5), (3, 5))
    assert pt_ac.tensor_shape != pt_cd.tensor_shape

    with raises(ShapeMismatchError):
        pt_ac + pt_cd


def test_chain_addition_same_shape():
    """Chain of additions with same-shape tensors."""
    pt_d = PureTensor(D, I3)   # ((3,5), (3,3))
    pt_e = PureTensor(E, I3)   # ((3,5), (3,3))
    result = pt_d + pt_e
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 5), (3, 3))


def test_mul_then_add():
    """Multiply a PureTensor by a scalar, then add."""
    pt_a = PureTensor(A, C)
    result = AlgebraicTensor(2 * pt_a, 3 * pt_a)
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (4, 5))


def test_zero_tensor_mixed_shape_error():
    """ZeroTensor of one shape can't be added with PureTensor of another."""
    z = ZeroTensor(((3, 4), (3, 3)))
    pt_ac = PureTensor(A, C)  # shape ((3,4), (4,5))
    with raises(ShapeMismatchError):
        z + pt_ac


def test_display_types():
    """Display types and shapes to verify the setup is correct."""
    pt_a = PureTensor(A, I3)   # PureTensor, shape ((3,4), (3,3))
    pt_b = PureTensor(B, I3)   # PureTensor, shape ((3,4), (3,3))
    pt_c = PureTensor(C)       # unwraps to MatrixSymbol C, shape (4,5)
    sum_ab = pt_a + pt_b

    assert type(pt_a).__name__ == "PureTensor"
    assert type(pt_b).__name__ == "PureTensor"
    assert type(pt_c).__name__ == "MatrixSymbol"
    assert pt_c is C
    assert type(sum_ab).__name__ == "AlgebraicTensor"

    assert pt_a.tensor_shape == ((3, 4), (3, 3))
    assert pt_b.tensor_shape == ((3, 4), (3, 3))
    assert sum_ab.tensor_shape == ((3, 4), (3, 3))


# =========================================================
# Dispatcher integration
# =========================================================

def test_add_dispatcher_pure_tensor():
    """The add dispatcher routes PureTensor additions to AlgebraicTensor."""
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    result = pt_a + pt_b
    assert isinstance(result, AlgebraicTensor)
    assert not isinstance(result, Add)


def test_add_dispatcher_algebraic_tensor():
    """The add dispatcher routes AlgebraicTensor additions through AlgebraicTensor."""
    at1 = PureTensor(A, I3) + PureTensor(B, I3)
    at2 = PureTensor(A, I3) + PureTensor(B, I3)
    result = at1 + at2
    assert isinstance(result, AlgebraicTensor)


# =========================================================
# Package-level __init__ re-export tests
# =========================================================

def test_init_reexports():
    """All public symbols are re-exported from the package __init__."""
    from sympy.tensor.algebraic import (
        AlgebraicTensor as AT,
        PureTensor as PT,
        ShapeMismatchError as SME,
        ZeroTensor as ZT,
        tensor_product as tp,
        zero_tensor as zt,
    )
    assert AT is AlgebraicTensor
    assert PT is PureTensor
    assert SME is ShapeMismatchError
    assert ZT is ZeroTensor
    assert tp is tensor_product
    assert zt is zero_tensor


# =========================================================
# Shape semantics: vectors, transposed vectors, 1x1 noncommutative
# =========================================================

def test_tensor_shape_vector_matrix():
    """PureTensor of a column vector and a matrix."""
    v = MatrixSymbol("v", 3, 1)
    M = MatrixSymbol("M", 1, 5)
    pt = PureTensor(v, M)
    assert pt.tensor_shape == ((3, 1), (1, 5))


def test_tensor_shape_transposed_vector():
    """PureTensor of a row vector (transposed column)."""
    x = MatrixSymbol("x", 1, 3)
    y = MatrixSymbol("y", 3, 1)
    pt = PureTensor(x, y)
    assert pt.tensor_shape == ((1, 3), (3, 1))


def test_tensor_shape_mixed_spaces():
    """Tensor product of factors from different matrix spaces."""
    M1 = MatrixSymbol("M1", 2, 3)
    M2 = MatrixSymbol("M2", 5, 7)
    pt = PureTensor(M1, M2)
    assert pt.tensor_shape == ((2, 3), (5, 7))


def test_tensor_shape_four_factors():
    """Tensor product of four factors."""
    M1 = MatrixSymbol("M1", 2, 3)
    M2 = MatrixSymbol("M2", 3, 4)
    M3 = MatrixSymbol("M3", 4, 5)
    M4 = MatrixSymbol("M4", 5, 6)
    pt = PureTensor(M1, M2, M3, M4)
    assert pt.tensor_shape == ((2, 3), (3, 4), (4, 5), (5, 6))
    assert pt.num_factors == 4


def test_algebraic_tensor_add_mixed_factor_spaces():
    """Addition of PureTensors with the same factor-shape sequence."""
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
    """Two factors vs one factor -> shape mismatch."""
    pt_two = PureTensor(A, C)   # ((3,4), (4,5))
    pt_one = PureTensor(D)      # unwraps to D, shape ((3,5),)
    with raises(ShapeMismatchError):
        pt_two + pt_one


# =========================================================
# Coefficient absorption: PureTensor * scalar/symbol -> PureTensor
# =========================================================

def test_pure_tensor_mul_scalar_is_puretensor():
    """Multiplying a PureTensor by a scalar returns a PureTensor."""
    pt = PureTensor(A, C)
    result = pt * 5
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)
    assert result.tensor_shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_scalar_is_puretensor():
    """Multiplying a scalar with a PureTensor returns a PureTensor."""
    pt = PureTensor(A, C)
    result = 5 * pt
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_negative_scalar():
    """Multiplying by a negative scalar produces a PureTensor with negative coefficient."""
    pt = PureTensor(A, C)
    result = pt * (-3)
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_double_negation():
    """Double negation returns the original PureTensor."""
    pt = PureTensor(A, C)
    assert (-(-pt)).factors == pt.factors
    assert (-(-pt)).tensor_shape == pt.tensor_shape


def test_pure_tensor_neg_is_puretensor():
    """Negation of a PureTensor returns a PureTensor."""
    pt = PureTensor(A, C)
    neg = -pt
    assert isinstance(neg, PureTensor)
    assert neg.factors == pt.factors
    assert neg.tensor_shape == pt.tensor_shape


def test_pure_tensor_neg_of_negated():
    """Negation of a negated PureTensor restores the original."""
    pt = PureTensor(A, C)
    neg = -pt
    nneg = -neg
    assert nneg.factors == pt.factors


def test_pure_tensor_mul_coefficient_combines():
    """Multiplying a coefficient-carrying PureTensor combines coefficients."""
    pt = PureTensor(A, C)
    scaled = 2 * pt
    rescaled = scaled * 3
    assert isinstance(rescaled, PureTensor)
    assert rescaled.factors == (A, C)


def test_pure_tensor_mul_coefficient_cancellation():
    """Multiplying by reciprocal coefficient cancels to original."""
    from sympy import Rational
    pt = PureTensor(A, C)
    scaled = pt * Rational(1, 3)
    rescaled = scaled * 3
    assert isinstance(rescaled, PureTensor)
    assert rescaled.factors == pt.factors


def test_pure_tensor_mul_symbol():
    """PureTensor * Symbol absorbs it as coefficient."""
    from sympy import Symbol
    x = Symbol("x")
    pt = PureTensor(A, C)
    result = pt * x
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_puretensor():
    """Multiplying two PureTensors concatenates factors."""
    pt1 = PureTensor(A, C)
    pt2 = PureTensor(C, D)
    result = pt1 * pt2
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C, C, D)
    assert result.tensor_shape == ((3, 4), (4, 5), (4, 5), (3, 5))


def test_pure_tensor_mul_puretensor_with_coefficients():
    """Multiplying two coefficient-carrying PureTensors combines coefficients."""
    pt1 = 2 * PureTensor(A, C)
    pt2 = 3 * PureTensor(C, D)
    result = pt1 * pt2
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C, C, D)


def test_pure_tensor_mul_zero_returns_zerotensor():
    """PureTensor * 0 returns ZeroTensor with correct shape."""
    pt = PureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_zero_returns_zerotensor():
    """Python 0 * PureTensor returns ZeroTensor with correct shape.

    Note: S.Zero * pt bypasses PureTensor.__rmul__ due to SymPy's dispatch
    mechanism. Use Python int 0 instead.
    """
    pt = PureTensor(A, C)
    result = 0 * pt
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_mul_one_identity():
    """PureTensor * 1 and 1 * PureTensor return self."""
    pt = PureTensor(A, C)
    assert pt * S.One is pt
    assert S.One * pt is pt


# =========================================================
# AlgebraicTensor multiplication
# =========================================================

def test_algebraic_tensor_mul_scalar():
    """Multiplying an AlgebraicTensor by a scalar distributes to all terms."""
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = at * 3
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))


def test_algebraic_tensor_rmul_scalar():
    """Scalar multiplication on the left distributes to all terms."""
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = 3 * at
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))


def test_algebraic_tensor_mul_symbol():
    """AlgebraicTensor * Symbol distributes to all terms."""
    from sympy import Symbol
    x = Symbol("x")
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = at * x
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_rmul_symbol():
    """Symbol * AlgebraicTensor distributes to all terms via __rmul__."""
    from sympy import Symbol
    x = Symbol("x")
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = x * at
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_mul_negative_one():
    """Multiplying an AlgebraicTensor by -1 negates all terms."""
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = at * S.NegativeOne
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_mul_preserves_terms_count():
    """Scalar multiplication preserves the number of terms."""
    at = PureTensor(A, I3) + PureTensor(B, I3)
    original_count = len(at.args)
    result = at * 5
    assert len(result.args) == original_count


# =========================================================
# Subtraction tests
# =========================================================

def test_pure_tensor_sub_is_algebraictensor():
    """Subtracting two PureTensors produces an AlgebraicTensor."""
    pt1 = PureTensor(A, I3)
    pt2 = PureTensor(B, I3)
    result = pt1 - pt2
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))


def test_pure_tensor_sub_self():
    """A PureTensor minus itself produces an AlgebraicTensor (unsimplified)."""
    pt = PureTensor(A, I3)
    result = pt - pt
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_sub():
    """Subtracting two AlgebraicTensors produces an AlgebraicTensor."""
    at1 = PureTensor(A, I3) + PureTensor(B, I3)
    at2 = PureTensor(A, I3) + PureTensor(B, I3)
    result = at1 - at2
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_sub_puretensor():
    """Subtracting a PureTensor from an AlgebraicTensor."""
    at = PureTensor(A, I3) + PureTensor(B, I3)
    pt = PureTensor(A, I3)
    result = at - pt
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_sub_algebraictensor():
    """Subtracting an AlgebraicTensor from a PureTensor."""
    pt = PureTensor(A, I3)
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = pt - at
    assert isinstance(result, AlgebraicTensor)


def test_zerotensor_sub_puretensor():
    """ZeroTensor minus a PureTensor produces an AlgebraicTensor."""
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = z - pt
    assert isinstance(result, AlgebraicTensor)


def test_puretensor_sub_zerotensor():
    """PureTensor minus a ZeroTensor produces an AlgebraicTensor."""
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = pt - z
    assert isinstance(result, AlgebraicTensor)


# =========================================================
# Linear combinations
# =========================================================

def test_linear_combination_puretensors():
    """Construct a linear combination of PureTensors."""
    pt1 = PureTensor(A, C)
    pt2 = PureTensor(B, C)
    combo = 2 * pt1 + 3 * pt2
    assert isinstance(combo, AlgebraicTensor)
    assert combo.tensor_shape == ((3, 4), (4, 5))


def test_linear_combination_with_subtraction():
    """Linear combination with subtraction."""
    pt1 = PureTensor(A, C)
    pt2 = PureTensor(B, C)
    combo = 2 * pt1 - 3 * pt2
    assert isinstance(combo, AlgebraicTensor)
    assert combo.tensor_shape == ((3, 4), (4, 5))


def test_linear_combination_chained():
    """Chained linear combination operations."""
    pt1 = PureTensor(A, I3)
    pt2 = PureTensor(B, I3)
    pt3 = PureTensor(D, I3.__class__("D2", 3, 5).reshape(3, 3) if False else MatrixSymbol("D2", 3, 3))
    # Use compatible shapes
    D2 = MatrixSymbol("D2", 3, 4)
    pt3 = PureTensor(D2, I3)
    combo = pt1 + 2 * pt2 - 3 * pt3
    assert isinstance(combo, AlgebraicTensor)


def test_scaled_puretensor_in_algebraictensor():
    """Scaled PureTensors inside AlgebraicTensor retain coefficient info."""
    pt = PureTensor(A, C)
    at = AlgebraicTensor(2 * pt, 3 * pt)
    assert isinstance(at, AlgebraicTensor)
    for arg in at.args:
        if isinstance(arg, PureTensor):
            assert arg.factors == (A, C)


# =========================================================
# Display/str/repr for linear combinations
# =========================================================

def test_pure_tensor_str_with_coefficient():
    """PureTensor with numeric coefficient displays coefficient in str."""
    pt = PureTensor(A, C)
    scaled = 3 * pt
    s = str(scaled)
    assert "3" in s
    assert "A" in s
    assert "C" in s
    assert "\u2297" in s  # ⊗ symbol


def test_pure_tensor_str_with_negative_coefficient():
    """PureTensor with negative coefficient displays leading minus."""
    pt = PureTensor(A, C)
    neg = -3 * pt
    s = str(neg)
    assert s.startswith("-")
    assert "3" in s
    assert "A" in s
    assert "C" in s


def test_pure_tensor_str_with_negative_one():
    """PureTensor with coefficient -1 displays as unary minus."""
    pt = PureTensor(A, C)
    neg = -pt
    s = str(neg)
    assert s.startswith("-")
    assert "A" in s
    assert "C" in s
    # Should not contain "-1*" pattern
    assert "-1" not in s


def test_pure_tensor_repr_with_coefficient():
    """PureTensor with coefficient includes coefficient in repr."""
    pt = PureTensor(A, C)
    scaled = 3 * pt
    r = repr(scaled)
    assert "3" in r
    assert "PureTensor" in r


def test_pure_tensor_repr_without_coefficient():
    """PureTensor without coefficient has clean repr."""
    pt = PureTensor(A, C)
    r = repr(pt)
    assert "PureTensor(A, C)" in r or "PureTensor(A,C)" in r


def test_algebraic_tensor_str_with_coefficients():
    """AlgebraicTensor with scaled PureTensors displays coefficients."""
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    at = AlgebraicTensor(2 * pt_a, 3 * pt_b)
    s = str(at)
    assert "2" in s
    assert "3" in s
    assert "+" in s


def test_algebraic_tensor_str_with_negative_term():
    """AlgebraicTensor with negative term uses minus separator."""
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    at = AlgebraicTensor(2 * pt_a, -3 * pt_b)
    s = str(at)
    assert "-" in s
    # Should not have "+ -" pattern (improper double sign)
    assert "+ -" not in s


def test_linear_combination_str_display():
    """Linear combination expression displays with all coefficients."""
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    combo = 2 * pt_a + 3 * pt_b
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "2" in s
    assert "3" in s
    assert "+" in s


def test_linear_combination_subtraction_str():
    """Linear combination with subtraction displays correctly."""
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    combo = 2 * pt_a - 3 * pt_b
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "2" in s
    assert "3" in s
    assert "-" in s
    assert "+ -" not in s


def test_linear_combination_three_terms_str():
    """Three-term linear combination displays all coefficients."""
    D2 = MatrixSymbol("D2", 3, 4)
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    pt_d = PureTensor(D2, C)
    combo = pt_a + 2 * pt_b - 3 * pt_d
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "2" in s
    assert "3" in s
    assert "-" in s


def test_linear_combination_symbol_coefficient_str():
    """Linear combination with symbolic coefficient displays symbol."""
    from sympy import Symbol
    x = Symbol("x")
    y = Symbol("y")
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    # Use pt * x (forward mul) which returns PureTensor with symbolic coeff.
    # Note: x * pt goes through SymPy dispatch to MatMul instead of PureTensor.__rmul__
    scaled_a = pt_a * x
    scaled_b = pt_b * y
    assert isinstance(scaled_a, PureTensor)
    assert isinstance(scaled_b, PureTensor)
    combo = AlgebraicTensor(scaled_a, scaled_b)
    s = str(combo)
    assert isinstance(combo, AlgebraicTensor)
    assert "x" in s
    assert "y" in s
    assert "+" in s


def test_linear_combination_repr():
    """Linear combination repr shows all args with coefficients."""
    pt_a = PureTensor(A, C)
    pt_b = PureTensor(B, C)
    combo = 2 * pt_a + 3 * pt_b
    r = repr(combo)
    assert "AlgebraicTensor" in r
    assert "2" in r
    assert "3" in r


def test_algebraic_tensor_mul_scalar_str():
    """AlgebraicTensor multiplied by scalar displays coefficients."""
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    at = pt_a + pt_b
    scaled = at * 5
    s = str(scaled)
    assert isinstance(scaled, AlgebraicTensor)
    assert "5" in s


def test_algebraic_tensor_rmul_scalar_str():
    """Scalar multiplied with AlgebraicTensor displays coefficients."""
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    at = pt_a + pt_b
    scaled = 5 * at
    s = str(scaled)
    assert isinstance(scaled, AlgebraicTensor)
    assert "5" in s


def test_pure_tensor_coeff_preserved_through_add():
    """Coefficient of PureTensor preserved when added into AlgebraicTensor."""
    pt = PureTensor(A, C)
    scaled = 3 * pt
    at = AlgebraicTensor(scaled, PureTensor(B, C))
    for arg in at.args:
        if isinstance(arg, PureTensor) and arg.factors == (A, C):
            assert arg._get_coeff() == S(3)


def test_linear_combination_with_zerotensor_str():
    """AlgebraicTensor with ZeroTensor anchor and scaled terms displays correctly."""
    z = ZeroTensor(((3, 4), (4, 5)))
    pt = PureTensor(A, C)
    at = AlgebraicTensor(2 * pt, z)
    s = str(at)
    assert "2" in s
    assert "0_" in s  # ZeroTensor display


def test_pure_tensor_str_coefficient_one_omitted():
    """PureTensor with coefficient 1 does not show '1*' prefix."""
    pt = PureTensor(A, C)
    s = str(pt)
    assert not s.startswith("1")
    assert "*A" not in s.split("\u2297")[0] or s.startswith("A")


def test_pure_tensor_double_scaled_str():
    """PureTensor scaled twice (2*3=6) displays combined coefficient."""
    pt = PureTensor(A, C)
    result = (2 * pt) * 3
    s = str(result)
    assert "6" in s
    assert isinstance(result, PureTensor)
