from __future__ import annotations

from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)


def test_pure_tensor_neg():
    pt = AlgebraicPureTensor(A, C)
    neg = -pt
    assert isinstance(neg, AlgebraicPureTensor)
    assert neg.factors == pt.factors


def test_pure_tensor_neg_is_puretensor():
    pt = AlgebraicPureTensor(A, C)
    neg = -pt
    assert isinstance(neg, AlgebraicPureTensor)
    assert neg.factors == pt.factors
    assert neg.tensor_shape == pt.tensor_shape


def test_pure_tensor_neg_of_negated():
    pt = AlgebraicPureTensor(A, C)
    neg = -pt
    nneg = -neg
    assert nneg.factors == pt.factors


def test_pure_tensor_double_negation():
    pt = AlgebraicPureTensor(A, C)
    assert (-(-pt)).factors == pt.factors
    assert (-(-pt)).tensor_shape == pt.tensor_shape


def test_pure_tensor_mul_scalar():
    pt = AlgebraicPureTensor(A, C)
    result = pt * 3
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == pt.factors


def test_pure_tensor_rmul_scalar():
    pt = AlgebraicPureTensor(A, C)
    result = 3 * pt
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == pt.factors


def test_pure_tensor_mul_scalar_is_puretensor():
    pt = AlgebraicPureTensor(A, C)
    result = pt * 5
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, C)
    assert result.tensor_shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_scalar_is_puretensor():
    pt = AlgebraicPureTensor(A, C)
    result = 5 * pt
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_negative_scalar():
    pt = AlgebraicPureTensor(A, C)
    result = pt * (-3)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_one():
    pt = AlgebraicPureTensor(A, C)
    assert pt * S.One is pt
    assert S.One * pt is pt


def test_pure_tensor_mul_one_identity():
    pt = AlgebraicPureTensor(A, C)
    assert pt * S.One is pt
    assert S.One * pt is pt


def test_pure_tensor_mul_zero():
    pt = AlgebraicPureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_zero():
    pt = AlgebraicPureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_mul_zero_returns_zerotensor():
    pt = AlgebraicPureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_zero_returns_zerotensor():
    pt = AlgebraicPureTensor(A, C)
    result = 0 * pt
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_mul_coefficient_combines():
    pt = AlgebraicPureTensor(A, C)
    scaled = 2 * pt
    rescaled = scaled * 3
    assert isinstance(rescaled, AlgebraicPureTensor)
    assert rescaled.factors == (A, C)


def test_pure_tensor_mul_coefficient_cancellation():
    from sympy import Rational
    pt = AlgebraicPureTensor(A, C)
    scaled = pt * Rational(1, 3)
    rescaled = scaled * 3
    assert isinstance(rescaled, AlgebraicPureTensor)
    assert rescaled.factors == pt.factors


def test_pure_tensor_mul_symbol():
    from sympy import Symbol
    x = Symbol("x")
    pt = AlgebraicPureTensor(A, C)
    result = pt * x
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_puretensor():
    """* now performs composition (factor-wise matrix multiplication)."""
    # A(3,4), C(4,5), E(4,6), F(5,7) — inner dims match for composition
    E = MatrixSymbol("E", 4, 6)
    F = MatrixSymbol("F", 5, 7)
    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(E, F)
    result = pt1 * pt2
    assert isinstance(result, AlgebraicPureTensor)
    # Composition: (A@E) ⊗ (C@F)
    assert len(result.factors) == 2
    assert result.factors[0].shape == (3, 6)
    assert result.factors[1].shape == (4, 7)
    assert result.tensor_shape == ((3, 6), (4, 7))


def test_pure_tensor_mul_puretensor_with_coefficients():
    """Composition combines coefficients from both tensors."""
    E = MatrixSymbol("E", 4, 6)
    F = MatrixSymbol("F", 5, 7)
    pt1 = 2 * AlgebraicPureTensor(A, C)
    pt2 = 3 * AlgebraicPureTensor(E, F)
    result = pt1 * pt2
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 6
    assert len(result.factors) == 2
    assert result.factors[0].shape == (3, 6)
    assert result.factors[1].shape == (4, 7)
