from __future__ import annotations

from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.pure_tensor import PureTensor
from sympy.tensor.algebraic.zero_tensor import ZeroTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)


def test_pure_tensor_neg():
    pt = PureTensor(A, C)
    neg = -pt
    assert isinstance(neg, PureTensor)
    assert neg.factors == pt.factors


def test_pure_tensor_neg_is_puretensor():
    pt = PureTensor(A, C)
    neg = -pt
    assert isinstance(neg, PureTensor)
    assert neg.factors == pt.factors
    assert neg.tensor_shape == pt.tensor_shape


def test_pure_tensor_neg_of_negated():
    pt = PureTensor(A, C)
    neg = -pt
    nneg = -neg
    assert nneg.factors == pt.factors


def test_pure_tensor_double_negation():
    pt = PureTensor(A, C)
    assert (-(-pt)).factors == pt.factors
    assert (-(-pt)).tensor_shape == pt.tensor_shape


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


def test_pure_tensor_mul_scalar_is_puretensor():
    pt = PureTensor(A, C)
    result = pt * 5
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)
    assert result.tensor_shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_scalar_is_puretensor():
    pt = PureTensor(A, C)
    result = 5 * pt
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_negative_scalar():
    pt = PureTensor(A, C)
    result = pt * (-3)
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_one():
    pt = PureTensor(A, C)
    assert pt * S.One is pt
    assert S.One * pt is pt


def test_pure_tensor_mul_one_identity():
    pt = PureTensor(A, C)
    assert pt * S.One is pt
    assert S.One * pt is pt


def test_pure_tensor_mul_zero():
    pt = PureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_zero():
    pt = PureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_mul_zero_returns_zerotensor():
    pt = PureTensor(A, C)
    result = pt * S.Zero
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_rmul_zero_returns_zerotensor():
    pt = PureTensor(A, C)
    result = 0 * pt
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_pure_tensor_mul_coefficient_combines():
    pt = PureTensor(A, C)
    scaled = 2 * pt
    rescaled = scaled * 3
    assert isinstance(rescaled, PureTensor)
    assert rescaled.factors == (A, C)


def test_pure_tensor_mul_coefficient_cancellation():
    from sympy import Rational
    pt = PureTensor(A, C)
    scaled = pt * Rational(1, 3)
    rescaled = scaled * 3
    assert isinstance(rescaled, PureTensor)
    assert rescaled.factors == pt.factors


def test_pure_tensor_mul_symbol():
    from sympy import Symbol
    x = Symbol("x")
    pt = PureTensor(A, C)
    result = pt * x
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)


def test_pure_tensor_mul_puretensor():
    pt1 = PureTensor(A, C)
    pt2 = PureTensor(C, D)
    result = pt1 * pt2
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C, C, D)
    assert result.tensor_shape == ((3, 4), (4, 5), (4, 5), (3, 5))


def test_pure_tensor_mul_puretensor_with_coefficients():
    pt1 = 2 * PureTensor(A, C)
    pt2 = 3 * PureTensor(C, D)
    result = pt1 * pt2
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C, C, D)
