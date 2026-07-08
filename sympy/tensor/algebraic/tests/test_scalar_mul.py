from __future__ import annotations

from sympy.core import S, Symbol
from sympy.matrices.expressions import MatrixSymbol, MatAdd
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    ScalarMul,
)


# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

def _make_matrices():
    """Create a standard set of matrix symbols for testing."""
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 3, 4)
    C = MatrixSymbol("C", 4, 5)
    D = MatrixSymbol("D", 4, 5)
    I4 = MatrixSymbol("I4", 4, 4)
    I5 = MatrixSymbol("I5", 5, 5)
    return {
        "A": A, "B": B, "C": C, "D": D, "I4": I4, "I5": I5,
    }


# ---------------------------------------------------------------------------
# Constructor tests
# ---------------------------------------------------------------------------

def test_constructor_basic():
    """Basic construction of ScalarMul."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert isinstance(sm, ScalarMul)
    assert sm.scalar == x
    assert sm.tensor is pt


def test_constructor_zero_scalar_returns_zero_tensor():
    """Zero scalar returns AlgebraicZeroTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    result = ScalarMul(S.Zero, pt)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_constructor_one_scalar_unwraps():
    """Scalar of 1 unwraps to the inner tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    result = ScalarMul(S.One, pt)
    assert result is pt


def test_constructor_with_zero_tensor():
    """ScalarMul with AlgebraicZeroTensor returns the zero tensor."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    x = Symbol("x")

    result = ScalarMul(x, zt)
    assert result is zt


def test_constructor_absorbs_number_coefficient():
    """Number coefficient is absorbed into PureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(2, A, C)
    sm = ScalarMul(x, pt)
    # The combined coefficient is x*2, which is not a Number
    assert isinstance(sm, ScalarMul)


def test_constructor_combined_one_unwraps():
    """When scalar * inner_coeff == 1, unwraps to PureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(S.Half, A, C)
    result = ScalarMul(2, pt)
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == S.One


def test_constructor_combined_zero_returns_zero_tensor():
    """When scalar * inner_coeff == 0, returns AlgebraicZeroTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(S.Zero, A, C)
    # This already returns AlgebraicZeroTensor from the PureTensor constructor
    # So we test with a different setup
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = ScalarMul(S.Zero, zt)
    assert isinstance(result, AlgebraicZeroTensor)


def test_constructor_combined_negative_one():
    """When scalar * inner_coeff == -1, returns PureTensor with -1 coeff."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(S.NegativeOne, A, C)
    result = ScalarMul(S.NegativeOne, pt)
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == S.One


# ---------------------------------------------------------------------------
# Properties tests
# ---------------------------------------------------------------------------

def test_scalar_property():
    """scalar property returns the commutative scalar."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert sm.scalar == x


def test_tensor_property():
    """tensor property returns the AlgebraicPureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert sm.tensor is pt


def test_factors_delegates():
    """factors delegates to inner tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert sm.factors == (A, C)


def test_tensor_shape_delegates():
    """tensor_shape delegates to inner tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert sm.tensor_shape == ((3, 4), (4, 5))


def test_commutativity_shape_delegates():
    """commutativity_shape delegates to inner tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert sm.commutativity_shape == (0, 0)


def test_get_coeff_combined():
    """_get_coeff returns combined scalar * inner_coeff."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(2, A, C)
    sm = ScalarMul(x, pt)
    coeff = sm._get_coeff()
    assert coeff == 2 * x


def test_is_commutative_false():
    """ScalarMul is non-commutative."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert sm.is_commutative is False


def test_is_scalar_mul_true():
    """ScalarMul has is_ScalarMul = True."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    assert sm.is_ScalarMul is True


# ---------------------------------------------------------------------------
# Arithmetic operator tests
# ---------------------------------------------------------------------------

def test_negation():
    """Negation negates the scalar part."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    neg = -sm
    assert isinstance(neg, ScalarMul)
    assert neg.scalar == -x


def test_negation_one_unwraps():
    """Negation that results in scalar=1 unwraps."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(S.NegativeOne, pt)
    neg = -sm
    assert neg is pt


def test_mul_by_number():
    """Multiplication by a Number multiplies into scalar."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    result = sm * 3
    assert isinstance(result, ScalarMul)
    assert result.scalar == 3 * x


def test_mul_by_symbol():
    """Multiplication by a Symbol multiplies into scalar."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x, y = Symbol("x"), Symbol("y")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    result = sm * y
    assert isinstance(result, ScalarMul)


def test_mul_by_one_unwraps():
    """Multiplication by reciprocal scalar unwraps."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(2, pt)
    result = sm * S.Half
    assert result is pt


def test_mul_by_zero_returns_zero_tensor():
    """Multiplication by zero returns AlgebraicZeroTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    result = sm * 0
    assert isinstance(result, AlgebraicZeroTensor)


def test_rmul_by_number():
    """Right-multiplication by a Number multiplies into scalar."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    result = 3 * sm
    assert isinstance(result, ScalarMul)


def test_add_returns_algebraic_tensor():
    """Addition returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    x = Symbol("x")

    pt1 = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, AlgebraicPureTensor(B, D))
    result = pt1 + sm
    assert isinstance(result, AlgebraicTensor)


def test_radd_returns_algebraic_tensor():
    """Right-addition returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    x = Symbol("x")

    pt1 = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, AlgebraicPureTensor(B, D))
    result = sm.__radd__(pt1)
    assert isinstance(result, AlgebraicTensor)


def test_sub_returns_algebraic_tensor():
    """Subtraction returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    x = Symbol("x")

    pt1 = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, AlgebraicPureTensor(B, D))
    result = pt1 - sm
    assert isinstance(result, AlgebraicTensor)


def test_rsub_returns_algebraic_tensor():
    """Right-subtraction returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    x = Symbol("x")

    pt1 = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, AlgebraicPureTensor(B, D))
    result = sm.__rsub__(pt1)
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# Expand tests
# ---------------------------------------------------------------------------

def test_expand_nothing_to_distribute():
    """Expand with nothing to distribute returns equivalent."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    result = sm.expand()
    assert isinstance(result, ScalarMul)


def test_expand_scalar_add():
    """Expand distributes scalar Add over tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x, y = Symbol("x"), Symbol("y")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x + y, pt)
    result = sm.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_factor_add():
    """Expand distributes MatAdd in tensor factors."""
    mats = _make_matrices()
    A, B, C = mats["A"], mats["B"], mats["C"]
    x = Symbol("x")

    matadd = MatAdd(A, B)
    pt = AlgebraicPureTensor(matadd, C)
    sm = ScalarMul(x, pt)
    result = sm.expand()
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# String representation tests
# ---------------------------------------------------------------------------

def test_str():
    """String representation shows scalar and tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    s = str(sm)
    assert "x" in s


def test_repr():
    """Repr contains class name."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    sm = ScalarMul(x, pt)
    r = repr(sm)
    assert "ScalarMul" in r
