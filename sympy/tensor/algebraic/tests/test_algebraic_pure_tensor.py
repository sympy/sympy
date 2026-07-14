from __future__ import annotations

from sympy.core import S, Symbol, symbols
from sympy.core.numbers import Number
from sympy.matrices.expressions import MatrixSymbol, MatAdd, MatMul
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    algebraic_tensor_product,
    compose_algebraic_pure_tensors,
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
    E = MatrixSymbol("E", 5, 3)
    F = MatrixSymbol("F", 5, 3)
    I3 = MatrixSymbol("I3", 3, 3)
    I4 = MatrixSymbol("I4", 4, 4)
    I5 = MatrixSymbol("I5", 5, 5)
    return {
        "A": A, "B": B, "C": C, "D": D, "E": E, "F": F,
        "I3": I3, "I4": I4, "I5": I5,
    }


# ---------------------------------------------------------------------------
# Constructor tests
# ---------------------------------------------------------------------------

def test_constructor_basic():
    """Basic construction of AlgebraicPureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    assert isinstance(pt, AlgebraicPureTensor)
    assert pt.factors == (A, C)
    assert pt.shape == ((3, 4), (4, 5))


def test_constructor_with_number_coefficient():
    """Construction with a Number coefficient stores coefficient internally."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(2, A, C)
    assert isinstance(pt, AlgebraicPureTensor)
    assert pt.coeff == 2
    assert pt.factors == (A, C)


def test_constructor_with_negative_coefficient():
    """Construction with a negative Number coefficient stores coefficient internally."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(S.NegativeOne, A, C)
    assert isinstance(pt, AlgebraicPureTensor)
    assert pt.coeff == S.NegativeOne
    assert pt.factors == (A, C)


def test_constructor_symbolic_coefficient():
    """Symbolic (non-Number) coefficient stores coefficient internally."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    result = AlgebraicPureTensor(x, A, C)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x
    assert result.factors == (A, C)


def test_constructor_single_factor_unwraps():
    """Single factor with coefficient S.One unwraps to bare factor."""
    mats = _make_matrices()
    A = mats["A"]

    result = AlgebraicPureTensor(A)
    assert result is A
    assert not isinstance(result, AlgebraicPureTensor)


def test_constructor_single_factor_with_coeff():
    """Single factor with Number coefficient stores coefficient internally."""
    mats = _make_matrices()
    A = mats["A"]

    pt = AlgebraicPureTensor(2, A)
    assert isinstance(pt, AlgebraicPureTensor)
    assert pt.coeff == 2
    assert pt.factors == (A,)


def test_constructor_zero_coefficient_returns_zero_tensor():
    """Zero coefficient returns AlgebraicZeroTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    result = AlgebraicPureTensor(S.Zero, A, C)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_constructor_empty_raises():
    """Empty constructor raises ValueError."""
    raises(ValueError, lambda: AlgebraicPureTensor())


def test_constructor_scalar_only_raises():
    """Only a scalar (no tensor factors) raises ValueError."""
    raises(ValueError, lambda: AlgebraicPureTensor(S.One))


def test_constructor_number_as_factor_raises():
    """A Number used as a tensor factor (not coefficient) raises TypeError."""
    mats = _make_matrices()
    A = mats["A"]
    raises(TypeError, lambda: AlgebraicPureTensor(A, S.One))


def test_constructor_no_shape_raises():
    """A factor without .shape raises TypeError."""
    x = Symbol("x")
    raises(TypeError, lambda: AlgebraicPureTensor(x))


# ---------------------------------------------------------------------------
# Properties tests
# ---------------------------------------------------------------------------

def test_factors_property():
    """factors property returns tensor factors (excludes coefficient)."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(2, A, C)
    assert isinstance(pt, AlgebraicPureTensor)
    assert pt.factors == (A, C)

    pt2 = AlgebraicPureTensor(A, C)
    assert pt2.factors == (A, C)


def test_coeff_property():
    """coeff property returns the coefficient (S.One if none)."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt_no_coeff = AlgebraicPureTensor(A, C)
    assert pt_no_coeff.coeff == S.One

    pt_with_coeff = AlgebraicPureTensor(2, A, C)
    assert pt_with_coeff.coeff == 2

    pt_symbolic = AlgebraicPureTensor(Symbol("x"), A, C)
    assert pt_symbolic.coeff == Symbol("x")


def test_shape():
    """shape returns tuple of per-factor shapes."""
    mats = _make_matrices()
    A, C, E = mats["A"], mats["C"], mats["E"]

    pt = AlgebraicPureTensor(A, C, E)
    assert pt.shape == ((3, 4), (4, 5), (5, 3))


def test_commutativity_pattern_symbolic():
    """commutativity_pattern is 0 for symbolic matrix factors."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    assert pt.commutativity_pattern == (0, 0)


def test_commutativity_pattern_numeric():
    """commutativity_pattern is 1 for numeric matrix factors."""
    numeric = ImmutableDenseMatrix([[1, 2], [3, 4]])
    A = MatrixSymbol("A", 2, 3)

    pt = AlgebraicPureTensor(numeric, A)
    assert pt.commutativity_pattern == (1, 0)


def test_commutativity_pattern_all_numeric():
    """commutativity_pattern is all-1s when all factors are numeric."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[1, 0], [0, 1], [1, 1]])

    pt = AlgebraicPureTensor(M1, M2)
    assert pt.commutativity_pattern == (1, 1)


def test_is_commutative_false():
    """AlgebraicPureTensor is non-commutative."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    assert pt.is_commutative is False


def test_is_mul_false():
    """AlgebraicPureTensor has is_Mul = False."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    assert pt.is_Mul is False


# ---------------------------------------------------------------------------
# Arithmetic operator tests
# ---------------------------------------------------------------------------

def test_negation_basic():
    """Negation negates the coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(2, A, C)
    neg = -pt
    assert isinstance(neg, AlgebraicPureTensor)
    assert neg.coeff == S.NegativeOne * 2


def test_negation_no_coeff():
    """Negation of unit-coefficient tensor returns coefficient -1."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    neg = -pt
    assert isinstance(neg, AlgebraicPureTensor)
    assert neg.coeff == S.NegativeOne


def test_negation_double():
    """Double negation returns original."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(2, A, C)
    result = -(-pt)
    assert result.coeff == pt.coeff
    assert result.factors == pt.factors


def test_mul_commulative_number():
    """Multiplication by a Number scales the coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(3, A, C)
    result = pt * 2
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6


def test_mul_commulative_one():
    """Multiplication by 1 returns self."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    assert pt * 1 is pt


def test_mul_commulative_zero():
    """Multiplication by 0 returns AlgebraicZeroTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    result = pt * 0
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_mul_symbolic_returns_scaled():
    """Multiplication by a symbolic scales the coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    result = pt * x
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x


def test_rmul_commulative_number():
    """Right-multiplication by a Number scales the coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(3, A, C)
    result = 2 * pt
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6


def test_rmul_symbolic_returns_scaled():
    """Right-multiplication by a symbolic scales the coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    result = x * pt
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x


def test_add_returns_algebraic_tensor():
    """Addition of two PureTensors returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    result = pt1 + pt2
    assert isinstance(result, AlgebraicTensor)


def test_add_with_zero_tensor():
    """Addition with same-shape AlgebraicZeroTensor returns self."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = pt + zt
    assert result is pt


def test_sub_returns_algebraic_tensor():
    """Subtraction returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    result = pt1 - pt2
    assert isinstance(result, AlgebraicTensor)


def test_radd_returns_algebraic_tensor():
    """Right-addition returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    result = pt2.__radd__(pt1)
    assert isinstance(result, AlgebraicTensor)


def test_rsub_returns_algebraic_tensor():
    """Right-subtraction returns AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    result = pt2.__rsub__(pt1)
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# Composition tests
# ---------------------------------------------------------------------------

def test_compose_pure_tensors_basic():
    """Basic composition of two AlgebraicPureTensors."""
    mats = _make_matrices()
    A, C, I4, I5 = mats["A"], mats["C"], mats["I4"], mats["I5"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(I4, I5)
    result = compose_algebraic_pure_tensors(pt1, pt2)
    assert isinstance(result, AlgebraicPureTensor)
    assert len(result.factors) == 2


def test_compose_with_identity():
    """Composition with identity matrices returns equivalent tensor."""
    mats = _make_matrices()
    A, C, I4, I5 = mats["A"], mats["C"], mats["I4"], mats["I5"]

    pt = AlgebraicPureTensor(A, C)
    pt_id = AlgebraicPureTensor(I4, I5)
    result = compose_algebraic_pure_tensors(pt, pt_id)
    assert isinstance(result, AlgebraicPureTensor)


def test_compose_coefficients_combine():
    """Composition combines coefficients from both sides."""
    mats = _make_matrices()
    A, C, I4, I5 = mats["A"], mats["C"], mats["I4"], mats["I5"]

    pt1 = AlgebraicPureTensor(2, A, C)
    pt2 = AlgebraicPureTensor(3, I4, I5)
    result = compose_algebraic_pure_tensors(pt1, pt2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6


def test_compose_factor_count_mismatch_raises():
    """Composition with different factor counts raises ValueError."""
    mats = _make_matrices()
    A, C, E = mats["A"], mats["C"], mats["E"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(A, C, E)
    raises(ValueError, lambda: compose_algebraic_pure_tensors(pt1, pt2))


def test_compose_inner_dim_mismatch_raises():
    """Composition with incompatible inner dimensions raises ValueError."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(C, A)
    raises(ValueError, lambda: compose_algebraic_pure_tensors(pt1, pt2))


def test_compose_bare_matrix():
    """Composition accepts bare matrix objects as single-factor tensors."""
    mats = _make_matrices()
    A, I4 = mats["A"], mats["I4"]

    result = compose_algebraic_pure_tensors(A, I4)
    assert result is not None


def test_compose_with_zero_tensor():
    """Composition with AlgebraicZeroTensor returns zero tensor with composed shape."""
    mats = _make_matrices()
    A, C, I4, I5 = mats["A"], mats["C"], mats["I4"], mats["I5"]

    pt = AlgebraicPureTensor(A, C)  # shape ((3,4), (4,5))
    zt = AlgebraicZeroTensor(((4, 5), (5, 3)))
    # Composed shape: (3,3) from A(3,4)*zt_slot0(4,5)->(3,5) wait no...
    # pt factors: A(3,4), C(4,5); zt shape: ((4,5), (5,3))
    # Composed: (A_rows=3, zt_col0=5) and (C_rows=4, zt_col1=3) -> ((3,5), (4,3))
    result = compose_algebraic_pure_tensors(pt, zt)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 5), (4, 3))

    # Zero on the left
    zt2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt2 = AlgebraicPureTensor(I4, I5)  # shape ((4,4), (5,5))
    # Composed: (zt_row0=3, I4_cols=4) and (zt_row1=4, I5_cols=5) -> ((3,4), (4,5))
    result2 = compose_algebraic_pure_tensors(zt2, pt2)
    assert isinstance(result2, AlgebraicZeroTensor)
    assert result2.shape == ((3, 4), (4, 5))

    # Both zero tensors
    zt3 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    zt4 = AlgebraicZeroTensor(((4, 5), (5, 3)))
    result3 = compose_algebraic_pure_tensors(zt3, zt4)
    assert isinstance(result3, AlgebraicZeroTensor)
    assert result3.shape == ((3, 5), (4, 3))


def test_compose_invalid_type_raises():
    """Composition with invalid type raises TypeError."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    raises(TypeError, lambda: compose_algebraic_pure_tensors(pt, Symbol("x")))


# ---------------------------------------------------------------------------
# Scalar multiplication vs composition dispatch
# ---------------------------------------------------------------------------

def test_mul_dispatch_commulative():
    """Commutative operand goes to scalar multiplication path."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(A, C)
    result = pt * x
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x


def test_mul_dispatch_noncommutative():
    """Non-commutative operand goes to composition path."""
    mats = _make_matrices()
    A, C, I4, I5 = mats["A"], mats["C"], mats["I4"], mats["I5"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(I4, I5)
    result = pt1 * pt2
    assert isinstance(result, AlgebraicPureTensor)


# ---------------------------------------------------------------------------
# Expand tests
# ---------------------------------------------------------------------------

def test_expand_no_add():
    """Expand with no Add factors returns equivalent tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    result = pt.expand()
    assert isinstance(result, AlgebraicPureTensor)


def test_expand_with_matadd():
    """Expand distributes over MatAdd in factors."""
    mats = _make_matrices()
    A, B, C = mats["A"], mats["B"], mats["C"]

    matadd = MatAdd(A, B)
    pt = AlgebraicPureTensor(matadd, C)
    result = pt.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_coefficient_preserved():
    """Expand preserves the coefficient through distribution."""
    mats = _make_matrices()
    A, B, C = mats["A"], mats["B"], mats["C"]

    matadd = MatAdd(A, B)
    pt = AlgebraicPureTensor(2, matadd, C)
    result = pt.expand()
    assert isinstance(result, AlgebraicTensor)
    for term in result.args:
        if isinstance(term, AlgebraicPureTensor):
            assert term.coeff == 2


# ---------------------------------------------------------------------------
# String representation tests
# ---------------------------------------------------------------------------

def test_str_no_coeff():
    """String representation without coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    s = str(pt)
    assert "A" in s and "C" in s


def test_str_with_coeff():
    """String representation with coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(2, A, C)
    s = str(pt)
    assert "2" in s


def test_str_negative_coeff():
    """String representation with negative coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(S.NegativeOne, A, C)
    s = str(pt)
    assert s.startswith("-")


def test_repr():
    """srepr contains class name (repr == str per SymPy convention)."""
    from sympy.printing.repr import srepr
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    r = srepr(pt)
    assert "AlgebraicPureTensor" in r


# ---------------------------------------------------------------------------
# Convenience function tests
# ---------------------------------------------------------------------------

def test_algebraic_tensor_product():
    """algebraic_tensor_product is a convenience constructor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    result = algebraic_tensor_product(A, C)
    assert isinstance(result, AlgebraicPureTensor)


# ---------------------------------------------------------------------------
# has_zero_term test
# ---------------------------------------------------------------------------

def test_has_zero_term_false():
    """has_zero_term returns False for nonzero AlgebraicPureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    assert pt.has_zero_term() is False
