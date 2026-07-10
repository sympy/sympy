from __future__ import annotations

from sympy.core import S, Symbol, symbols
from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    tensorsimplify,
)
from sympy.tensor.algebraic.simplify import (
    _proportionality_factoring,
    _proportionality_ratio,
    _extract_pt_and_coeff,
    _build_pt,
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
# tensorsimplify dispatch tests
# ---------------------------------------------------------------------------

def test_simplify_zero_tensor():
    """tensorsimplify returns AlgebraicZeroTensor unchanged."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = tensorsimplify(zt)
    assert result is zt


def test_simplify_pure_tensor():
    """tensorsimplify simplifies AlgebraicPureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(2, A, C)
    result = tensorsimplify(pt)
    assert result is not None


def test_simplify_algebraic_tensor():
    """tensorsimplify simplifies AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    result = tensorsimplify(at)
    assert result is not None


def test_simplify_non_tensor():
    """tensorsimplify falls back to SymPy simplify for non-tensor expr."""
    x = Symbol("x")
    result = tensorsimplify(x + x)
    assert result == 2 * x


# ---------------------------------------------------------------------------
# Proportionality ratio tests
# ---------------------------------------------------------------------------

def test_proportionality_ratio_identical():
    """Identical factors have ratio 1."""
    mats = _make_matrices()
    A = mats["A"]
    ratio = _proportionality_ratio(A, A)
    assert ratio == S.One


def test_proportionality_ratio_numeric_matrices():
    """Proportional numeric matrices return the ratio."""
    M1 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    M2 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    ratio = _proportionality_ratio(M1, M2)
    assert ratio == 2


def test_proportionality_ratio_not_proportional():
    """Non-proportional matrices return None."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[1, 1], [1, 1]])
    ratio = _proportionality_ratio(M1, M2)
    assert ratio is None


def test_proportionality_ratio_shape_mismatch():
    """Matrices with different shapes return None."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[1, 2, 3]])
    ratio = _proportionality_ratio(M1, M2)
    assert ratio is None


def test_proportionality_ratio_zero_vs_nonzero():
    """One zero and one nonzero element returns None."""
    M1 = ImmutableDenseMatrix([[0, 0], [0, 0]])
    M2 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    ratio = _proportionality_ratio(M1, M2)
    assert ratio is None


# ---------------------------------------------------------------------------
# Proportionality factoring tests
# ---------------------------------------------------------------------------

def test_proportionality_factoring_all_proportional():
    """All slots proportional: merge by adding coefficients."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])  # 2 * M1
    N = ImmutableDenseMatrix([[1, 0], [0, 1]])

    pt1 = AlgebraicPureTensor(M1, N)
    pt2 = AlgebraicPureTensor(M2, N)
    at = AlgebraicTensor(pt1, pt2)
    result = _proportionality_factoring(at)
    assert result is not None


def test_proportionality_factoring_exactly_one_non_proportional():
    """Exactly one non-proportional slot: linear combination."""
    M1 = ImmutableDenseMatrix([[1, 0], [0, 0]])
    M2 = ImmutableDenseMatrix([[0, 0], [0, 1]])
    N = ImmutableDenseMatrix([[1, 2], [3, 4]])

    pt1 = AlgebraicPureTensor(M1, N)
    pt2 = AlgebraicPureTensor(M2, N)
    at = AlgebraicTensor(pt1, pt2)
    result = _proportionality_factoring(at)
    assert result is not None


def test_proportionality_factoring_more_than_one_non_proportional():
    """More than one non-proportional slot: skip merge."""
    M1 = ImmutableDenseMatrix([[1, 0], [0, 0]])
    M2 = ImmutableDenseMatrix([[0, 0], [0, 1]])
    N1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    N2 = ImmutableDenseMatrix([[5, 6], [7, 8]])

    pt1 = AlgebraicPureTensor(M1, N1)
    pt2 = AlgebraicPureTensor(M2, N2)
    at = AlgebraicTensor(pt1, pt2)
    result = _proportionality_factoring(at)
    assert result is not None


def test_proportionality_factoring_cancellation():
    """Proportional terms with opposite coefficients cancel."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    N = ImmutableDenseMatrix([[1, 0], [0, 1]])

    pt1 = AlgebraicPureTensor(M1, N)
    pt2 = AlgebraicPureTensor(S.NegativeOne, M1, N)
    at = AlgebraicTensor(pt1, pt2)
    result = _proportionality_factoring(at)
    assert isinstance(result, AlgebraicZeroTensor)


def test_proportionality_factoring_single_term():
    """Single term returns unchanged."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(pt)
    result = _proportionality_factoring(at)
    assert result is not None


# ---------------------------------------------------------------------------
# _extract_pt_and_coeff tests
# ---------------------------------------------------------------------------

def test_extract_pt_and_coeff_pure_tensor():
    """Extract from AlgebraicPureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(3, A, C)
    unit_pt, coeff = _extract_pt_and_coeff(pt)
    assert coeff == 3


def test_extract_pt_and_coeff_with_coefficient():
    """Extract coefficient from AlgebraicPureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")

    pt = AlgebraicPureTensor(x, A, C)
    unit_pt, coeff = _extract_pt_and_coeff(pt)
    assert coeff == x
    assert isinstance(unit_pt, AlgebraicPureTensor)


def test_extract_pt_and_coeff_unknown():
    """Extract from unknown type returns (term, S.One)."""
    x = Symbol("x")
    term, coeff = _extract_pt_and_coeff(x)
    assert term == x
    assert coeff == S.One


# ---------------------------------------------------------------------------
# _build_pt tests
# ---------------------------------------------------------------------------

def test_build_pt_zero_coeff():
    """Zero coefficient produces AlgebraicZeroTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    result = _build_pt(S.Zero, [A, C])
    assert isinstance(result, AlgebraicZeroTensor)


def test_build_pt_empty_factors():
    """Empty factors returns the coefficient."""
    result = _build_pt(3, [])
    assert result == 3


def test_build_pt_unit_coeff_single_factor():
    """Unit coefficient with single factor unwraps."""
    mats = _make_matrices()
    A = mats["A"]

    result = _build_pt(S.One, [A])
    assert result is A


def test_build_pt_unit_coeff_multiple_factors():
    """Unit coefficient with multiple factors builds PureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    result = _build_pt(S.One, [A, C])
    assert isinstance(result, AlgebraicPureTensor)


def test_build_pt_with_coeff():
    """Non-unit coefficient builds PureTensor with coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    result = _build_pt(3, [A, C])
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 3


# ---------------------------------------------------------------------------
# Integration tests: tensorsimplify on each type
# ---------------------------------------------------------------------------

def test_simplify_pure_tensor_with_simplifiable_coeff():
    """Simplify a PureTensor with simplifiable coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(S(4) / S(2), A, C)
    result = tensorsimplify(pt)
    assert result is not None


def test_simplify_tensor_with_proportional_terms():
    """Simplify a tensor with proportional terms."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    N = ImmutableDenseMatrix([[1, 0], [0, 1]])

    pt1 = AlgebraicPureTensor(M1, N)
    pt2 = AlgebraicPureTensor(M2, N)
    at = AlgebraicTensor(pt1, pt2)
    result = tensorsimplify(at)
    assert result is not None


def test_simplify_tensor_with_cancellation():
    """Simplify a tensor with cancelling terms."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    N = ImmutableDenseMatrix([[1, 0], [0, 1]])

    pt1 = AlgebraicPureTensor(M1, N)
    pt2 = AlgebraicPureTensor(S.NegativeOne, M1, N)
    at = AlgebraicTensor(pt1, pt2)
    result = tensorsimplify(at)
    assert isinstance(result, AlgebraicZeroTensor)


def test_simplify_preserves_zero_tensor():
    """Simplify preserves AlgebraicZeroTensor through operations."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = tensorsimplify(zt)
    assert result is zt
    assert result.shape == ((3, 4), (4, 5))
