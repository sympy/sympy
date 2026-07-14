from __future__ import annotations

from sympy.core import S, Symbol, symbols
from sympy.core.numbers import Integer
from sympy.matrices.expressions import MatrixSymbol
from sympy.printing.pretty import pretty as _pretty

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    algebraic_tensor_product,
    algebraic_zero_tensor,
    tensorsimplify,
)


def pretty(expr):
    """ASCII pretty-printing."""
    return _pretty(expr, use_unicode=False, wrap_line=False)


def upretty(expr):
    """Unicode pretty-printing."""
    return _pretty(expr, use_unicode=True, wrap_line=False)


# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------

def _make_matrices():
    """Create a standard set of matrix symbols for testing."""
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 3, 4)
    C = MatrixSymbol("C", 4, 5)
    D = MatrixSymbol("D", 4, 5)
    return {"A": A, "B": B, "C": C, "D": D}


# ---------------------------------------------------------------------------
# AlgebraicPureTensor pretty printing
# ---------------------------------------------------------------------------

def test_pretty_pure_tensor_simple_ascii():
    """ASCII: A x C for a two-factor pure tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    res = pretty(pt)
    assert "A" in res
    assert "C" in res
    assert "x" in res


def test_pretty_pure_tensor_simple_unicode():
    """Unicode: A with circled-times C for a two-factor pure tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    res = upretty(pt)
    assert "A" in res
    assert "C" in res
    assert "\u2a02" in res


def test_pretty_pure_tensor_single_factor_ascii():
    """Single factor renders without tensor product symbol."""
    mats = _make_matrices()
    A = mats["A"]
    pt = AlgebraicPureTensor(A)
    res = pretty(pt)
    assert "A" in res
    assert "x" not in res


def test_pretty_pure_tensor_single_factor_unicode():
    """Single factor renders without tensor product symbol."""
    mats = _make_matrices()
    A = mats["A"]
    pt = AlgebraicPureTensor(A)
    res = upretty(pt)
    assert "A" in res
    assert "\u2a02" not in res


def test_pretty_pure_tensor_with_coeff_ascii():
    """Numeric coefficient appears before the tensor product."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(2, A, C)
    res = pretty(pt)
    assert "2" in res
    assert "A" in res
    assert "C" in res
    assert "x" in res


def test_pretty_pure_tensor_with_coeff_unicode():
    """Numeric coefficient appears before the tensor product."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(2, A, C)
    res = upretty(pt)
    assert "2" in res
    assert "A" in res
    assert "C" in res
    assert "\u2a02" in res


def test_pretty_pure_tensor_negative_one_coeff_ascii():
    """Coefficient -1 renders as a leading minus sign."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(S.NegativeOne, A, C)
    res = pretty(pt)
    assert "-" in res
    assert "A" in res
    assert "C" in res


def test_pretty_pure_tensor_negative_one_coeff_unicode():
    """Coefficient -1 renders as a leading minus sign."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(S.NegativeOne, A, C)
    res = upretty(pt)
    assert "-" in res
    assert "A" in res
    assert "C" in res
    assert "\u2a02" in res


def test_pretty_pure_tensor_symbolic_coeff_ascii():
    """Symbolic coefficient renders before the tensor product."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")
    pt = AlgebraicPureTensor(x, A, C)
    res = pretty(pt)
    assert "x" in res
    assert "A" in res
    assert "C" in res


def test_pretty_pure_tensor_symbolic_coeff_unicode():
    """Symbolic coefficient renders before the tensor product."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    x = Symbol("x")
    pt = AlgebraicPureTensor(x, A, C)
    res = upretty(pt)
    assert "x" in res
    assert "A" in res
    assert "C" in res
    assert "\u2a02" in res


def test_pretty_pure_tensor_three_factors_ascii():
    """Three-factor tensor renders with two x delimiters."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    E = MatrixSymbol("E", 5, 3)
    pt = AlgebraicPureTensor(A, C, E)
    res = pretty(pt)
    assert "A" in res
    assert "C" in res
    assert "E" in res
    assert res.count("x") == 2


def test_pretty_pure_tensor_three_factors_unicode():
    """Three-factor tensor renders with two circled-times delimiters."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    E = MatrixSymbol("E", 5, 3)
    pt = AlgebraicPureTensor(A, C, E)
    res = upretty(pt)
    assert "A" in res
    assert "C" in res
    assert "E" in res
    assert res.count("\u2a02") == 2


def test_pretty_pure_tensor_coeff_one_no_prefix_ascii():
    """Coefficient of 1 does not produce a prefix."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    res = pretty(pt)
    assert "1" not in res


def test_pretty_pure_tensor_coeff_one_no_prefix_unicode():
    """Coefficient of 1 does not produce a prefix."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    res = upretty(pt)
    assert "1" not in res


# ---------------------------------------------------------------------------
# AlgebraicZeroTensor pretty printing
# ---------------------------------------------------------------------------

def test_pretty_zero_tensor_simple_ascii():
    """Zero tensor renders 0 with subscript showing shape."""
    zt = algebraic_zero_tensor(((3, 4),))
    res = pretty(zt)
    assert "0" in res
    assert "3" in res
    assert "4" in res


def test_pretty_zero_tensor_simple_unicode():
    """Zero tensor renders 0 with subscript showing shape."""
    zt = algebraic_zero_tensor(((3, 4),))
    res = upretty(zt)
    assert "0" in res
    assert "3" in res
    assert "4" in res
    assert "\u00d7" in res


def test_pretty_zero_tensor_two_factors_ascii():
    """Zero tensor with two factors shows both shapes."""
    zt = algebraic_zero_tensor(((3, 4), (4, 5)))
    res = pretty(zt)
    assert "0" in res
    assert "3" in res
    assert "4" in res
    assert "5" in res


def test_pretty_zero_tensor_two_factors_unicode():
    """Zero tensor with two factors shows both shapes."""
    zt = algebraic_zero_tensor(((3, 4), (4, 5)))
    res = upretty(zt)
    assert "0" in res
    assert "3" in res
    assert "4" in res
    assert "5" in res
    assert "\u00d7" in res


# ---------------------------------------------------------------------------
# AlgebraicTensor pretty printing
# ---------------------------------------------------------------------------

def test_pretty_tensor_sum_two_terms_ascii():
    """Two-term sum renders with + between terms."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(AlgebraicPureTensor(A, C), AlgebraicPureTensor(B, D))
    res = pretty(at)
    assert "A" in res
    assert "C" in res
    assert "B" in res
    assert "D" in res
    assert "+" in res


def test_pretty_tensor_sum_two_terms_unicode():
    """Two-term sum renders with + and circled-times."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(AlgebraicPureTensor(A, C), AlgebraicPureTensor(B, D))
    res = upretty(at)
    assert "A" in res
    assert "C" in res
    assert "B" in res
    assert "D" in res
    assert "+" in res
    assert "\u2a02" in res


def test_pretty_tensor_sum_with_negative_ascii():
    """Negative coefficient term renders with - sign."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(A, C),
        AlgebraicPureTensor(S.NegativeOne, B, D)
    )
    res = pretty(at)
    assert "-" in res
    assert "A" in res
    assert "B" in res


def test_pretty_tensor_sum_with_negative_unicode():
    """Negative coefficient term renders with - sign."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(A, C),
        AlgebraicPureTensor(S.NegativeOne, B, D)
    )
    res = upretty(at)
    assert "-" in res
    assert "A" in res
    assert "B" in res


def test_pretty_tensor_sum_with_coeff_ascii():
    """Term with numeric coefficient renders correctly."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(A, C),
        AlgebraicPureTensor(2, B, D)
    )
    res = pretty(at)
    assert "2" in res
    assert "+" in res
    assert "A" in res
    assert "B" in res


def test_pretty_tensor_sum_with_coeff_unicode():
    """Term with numeric coefficient renders correctly."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(A, C),
        AlgebraicPureTensor(2, B, D)
    )
    res = upretty(at)
    assert "2" in res
    assert "+" in res
    assert "A" in res
    assert "B" in res


def test_pretty_tensor_all_cancelled_ascii():
    """When all terms cancel via tensorsimplify, renders as zero tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(A, C),
        AlgebraicPureTensor(S.NegativeOne, A, C)
    )
    simplified = tensorsimplify(at)
    res = pretty(simplified)
    assert "0" in res


def test_pretty_tensor_all_cancelled_unicode():
    """When all terms cancel via tensorsimplify, renders as zero tensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(A, C),
        AlgebraicPureTensor(S.NegativeOne, A, C)
    )
    simplified = tensorsimplify(at)
    res = upretty(simplified)
    assert "0" in res


def test_pretty_tensor_single_term_unwraps_ascii():
    """Single-term AlgebraicTensor unwraps to PureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    at = AlgebraicTensor(AlgebraicPureTensor(A, C))
    res = pretty(at)
    assert "A" in res
    assert "C" in res
    assert "+" not in res


def test_pretty_tensor_single_term_unwraps_unicode():
    """Single-term AlgebraicTensor unwraps to PureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    at = AlgebraicTensor(AlgebraicPureTensor(A, C))
    res = upretty(at)
    assert "A" in res
    assert "C" in res
    assert "+" not in res


def test_pretty_tensor_leading_negative_ascii():
    """First term with negative coefficient renders with leading minus."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(S.NegativeOne, A, C),
        AlgebraicPureTensor(B, D)
    )
    res = pretty(at)
    assert "-" in res
    assert "A" in res
    assert "B" in res


def test_pretty_tensor_leading_negative_unicode():
    """First term with negative coefficient renders with leading minus."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(S.NegativeOne, A, C),
        AlgebraicPureTensor(B, D)
    )
    res = upretty(at)
    assert "-" in res
    assert "A" in res
    assert "B" in res


def test_pretty_tensor_leading_negative_coeff_ascii():
    """First term with negative numeric coefficient renders correctly."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(S.NegativeOne * 3, A, C),
        AlgebraicPureTensor(B, D)
    )
    res = pretty(at)
    assert "-" in res
    assert "3" in res
    assert "A" in res


def test_pretty_tensor_leading_negative_coeff_unicode():
    """First term with negative numeric coefficient renders correctly."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    at = AlgebraicTensor(
        AlgebraicPureTensor(S.NegativeOne * 3, A, C),
        AlgebraicPureTensor(B, D)
    )
    res = upretty(at)
    assert "-" in res
    assert "3" in res
    assert "A" in res


# ---------------------------------------------------------------------------
# Integration: arithmetic produces correct pretty output
# ---------------------------------------------------------------------------

def test_pretty_arithmetic_add_ascii():
    """A⊗C + B⊗D via arithmetic produces correct pretty output."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    expr = pt1 + pt2
    res = pretty(expr)
    assert "A" in res
    assert "B" in res


def test_pretty_arithmetic_add_unicode():
    """A⊗C + B⊗D via arithmetic produces correct pretty output."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    expr = pt1 + pt2
    res = upretty(expr)
    assert "A" in res
    assert "B" in res
    assert "+" in res


def test_pretty_arithmetic_sub_ascii():
    """A⊗C - B⊗D via arithmetic produces correct pretty output."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    expr = pt1 - pt2
    res = pretty(expr)
    assert "-" in res
    assert "A" in res
    assert "B" in res


def test_pretty_arithmetic_sub_unicode():
    """A⊗C - B⊗D via arithmetic produces correct pretty output."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    expr = pt1 - pt2
    res = upretty(expr)
    assert "-" in res
    assert "A" in res
    assert "B" in res


def test_pretty_scalar_mul_ascii():
    """2*(A⊗C) renders with coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    expr = 2 * pt
    res = pretty(expr)
    assert "2" in res
    assert "A" in res
    assert "C" in res


def test_pretty_scalar_mul_unicode():
    """2*(A⊗C) renders with coefficient."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    expr = 2 * pt
    res = upretty(expr)
    assert "2" in res
    assert "A" in res
    assert "C" in res


def test_pretty_negation_ascii():
    """-(A⊗C) renders with leading minus."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    expr = -pt
    res = pretty(expr)
    assert "-" in res
    assert "A" in res
    assert "C" in res


def test_pretty_negation_unicode():
    """-(A⊗C) renders with leading minus."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]
    pt = AlgebraicPureTensor(A, C)
    expr = -pt
    res = upretty(expr)
    assert "-" in res
    assert "A" in res
    assert "C" in res
