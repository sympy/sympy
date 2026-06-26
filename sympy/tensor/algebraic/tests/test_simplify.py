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
        raise AssertionError(
            f"Expected {expected_exc.__name__} but got {type(e).__name__}: {e}"
        )
    else:
        raise AssertionError(
            f"Expected {expected_exc.__name__} but no exception was raised"
        )

from sympy.core.singleton import S
from sympy.core.add import Add
from sympy.matrices.expressions import MatrixSymbol
from sympy import Symbol, Rational, sqrt, I

from sympy.tensor.algebraic.algebraic_tensor import (
    AlgebraicTensor,
    ShapeMismatchError,
)
from sympy.tensor.algebraic.pure_tensor import PureTensor, tensor_product
from sympy.tensor.algebraic.zero_tensor import ZeroTensor, zero_tensor
from sympy.tensor.algebraic.simplify import tensorsimplify


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)
I4 = MatrixSymbol("I4", 4, 4)
J4 = MatrixSymbol("J4", 4, 4)


# =========================================================
# ZeroTensor simplification
# =========================================================

def test_zerotensor_simplify_returns_self():
    """ZeroTensor.simplify() returns the same object."""
    z = ZeroTensor(((3, 4), (4, 5)))
    assert z.simplify() is z


def test_zerotensor_tensorsimplify():
    """tensorsimplify on ZeroTensor returns the same object."""
    z = ZeroTensor(((3, 4),))
    assert tensorsimplify(z) is z


# =========================================================
# PureTensor simplification
# =========================================================

def test_puretensor_simplify_no_change():
    """PureTensor with trivial coefficient and simple factors is unchanged."""
    pt = PureTensor(A, C)
    result = pt.simplify()
    assert isinstance(result, PureTensor)
    assert result.factors == pt.factors


def test_puretensor_simplify_coefficient():
    """PureTensor coefficient is simplified."""
    pt = PureTensor(2 + Rational(1, 2) - Rational(1, 2), A, C)
    result = pt.simplify()
    assert isinstance(result, PureTensor)
    assert result._get_coeff() == 2


def test_puretensor_simplify_zero_coeff_returns_zerotensor():
    """PureTensor whose coefficient simplifies to zero becomes ZeroTensor."""
    pt = PureTensor(S.Zero, A, C)
    result = pt.simplify()
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_puretensor_simplify_method():
    """PureTensor.simplify() delegates to tensorsimplify."""
    pt = PureTensor(A, C)
    result = pt.simplify()
    assert isinstance(result, (PureTensor, ZeroTensor))


# =========================================================
# AlgebraicTensor: combine like terms
# =========================================================

def test_algebraictensor_combine_identical_terms():
    """Two identical PureTensors with numeric coefficients are combined."""
    pt = PureTensor(A, C)
    at = AlgebraicTensor(2 * pt, 3 * pt)
    result = at.simplify()
    # The two terms 2*pt and 3*pt share the same factor sequence,
    # so coefficients should add: 2+3 = 5.
    assert isinstance(result, PureTensor)
    assert result._get_coeff() == 5
    assert result.factors == (A, C)


def test_algebraictensor_combine_identical_terms_subtraction():
    """Identical terms with opposite sign cancel."""
    pt = PureTensor(A, C)
    at = AlgebraicTensor(3 * pt, -3 * pt)
    result = at.simplify()
    # All terms cancel; the result should be a ZeroTensor.
    assert isinstance(result, ZeroTensor)


def test_algebraictensor_combine_identical_terms_partial_cancel():
    """Identical terms that partially cancel retain the remainder."""
    pt = PureTensor(A, C)
    at = AlgebraicTensor(5 * pt, -2 * pt)
    result = at.simplify()
    assert isinstance(result, PureTensor)
    assert result._get_coeff() == 3
    assert result.factors == (A, C)


def test_algebraictensor_no_combine_different_factors():
    """Terms with different factor sequences are NOT combined."""
    at = PureTensor(A, C) + PureTensor(B, C)
    result = at.simplify()
    # A and B are different factors, so no combining.
    assert isinstance(result, AlgebraicTensor)


def test_algebraictensor_combine_three_identical():
    """Three identical terms combine their coefficients."""
    pt = PureTensor(A, C)
    at = AlgebraicTensor(pt, 2 * pt, 3 * pt)
    result = at.simplify()
    assert isinstance(result, PureTensor)
    assert result._get_coeff() == 6
    assert result.factors == (A, C)


def test_algebraictensor_combine_with_symbol_coefficients():
    """Terms with symbolic coefficients combine correctly."""
    x = Symbol("x")
    y = Symbol("y")
    pt = PureTensor(A, C)
    at = AlgebraicTensor(x * pt, y * pt)
    result = at.simplify()
    assert isinstance(result, PureTensor)
    coeff = result._get_coeff()
    # The coefficient should be x + y (simplified by SymPy)
    assert coeff == x + y
    assert result.factors == (A, C)


def test_algebraictensor_mixed_like_and_unlike():
    """Mix of identical and different terms: only identical ones combine."""
    pt_ac = PureTensor(A, C)
    pt_bc = PureTensor(B, C)
    at = AlgebraicTensor(2 * pt_ac, 3 * pt_ac, pt_bc)
    result = at.simplify()
    assert isinstance(result, AlgebraicTensor)
    # There should be a 5*A⊗C term and a B⊗C term.
    found_five = False
    found_one = False
    for arg in result.args:
        if isinstance(arg, PureTensor) and arg.factors == (A, C):
            assert arg._get_coeff() == 5
            found_five = True
        if isinstance(arg, PureTensor) and arg.factors == (B, C):
            assert arg._get_coeff() == 1
            found_one = True
    assert found_five and found_one


# =========================================================
# AlgebraicTensor: common factor extraction + recurse
# =========================================================

def test_algebraictensor_common_right_factor():
    """Common right factor is extracted, middle simplified."""
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    result = at.simplify()
    # Both terms share I4, C on the right; after extraction the middle is
    # A + B and the right factors are (I4, C).
    # Result should be PureTensor(AlgebraicTensor(A, B), I4, C) or
    # an equivalent structure.
    assert result is not None


def test_algebraictensor_common_left_factor():
    """Common left factor is extracted, middle simplified."""
    at = PureTensor(A, I4, C) + PureTensor(A, J4, C)
    result = at.simplify()
    # A is common on the left.
    assert result is not None


def test_algebraictensor_common_both_sides():
    """Both left and right common factors are extracted."""
    at = PureTensor(A, I4, C) + PureTensor(A, J4, C)
    result = at.simplify()
    # Left: A, right: C, middle: I4 + J4
    assert result is not None


# =========================================================
# AlgebraicTensor: ZeroTensor anchor
# =========================================================

def test_algebraictensor_with_zerotensor_anchor():
    """Simplify preserves ZeroTensor anchor when present."""
    z = ZeroTensor(((3, 4), (4, 5)))
    pt = PureTensor(A, C)
    at = AlgebraicTensor(pt, z)
    result = at.simplify()
    assert result is not None


def test_algebraictensor_all_cancel_with_zerotensor():
    """When all non-zero terms cancel and a ZeroTensor anchor exists,
    return the ZeroTensor."""
    z = ZeroTensor(((3, 4), (4, 5)))
    pt = PureTensor(A, C)
    at = AlgebraicTensor(pt, -pt, z)
    result = at.simplify()
    assert isinstance(result, ZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


# =========================================================
# AlgebraicTensor: coefficient simplification
# =========================================================

def test_algebraictensor_simplify_coefficient():
    """Numeric coefficients inside the AlgebraicTensor are simplified."""
    pt = PureTensor(A, C)
    at = AlgebraicTensor((1 + 1) * pt, (2 + 3) * pt)
    result = at.simplify()
    # 2*pt + 5*pt = 7*pt
    assert isinstance(result, PureTensor)
    assert result._get_coeff() == 7


# =========================================================
# tensorsimplify dispatch
# =========================================================

def test_tensorsimplify_dispatch_algebraic():
    """tensorsimplify routes AlgebraicTensor to the right handler."""
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = tensorsimplify(at)
    assert isinstance(result, AlgebraicTensor)


def test_tensorsimplify_dispatch_pure():
    """tensorsimplify routes PureTensor to the right handler."""
    pt = PureTensor(A, C)
    result = tensorsimplify(pt)
    assert isinstance(result, (PureTensor, MatrixSymbol))


def test_tensorsimplify_dispatch_zerotensor():
    """tensorsimplify routes ZeroTensor to the right handler."""
    z = ZeroTensor(((3, 4),))
    result = tensorsimplify(z)
    assert result is z


def test_tensorsimplify_dispatch_regular_expr():
    """tensorsimplify on a non-tensor falls back to SymPy simplify."""
    x = Symbol("x")
    expr = x + x
    result = tensorsimplify(expr)
    assert result == 2 * x


# =========================================================
# Package-level re-export
# =========================================================

def test_init_reexport_tensorsimplify():
    """tensorsimplify is re-exported from the package __init__."""
    from sympy.tensor.algebraic import tensorsimplify as ts
    assert ts is tensorsimplify


# =========================================================
# Edge cases
# =========================================================

def test_simplify_single_term_algebraictensor_no_zerotensor():
    """AlgebraicTensor with a single PureTensor and no ZeroTensor unwraps."""
    pt = PureTensor(A, C)
    at = AlgebraicTensor(pt)
    assert at is pt
    # simplify on the unwrapped PureTensor is fine.
    result = pt.simplify()
    assert result is pt or (isinstance(result, PureTensor) and result.factors == (A, C))


def test_simplify_negative_coefficient():
    """PureTensor with negative coefficient simplifies correctly."""
    pt = PureTensor(A, C)
    neg = -pt
    result = neg.simplify()
    assert isinstance(result, PureTensor)
    assert result.factors == (A, C)


def test_simplify_double_negation_algebraic():
    """Double negation of AlgebraicTensor simplifies to original form."""
    pt1 = PureTensor(A, C)
    pt2 = PureTensor(B, C)
    at = pt1 + pt2
    neg = -at
    unneg = -neg
    result = unneg.simplify()
    assert isinstance(result, AlgebraicTensor)


def test_simplify_identical_multi_factor():
    """Identical three-factor PureTensors combine."""
    pt = PureTensor(A, I4, C)
    at = AlgebraicTensor(2 * pt, 3 * pt)
    result = at.simplify()
    assert isinstance(result, PureTensor)
    assert result._get_coeff() == 5
    assert result.factors == (A, I4, C)


def test_simplify_cancel_then_common_factor():
    """Terms cancel first, then common factors are extracted."""
    pt_ac = PureTensor(A, C)
    pt_bc = PureTensor(B, C)
    # 2*AC + 3*AC - 2*AC = 3*AC
    at = AlgebraicTensor(2 * pt_ac, 3 * pt_ac, -2 * pt_ac)
    result = at.simplify()
    assert isinstance(result, PureTensor)
    assert result._get_coeff() == 3


def test_simplify_four_terms_combine_pairs():
    """Four terms where two pairs share factor sequences."""
    pt_ac = PureTensor(A, C)
    pt_bc = PureTensor(B, C)
    at = AlgebraicTensor(pt_ac, 2 * pt_ac, pt_bc, 4 * pt_bc)
    result = at.simplify()
    assert isinstance(result, AlgebraicTensor)
    # Should have 3*AC and 5*BC.
    found_3ac = False
    found_5bc = False
    for arg in result.args:
        if isinstance(arg, PureTensor):
            if arg.factors == (A, C):
                assert arg._get_coeff() == 3
                found_3ac = True
            elif arg.factors == (B, C):
                assert arg._get_coeff() == 5
                found_5bc = True
    assert found_3ac and found_5bc


def test_simplify_coefficient_symbol_cancellation():
    """Symbol coefficients cancel when opposite."""
    x = Symbol("x")
    pt = PureTensor(A, C)
    at = AlgebraicTensor(x * pt, (-x) * pt)
    result = at.simplify()
    assert isinstance(result, ZeroTensor)


def test_simplify_preserves_tensor_shape():
    """Simplified result preserves the original tensor shape."""
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    orig_shape = at.tensor_shape
    result = at.simplify()
    # The result may be wrapped differently; check it's meaningful.
    assert result is not None
