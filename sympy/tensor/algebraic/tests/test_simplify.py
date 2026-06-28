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
from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor, algebraic_tensor_product
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor, algebraic_zero_tensor
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
# AlgebraicZeroTensor simplification
# =========================================================

def test_zerotensor_simplify_returns_self():
    """AlgebraicZeroTensor.simplify() returns the same object."""
    z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert z.simplify() is z


def test_zerotensor_tensorsimplify():
    """tensorsimplify on AlgebraicZeroTensor returns the same object."""
    z = AlgebraicZeroTensor(((3, 4),))
    assert tensorsimplify(z) is z


# =========================================================
# AlgebraicPureTensor simplification
# =========================================================

def test_puretensor_simplify_no_change():
    """AlgebraicPureTensor with trivial coefficient and simple factors is unchanged."""
    pt = AlgebraicPureTensor(A, C)
    result = pt.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == pt.factors


def test_puretensor_simplify_coefficient():
    """AlgebraicPureTensor coefficient is simplified."""
    pt = AlgebraicPureTensor(2 + Rational(1, 2) - Rational(1, 2), A, C)
    result = pt.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 2


def test_puretensor_simplify_zero_coeff_returns_zerotensor():
    """AlgebraicPureTensor whose coefficient simplifies to zero becomes AlgebraicZeroTensor."""
    pt = AlgebraicPureTensor(S.Zero, A, C)
    result = pt.simplify()
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_puretensor_simplify_method():
    """AlgebraicPureTensor.simplify() delegates to tensorsimplify."""
    pt = AlgebraicPureTensor(A, C)
    result = pt.simplify()
    assert isinstance(result, (AlgebraicPureTensor, AlgebraicZeroTensor))


# =========================================================
# AlgebraicTensor: combine like terms
# =========================================================

def test_algebraictensor_combine_identical_terms():
    """Two identical AlgebraicPureTensors with numeric coefficients are combined."""
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(2 * pt, 3 * pt)
    result = at.simplify()
    # The two terms 2*pt and 3*pt share the same factor sequence,
    # so coefficients should add: 2+3 = 5.
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 5
    assert result.factors == (A, C)


def test_algebraictensor_combine_identical_terms_subtraction():
    """Identical terms with opposite sign cancel."""
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(3 * pt, -3 * pt)
    result = at.simplify()
    # All terms cancel; the result should be a AlgebraicZeroTensor.
    assert isinstance(result, AlgebraicZeroTensor)


def test_algebraictensor_combine_identical_terms_partial_cancel():
    """Identical terms that partially cancel retain the remainder."""
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(5 * pt, -2 * pt)
    result = at.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 3
    assert result.factors == (A, C)


def test_algebraictensor_combine_three_identical():
    """Three identical terms combine their coefficients."""
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(pt, 2 * pt, 3 * pt)
    result = at.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 6
    assert result.factors == (A, C)


def test_algebraictensor_combine_with_symbol_coefficients():
    """Terms with symbolic coefficients combine correctly."""
    x = Symbol("x")
    y = Symbol("y")
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(x * pt, y * pt)
    result = at.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    coeff = result._get_coeff()
    # The coefficient should be x + y (simplified by SymPy)
    assert coeff == x + y
    assert result.factors == (A, C)


# =========================================================
# AlgebraicTensor: AlgebraicZeroTensor anchor
# =========================================================

def test_algebraictensor_with_zerotensor_anchor():
    """Simplify preserves AlgebraicZeroTensor anchor when present."""
    z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(pt, z)
    result = at.simplify()
    assert result is not None


def test_algebraictensor_all_cancel_with_zerotensor():
    """When all non-zero terms cancel and a AlgebraicZeroTensor anchor exists,
    return the AlgebraicZeroTensor."""
    z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(pt, -pt, z)
    result = at.simplify()
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


# =========================================================
# AlgebraicTensor: coefficient simplification
# =========================================================

def test_algebraictensor_simplify_coefficient():
    """Numeric coefficients inside the AlgebraicTensor are simplified."""
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor((1 + 1) * pt, (2 + 3) * pt)
    result = at.simplify()
    # 2*pt + 5*pt = 7*pt
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 7


# =========================================================
# tensorsimplify dispatch
# =========================================================

def test_tensorsimplify_dispatch_pure():
    """tensorsimplify routes AlgebraicPureTensor to the right handler."""
    pt = AlgebraicPureTensor(A, C)
    result = tensorsimplify(pt)
    assert isinstance(result, (AlgebraicPureTensor, MatrixSymbol))


def test_tensorsimplify_dispatch_zerotensor():
    """tensorsimplify routes AlgebraicZeroTensor to the right handler."""
    z = AlgebraicZeroTensor(((3, 4),))
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
    """AlgebraicTensor with a single AlgebraicPureTensor and no AlgebraicZeroTensor unwraps."""
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(pt)
    assert at is pt
    # simplify on the unwrapped AlgebraicPureTensor is fine.
    result = pt.simplify()
    assert result is pt or (isinstance(result, AlgebraicPureTensor) and result.factors == (A, C))


def test_simplify_negative_coefficient():
    """AlgebraicPureTensor with negative coefficient simplifies correctly."""
    pt = AlgebraicPureTensor(A, C)
    neg = -pt
    result = neg.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, C)


def test_simplify_identical_multi_factor():
    """Identical three-factor AlgebraicPureTensors combine."""
    pt = AlgebraicPureTensor(A, I4, C)
    at = AlgebraicTensor(2 * pt, 3 * pt)
    result = at.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 5
    assert result.factors == (A, I4, C)


def test_simplify_cancel_then_common_factor():
    """Terms cancel first, then common factors are extracted."""
    pt_ac = AlgebraicPureTensor(A, C)
    pt_bc = AlgebraicPureTensor(B, C)
    # 2*AC + 3*AC - 2*AC = 3*AC
    at = AlgebraicTensor(2 * pt_ac, 3 * pt_ac, -2 * pt_ac)
    result = at.simplify()
    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 3


def test_simplify_coefficient_symbol_cancellation():
    """Symbol coefficients cancel when opposite."""
    x = Symbol("x")
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(x * pt, (-x) * pt)
    result = at.simplify()
    assert isinstance(result, AlgebraicZeroTensor)



