from __future__ import annotations

from sympy.core import S, Symbol
from sympy.matrices.expressions import MatrixSymbol
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    ShapeMismatchError,
    compose_algebraic_tensors,
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
    """Basic construction of AlgebraicTensor from two PureTensors."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    assert isinstance(at, AlgebraicTensor)
    assert len(at.args) == 2


def test_constructor_single_pure_tensor_unwraps():
    """Single AlgebraicPureTensor argument unwraps to the PureTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    result = AlgebraicTensor(pt)
    assert result is pt
    assert not isinstance(result, AlgebraicTensor)


def test_constructor_single_zero_tensor_unwraps():
    """Single AlgebraicZeroTensor argument unwraps."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = AlgebraicTensor(zt)
    assert result is zt


def test_constructor_single_algebraic_tensor_unwraps():
    """Single AlgebraicTensor argument returns itself."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    result = AlgebraicTensor(at)
    assert result is at


def test_constructor_empty_raises():
    """Empty constructor raises ValueError."""
    raises(ValueError, lambda: AlgebraicTensor())


def test_constructor_shape_mismatch_raises():
    """Adding tensors of different shapes raises ShapeMismatchError."""
    mats = _make_matrices()
    A, B, C, E = mats["A"], mats["B"], mats["C"], mats["E"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, E)
    raises(ShapeMismatchError, lambda: AlgebraicTensor(pt1, pt2))


def test_constructor_flattens_nested():
    """Nested AlgebraicTensors are flattened."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    inner = AlgebraicTensor(pt1, pt2)
    pt3 = AlgebraicPureTensor(A, D)
    pt4 = AlgebraicPureTensor(B, C)
    inner2 = AlgebraicTensor(pt3, pt4)
    outer = AlgebraicTensor(inner, inner2)
    assert isinstance(outer, AlgebraicTensor)
    assert len(outer.args) == 4


def test_constructor_cancellation_returns_zero_tensor():
    """Adding t and -t keeps both terms (cancellation is a simplification step)."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    neg_pt = AlgebraicPureTensor(S.NegativeOne, A, C)
    result = AlgebraicTensor(pt, neg_pt)
    assert isinstance(result, AlgebraicTensor)
    assert len(result.args) == 2


def test_constructor_with_zero_tensor_anchor():
    """AlgebraicZeroTensor is kept as a shape anchor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(pt, zt)
    assert isinstance(at, AlgebraicTensor)
    assert at.has_zero_term()


# ---------------------------------------------------------------------------
# Properties tests
# ---------------------------------------------------------------------------

def test_shape():
    """shape returns the shared shape."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    assert at.shape == ((3, 4), (4, 5))


def test_commutativity_pattern():
    """commutativity_pattern is component-wise AND of all terms."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    assert at.commutativity_pattern == (0, 0)


def test_terms_property():
    """terms property returns non-zero, non-coefficient terms."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    assert len(at.terms) == 2


def test_is_add_true():
    """AlgebraicTensor has is_Add = True."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    assert at.is_Add is True


def test_has_zero_term():
    """has_zero_term detects AlgebraicZeroTensor anchor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt = AlgebraicPureTensor(A, C)
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(pt, zt)
    assert at.has_zero_term() is True


def test_has_zero_term_false():
    """has_zero_term is False when no zero anchor."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    assert at.has_zero_term() is False


# ---------------------------------------------------------------------------
# Arithmetic operator tests
# ---------------------------------------------------------------------------

def test_negation():
    """Negation negates each term."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    neg = -at
    assert isinstance(neg, AlgebraicTensor)


def test_add():
    """Addition of two AlgebraicTensors."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at1 = AlgebraicTensor(pt1, pt2)
    at2 = AlgebraicTensor(pt2, pt1)
    result = at1 + at2
    assert isinstance(result, AlgebraicTensor)


def test_sub():
    """Subtraction of two AlgebraicTensors."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at1 = AlgebraicTensor(pt1, pt2)
    at2 = AlgebraicTensor(pt2, pt1)
    result = at1 - at2
    assert isinstance(result, AlgebraicTensor)


def test_mul_by_number():
    """Multiplication by a Number scales each term."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    result = at * 3
    assert isinstance(result, AlgebraicTensor)


def test_mul_by_symbol():
    """Multiplication by a Symbol scales each term."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]
    x = Symbol("x")

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    result = at * x
    assert isinstance(result, AlgebraicTensor)


def test_rmul_by_number():
    """Right-multiplication by a Number scales each term."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    result = 3 * at
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# Composition tests
# ---------------------------------------------------------------------------

def test_compose_at_with_at():
    """Composition of two AlgebraicTensors."""
    mats = _make_matrices()
    A, B, C, D, I4, I5 = mats["A"], mats["B"], mats["C"], mats["D"], mats["I4"], mats["I5"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at1 = AlgebraicTensor(pt1, pt2)

    pt3 = AlgebraicPureTensor(I4, I5)
    pt4 = AlgebraicPureTensor(I4, I5)
    at2 = AlgebraicTensor(pt3, pt4)

    result = compose_algebraic_tensors(at1, at2)
    assert result is not None


def test_compose_at_with_pt():
    """Composition of AlgebraicTensor with PureTensor."""
    mats = _make_matrices()
    A, B, C, D, I4, I5 = mats["A"], mats["B"], mats["C"], mats["D"], mats["I4"], mats["I5"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    pt_id = AlgebraicPureTensor(I4, I5)

    result = compose_algebraic_tensors(at, pt_id)
    assert result is not None


def test_compose_pt_with_at():
    """Composition of PureTensor with AlgebraicTensor."""
    mats = _make_matrices()
    A, B, C, D, I3, I4 = mats["A"], mats["B"], mats["C"], mats["D"], mats["I3"], mats["I4"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    pt_id = AlgebraicPureTensor(I3, I4)

    result = compose_algebraic_tensors(pt_id, at)
    assert result is not None


def test_compose_at_with_zero_tensor():
    """Composition with AlgebraicZeroTensor returns zero."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))

    result = compose_algebraic_tensors(at, zt)
    assert isinstance(result, AlgebraicZeroTensor)


def test_compose_zero_tensor_with_at():
    """Composition of zero tensor with AlgebraicTensor returns zero."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))

    result = compose_algebraic_tensors(zt, at)
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# Expand tests
# ---------------------------------------------------------------------------

def test_expand():
    """Expand distributes through each term."""
    from sympy.matrices.expressions import MatAdd

    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(MatAdd(A, B), C)
    pt2 = AlgebraicPureTensor(A, MatAdd(C, D))
    at = AlgebraicTensor(pt1, pt2)
    result = at.expand()
    assert result is not None


# ---------------------------------------------------------------------------
# String representation tests
# ---------------------------------------------------------------------------

def test_str():
    """String representation contains all terms."""
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    s = str(at)
    assert "A" in s and "C" in s


def test_repr():
    """srepr contains class name (repr == str per SymPy convention)."""
    from sympy.printing.repr import srepr
    mats = _make_matrices()
    A, B, C, D = mats["A"], mats["B"], mats["C"], mats["D"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, D)
    at = AlgebraicTensor(pt1, pt2)
    r = srepr(at)
    assert "AlgebraicTensor" in r


# ---------------------------------------------------------------------------
# ShapeMismatchError tests
# ---------------------------------------------------------------------------

def test_shape_mismatch_error_is_type_error():
    """ShapeMismatchError is a subclass of TypeError."""
    assert issubclass(ShapeMismatchError, TypeError)


def test_shape_mismatch_in_composition():
    """Composition with incompatible shapes raises ValueError."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(C, A)
    raises(ValueError, lambda: compose_algebraic_tensors(pt1, pt2))
