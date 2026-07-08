from __future__ import annotations

from sympy.core import S, Symbol
from sympy.core.sympify import sympify
from sympy.matrices.expressions import MatrixSymbol
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    algebraic_zero_tensor,
)


# ---------------------------------------------------------------------------
# Constructor tests
# ---------------------------------------------------------------------------

def test_constructor_basic():
    """Basic construction of AlgebraicZeroTensor."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert isinstance(zt, AlgebraicZeroTensor)
    assert zt.shape == ((3, 4), (4, 5))


def test_constructor_bare_pair_wrapped():
    """Bare (m, n) pair is wrapped as ((m, n),)."""
    zt = AlgebraicZeroTensor((3, 4))
    assert zt.shape == ((3, 4),)


def test_constructor_list_wrapped():
    """List input is converted to tuple of tuples."""
    zt = AlgebraicZeroTensor([(3, 4)])
    assert zt.shape == ((3, 4),)


def test_constructor_tuple_unchanged():
    """Tuple of tuples is preserved."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert zt.shape == ((3, 4), (4, 5))


def test_constructor_empty_shape():
    """Empty shape is allowed."""
    zt = AlgebraicZeroTensor(())
    assert zt.shape == ()


def test_equality_same_shape():
    """Two AlgebraicZeroTensors with same shape are equal."""
    zt1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    zt2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert zt1 == zt2


def test_inequality_different_shape():
    """AlgebraicZeroTensors with different shapes are not equal."""
    zt1 = AlgebraicZeroTensor(((3, 4),))
    zt2 = AlgebraicZeroTensor(((2, 2),))
    assert zt1 != zt2


def test_hash_consistent():
    """Hash is consistent with equality."""
    zt1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    zt2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert hash(zt1) == hash(zt2)


def test_pickle_roundtrip():
    """__getnewargs__ allows pickle reconstruction."""
    import pickle
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    data = pickle.dumps(zt)
    zt2 = pickle.loads(data)
    assert zt2 == zt
    assert zt2.shape == zt.shape


# ---------------------------------------------------------------------------
# Properties tests
# ---------------------------------------------------------------------------

def test_shape_property():
    """shape returns the stored tensor shape."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert zt.shape == ((3, 4), (4, 5))


def test_tensor_shape_alias():
    """tensor_shape is an alias for shape."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert zt.tensor_shape == zt.shape


def test_commutativity_shape():
    """commutativity_shape is all-1s for zero tensor."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert zt.commutativity_shape == (1, 1)


def test_commutativity_shape_single():
    """commutativity_shape for single-factor zero tensor."""
    zt = AlgebraicZeroTensor(((3, 4),))
    assert zt.commutativity_shape == (1,)


def test_is_commutative_true():
    """AlgebraicZeroTensor is commutative."""
    zt = AlgebraicZeroTensor(((3, 4),))
    assert zt.is_commutative is True


def test_is_zero_none():
    """is_zero is None (not True) to prevent Mul collapse."""
    zt = AlgebraicZeroTensor(((3, 4),))
    assert zt.is_zero is None


def test_free_symbols_empty():
    """free_symbols returns empty set."""
    zt = AlgebraicZeroTensor(((3, 4),))
    assert zt.free_symbols == set()


def test_sympify_identity():
    """sympify returns the object unchanged."""
    zt = AlgebraicZeroTensor(((3, 4),))
    assert sympify(zt) is zt


# ---------------------------------------------------------------------------
# Arithmetic operator tests
# ---------------------------------------------------------------------------

def test_negation_returns_self():
    """Negation of zero tensor returns self."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert -zt is zt


def test_add_returns_other():
    """Addition returns the other operand (identity)."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(A, C)
    result = zt + pt
    assert result is pt


def test_radd_returns_other():
    """Right-addition returns the other operand (identity)."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(A, C)
    result = pt.__radd__(zt)
    assert result is pt


def test_sub_returns_algebraic_tensor():
    """Subtraction delegates to AlgebraicTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(A, C)
    result = zt - pt
    assert isinstance(result, AlgebraicTensor)


def test_rsub_returns_algebraic_tensor():
    """Right-subtraction delegates to AlgebraicTensor."""
    mats = _make_matrices()
    A, C = mats["A"], mats["C"]

    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(A, C)
    result = pt.__rsub__(zt)
    assert isinstance(result, AlgebraicTensor)


def test_mul_by_number_returns_self():
    """Multiplication by a Number returns self."""
    zt = AlgebraicZeroTensor(((3, 4),))
    result = zt * 5
    assert result is zt


def test_mul_by_symbol_returns_self():
    """Multiplication by a Symbol returns self."""
    zt = AlgebraicZeroTensor(((3, 4),))
    x = Symbol("x")
    result = zt * x
    assert result is zt


def test_rmul_by_number_returns_self():
    """Right-multiplication by a Number returns self."""
    zt = AlgebraicZeroTensor(((3, 4),))
    result = 5 * zt
    assert result is zt


def test_rmul_by_symbol_returns_self():
    """Right-multiplication by a Symbol returns self."""
    zt = AlgebraicZeroTensor(((3, 4),))
    x = Symbol("x")
    result = x * zt
    assert result is zt


def test_bool_false():
    """Boolean conversion returns False."""
    zt = AlgebraicZeroTensor(((3, 4),))
    assert bool(zt) is False


def test_if_statement_false():
    """AlgebraicZeroTensor evaluates as False in conditionals."""
    zt = AlgebraicZeroTensor(((3, 4),))
    if zt:
        raise AssertionError("Should not reach here")


# ---------------------------------------------------------------------------
# Expand tests
# ---------------------------------------------------------------------------

def test_expand_returns_self():
    """Expand returns self (zero tensor is already expanded)."""
    zt = AlgebraicZeroTensor(((3, 4),))
    result = zt.expand()
    assert result is zt


# ---------------------------------------------------------------------------
# Copy tests
# ---------------------------------------------------------------------------

def test_copy_returns_self():
    """copy returns self (immutable atom)."""
    zt = AlgebraicZeroTensor(((3, 4),))
    assert zt.copy() is zt


# ---------------------------------------------------------------------------
# String representation tests
# ---------------------------------------------------------------------------

def test_str():
    """String representation shows shape."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    s = str(zt)
    assert "3" in s and "4" in s


def test_repr():
    """Repr contains class name."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    r = repr(zt)
    assert "AlgebraicZeroTensor" in r


# ---------------------------------------------------------------------------
# Convenience function tests
# ---------------------------------------------------------------------------

def test_algebraic_zero_tensor():
    """algebraic_zero_tensor is a convenience constructor."""
    zt = algebraic_zero_tensor(((3, 4), (4, 5)))
    assert isinstance(zt, AlgebraicZeroTensor)
    assert zt.shape == ((3, 4), (4, 5))


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _make_matrices():
    """Create a standard set of matrix symbols for testing."""
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 3, 4)
    C = MatrixSymbol("C", 4, 5)
    D = MatrixSymbol("D", 4, 5)
    return {"A": A, "B": B, "C": C, "D": D}
