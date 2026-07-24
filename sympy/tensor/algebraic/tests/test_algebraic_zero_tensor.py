"""Tests for algebraic_zero_tensor.py"""
from sympy.core.symbol import Symbol
from sympy.matrices.expressions import MatrixSymbol, ZeroMatrix
from sympy.testing.pytest import raises
from sympy.tensor.algebraic import AlgebraicZeroTensor


# ---------------------------------------------------------------------------
# Multi-factor AlgebraicZeroTensor tests
# ---------------------------------------------------------------------------

def test_algebraic_zero_tensor_constructor():
    """Test construction of AlgebraicZeroTensor with multi-factor shapes."""
    # Multiple factor shapes
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z2.shape == ((3, 4), (4, 5))

    # Multi-factor shape via direct constructor
    Z6 = AlgebraicZeroTensor(((2, 3), (3, 4)))
    assert Z6.shape == ((2, 3), (3, 4))


def test_algebraic_zero_tensor_hash_equality():
    """Test hashing and equality of AlgebraicZeroTensor."""
    Z1 = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor((3, 4))
    Z3 = AlgebraicZeroTensor((4, 3))

    assert Z1 == Z2
    assert Z1 != Z3
    assert hash(Z1) == hash(Z2)
    assert hash(Z1) != hash(Z3)

    # Multi-factor
    Z4 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z5 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z6 = AlgebraicZeroTensor(((3, 4), (5, 4)))

    assert Z4 == Z5
    assert Z4 != Z6
    assert hash(Z4) == hash(Z5)
    assert hash(Z4) != hash(Z6)


def test_algebraic_zero_tensor_negation():
    """Test negation of multi-factor AlgebraicZeroTensor."""
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert -Z2 is Z2


def test_algebraic_zero_tensor_addition():
    """Test addition with AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))

    from sympy.tensor.algebraic.algebraic_tensor import ShapeMismatchError

    # Addition with scalar should raise ShapeMismatchError
    raises(ShapeMismatchError, lambda: Z + 5)
    raises(ShapeMismatchError, lambda: 5 + Z)

    # Addition with same-shape zero tensor
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert (Z + Z2).shape == Z.shape

    # Shape mismatch should raise
    Z3 = AlgebraicZeroTensor((4, 5))
    raises(ShapeMismatchError, lambda: Z + Z3)

    # Multi-factor shape mismatch
    Z4 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z5 = AlgebraicZeroTensor(((3, 4), (5, 4)))
    raises(ShapeMismatchError, lambda: Z4 + Z5)


def test_algebraic_zero_tensor_subtraction():
    """Test subtraction with multi-factor AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))

    from sympy.tensor.algebraic.algebraic_tensor import ShapeMismatchError

    # Subtraction with scalar should raise ShapeMismatchError
    raises(ShapeMismatchError, lambda: Z - 5)
    raises(ShapeMismatchError, lambda: 5 - Z)

    # Z - Z = zero tensor of same shape
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z - Z2
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == Z.shape


def test_algebraic_zero_tensor_multiplication():
    """Test multiplication with multi-factor AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))

    # Scalar multiplication returns self
    assert 5 * Z is Z
    assert Z * 5 is Z
    assert 0 * Z is Z

    # Commutative symbol multiplication
    x = Symbol('x')
    assert x * Z is Z
    assert Z * x is Z


def test_algebraic_zero_tensor_properties():
    """Test properties of multi-factor AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))

    # is_zero should be True
    assert Z.is_zero is True

    # is_commutative should be True
    assert Z.is_commutative is True

    # is_AlgebraicZeroTensor should be True
    assert Z.is_AlgebraicZeroTensor is True

    # free_symbols should be empty
    assert Z.free_symbols == set()

    # args should be empty (atomic behavior)
    assert Z.args == ()

    # commutativity_pattern should be all 1s
    assert Z.commutativity_pattern == (1, 1)


def test_algebraic_zero_tensor_conjugate():
    """Test conjugation of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4), (4, 5))
    assert Z.conjugate() is Z
    assert Z._eval_conjugate() is Z


def test_algebraic_zero_tensor_transpose():
    """Test transpose of multi-factor AlgebraicZeroTensor."""
    Z2 = AlgebraicZeroTensor(((1, 2), (3, 4)))
    Z2T = Z2.T
    assert Z2T.shape == ((2, 1), (4, 3))

    # Double transpose returns self
    assert Z2.T.T.shape == Z2.shape


def test_algebraic_zero_tensor_diff():
    """Test differentiation of multi-factor AlgebraicZeroTensor."""
    from sympy.abc import x, y
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))

    # Derivative of zero is zero
    assert Z.diff(x) is Z
    assert Z.diff(x, x) is Z
    assert Z.diff(x, y) is Z


def test_algebraic_zero_tensor_doit():
    """Test doit of multi-factor AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.doit() is Z
    assert Z.doit(deep=False) is Z


def test_algebraic_zero_tensor_expand():
    """Test expand of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4,5)))
    assert Z.expand() is Z


def test_algebraic_zero_tensor_bool():
    """Test boolean conversion of multi-factor AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert bool(Z) is False


def test_algebraic_zero_tensor_pickle():
    """Test pickling of AlgebraicZeroTensor."""
    import pickle
    Z = AlgebraicZeroTensor((3, 4))
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))

    Z_pickled = pickle.loads(pickle.dumps(Z))
    assert Z_pickled == Z
    assert Z_pickled.shape == Z.shape

    Z2_pickled = pickle.loads(pickle.dumps(Z2))
    assert Z2_pickled == Z2
    assert Z2_pickled.shape == Z2.shape


# ---------------------------------------------------------------------------
# Single-factor dispatch to ZeroMatrix
# ---------------------------------------------------------------------------
# When AlgebraicZeroTensor is constructed with a single-factor shape,
# __new__ delegates to ZeroMatrix. These tests verify the dispatch and
# the boundary behaviour between ZeroMatrix and the algebraic tensor types.
# ---------------------------------------------------------------------------

def test_single_factor_dispatch_to_zero_matrix():
    """Single-factor AlgebraicZeroTensor dispatches to ZeroMatrix."""
    Z = AlgebraicZeroTensor((3, 4))
    assert isinstance(Z, ZeroMatrix)
    assert Z.shape == (3, 4)

    Z2 = AlgebraicZeroTensor([3, 4])
    assert isinstance(Z2, ZeroMatrix)
    assert Z2.shape == (3, 4)

    Z3 = AlgebraicZeroTensor([(3, 4)])
    assert isinstance(Z3, ZeroMatrix)
    assert Z3.shape == (3, 4)

    Z4 = AlgebraicZeroTensor((2, 3))
    assert isinstance(Z4, ZeroMatrix)
    assert Z4.shape == (2, 3)


def test_single_factor_dispatch_properties():
    """Dispatched ZeroMatrix does NOT carry AlgebraicZeroTensor attributes."""
    Z = AlgebraicZeroTensor((3, 4))
    assert isinstance(Z, ZeroMatrix)

    # ZeroMatrix does NOT carry AlgebraicZeroTensor-specific attributes
    assert not getattr(Z, 'is_AlgebraicZeroTensor', False)
    assert not hasattr(Z, 'commutativity_pattern')


def test_single_factor_dispatch_subtraction():
    """Subtraction with dispatched ZeroMatrix rejects scalars via TypeError."""
    Z = AlgebraicZeroTensor((3, 4))

    # ZeroMatrix.__sub__ routes through MatAdd, which rejects scalars
    raises(TypeError, lambda: Z - 5)
    raises(TypeError, lambda: 5 - Z)


def test_zeromatrix_cannot_add_to_multifactor_zero_tensor():
    """ZeroMatrix cannot be added to AlgebraicZeroTensor — TypeError from MatAdd."""
    zm = AlgebraicZeroTensor((3, 4))  # dispatches to ZeroMatrix
    azt_multi = AlgebraicZeroTensor(((3, 4), (4, 5)))

    # ZeroMatrix.__add__ routes through MatAdd, which rejects non-MatrixExpr
    raises(TypeError, lambda: zm + azt_multi)
    raises(TypeError, lambda: zm - azt_multi)

    # Reverse: AlgebraicZeroTensor.__radd__ validates shape, sees flat (3,4)
    # vs nested ((3,4),(4,5)) — ShapeMismatchError
    from sympy.tensor.algebraic.algebraic_tensor import ShapeMismatchError
    raises(ShapeMismatchError, lambda: azt_multi + zm)
    raises(ShapeMismatchError, lambda: azt_multi - zm)


def test_zeromatrix_compose_with_zero_tensor():
    """ZeroMatrix * AlgebraicZeroTensor routes through MatMul; reverse raises TypeError."""
    zm = AlgebraicZeroTensor((3, 4))  # dispatches to ZeroMatrix
    azt = AlgebraicZeroTensor(((3, 4), (4, 5)))

    # ZeroMatrix.__mul__ routes through MatMul (not compose_algebraic_tensors),
    # producing a MatMul result rather than raising an error.
    result = zm * azt
    assert result.shape == zm.shape  # MatMul simplifies zero * anything -> zero

    # Reverse: AlgebraicZeroTensor.__mul__ rejects ZeroMatrix with TypeError
    raises(TypeError, lambda: azt * zm)


# ---------------------------------------------------------------------------
# subs tests
# ---------------------------------------------------------------------------

def test_algebraic_zero_tensor_subs():
    """Subs on AlgebraicZeroTensor returns self (no symbols to replace)."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    x = Symbol('x')
    result = Z.subs(x, 5)
    assert result == Z


def test_algebraic_zero_tensor_subs_empty_free_symbols():
    """AlgebraicZeroTensor has no free_symbols, so any subs is a no-op."""
    Z = AlgebraicZeroTensor(((2, 3), (3, 4)))
    x, y = Symbol('x'), Symbol('y')
    result = Z.subs({x: 1, y: 2})
    assert result == Z
    assert Z.free_symbols == set()


def test_algebraic_zero_tensor_free_symbols_empty():
    """AlgebraicZeroTensor always has empty free_symbols."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z.free_symbols == set()


def test_algebraic_zero_tensor_free_symbols_single_factor():
    """Single-factor zero tensor (ZeroMatrix) also has empty free_symbols."""
    Z = AlgebraicZeroTensor((3, 4))  # dispatches to ZeroMatrix
    assert Z.free_symbols == set()


def test_algebraic_zero_tensor_simplify():
    """simplify on AlgebraicZeroTensor returns self."""
    from sympy import simplify
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert simplify(Z) is Z
    assert Z._eval_simplify() is Z


