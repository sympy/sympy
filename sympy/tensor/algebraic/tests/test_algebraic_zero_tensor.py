"""Tests for algebraic_zero_tensor.py"""
from sympy.core.symbol import Symbol
from sympy.matrices.expressions import MatrixSymbol
from sympy.testing.pytest import raises
from sympy.tensor.algebraic import AlgebraicZeroTensor, algebraic_zero_tensor


def test_algebraic_zero_tensor_constructor():
    """Test construction of AlgebraicZeroTensor."""
    # Single factor shape
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.shape == ((3, 4),)
    
    # Multiple factor shapes
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z2.shape == ((3, 4), (4, 5))
    
    # List input
    Z3 = AlgebraicZeroTensor([(3, 4)])
    assert Z3.shape == ((3, 4),)
    
    # Bare (m, n) or [m, n] -> wrapped as ((m, n),)
    Z4 = AlgebraicZeroTensor([3, 4])
    assert Z4.shape == ((3, 4),)
    
    # Convenience function
    Z5 = algebraic_zero_tensor((2, 3))
    assert Z5.shape == ((2, 3),)
    
    Z6 = algebraic_zero_tensor(((2, 3), (3, 4)))
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
    """Test negation of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    assert -Z is Z
    assert (-Z).shape == Z.shape
    
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert -Z2 is Z2


def test_algebraic_zero_tensor_addition():
    """Test addition with AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    
    # Addition with scalar (should return the other operand)
    assert Z + 5 == 5
    assert 5 + Z == 5
    
    # Addition with same-shape zero tensor
    Z2 = AlgebraicZeroTensor((3, 4))
    assert (Z + Z2).shape == Z.shape
    
    # Shape mismatch should raise
    Z3 = AlgebraicZeroTensor((4, 5))
    from sympy.tensor.algebraic.algebraic_tensor import ShapeMismatchError
    raises(ShapeMismatchError, lambda: Z + Z3)
    
    # Multi-factor shape mismatch
    Z4 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z5 = AlgebraicZeroTensor(((3, 4), (5, 4)))
    raises(ShapeMismatchError, lambda: Z4 + Z5)


def test_algebraic_zero_tensor_subtraction():
    """Test subtraction with AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    
    # Z - x = -x
    assert (Z - 5) == -5
    assert (5 - Z) == 5
    
    # Z - Z = zero tensor of same shape
    Z2 = AlgebraicZeroTensor((3, 4))
    result = Z - Z2
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == Z.shape


def test_algebraic_zero_tensor_multiplication():
    """Test multiplication with AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    
    # Scalar multiplication returns self
    assert 5 * Z is Z
    assert Z * 5 is Z
    assert 0 * Z is Z
    
    # Commutative symbol multiplication
    x = Symbol('x')
    assert x * Z is Z
    assert Z * x is Z
    
    # Non-commutative should delegate to composition
    A = MatrixSymbol('A', 3, 4)
    result = Z * A
    # Should return a composed zero tensor
    assert isinstance(result, AlgebraicZeroTensor)


def test_algebraic_zero_tensor_properties():
    """Test properties of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    
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
    assert Z.commutativity_pattern == (1,)
    
    # Multi-factor commutativity_pattern
    Z2 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert Z2.commutativity_pattern == (1, 1)


def test_algebraic_zero_tensor_conjugate():
    """Test conjugation of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.conjugate() is Z
    assert Z._eval_conjugate() is Z


def test_algebraic_zero_tensor_transpose():
    """Test transpose of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    ZT = Z.T
    assert ZT.shape == ((4, 3),)
    
    # Multi-factor transpose
    Z2 = AlgebraicZeroTensor(((1, 2), (3, 4)))
    Z2T = Z2.T
    assert Z2T.shape == ((2, 1), (4, 3))
    
    # Double transpose returns self
    assert Z.T.T.shape == Z.shape


def test_algebraic_zero_tensor_diff():
    """Test differentiation of AlgebraicZeroTensor."""
    from sympy.abc import x, y
    Z = AlgebraicZeroTensor((3, 4))
    
    # Derivative of zero is zero
    assert Z.diff(x) is Z
    assert Z.diff(x, x) is Z
    assert Z.diff(x, y) is Z


def test_algebraic_zero_tensor_doit():
    """Test doit of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.doit() is Z
    assert Z.doit(deep=False) is Z


def test_algebraic_zero_tensor_expand():
    """Test expand of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
    assert Z.expand() is Z


def test_algebraic_zero_tensor_bool():
    """Test boolean conversion of AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor((3, 4))
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
