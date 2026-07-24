"""Tests for algebraic_pure_tensor.py"""
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.core.numbers import I
from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices import ImmutableDenseMatrix
from sympy.testing.pytest import raises
from sympy.tensor.algebraic import AlgebraicPureTensor, algebraic_tensor_product
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor


def test_algebraic_pure_tensor_constructor():
    """Test construction of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    # Two factors
    T = AlgebraicPureTensor(A, B)
    assert T.shape == ((3, 4), (4, 5))
    assert T.factors == (A, B)
    assert T.coeff == S.One
    
    # With numeric coefficient
    T2 = AlgebraicPureTensor(2, A, B)
    assert T2.coeff == 2
    assert T2.factors == (A, B)
    
    # With symbolic coefficient
    x = Symbol('x')
    T3 = AlgebraicPureTensor(x, A, B)
    assert T3.coeff == x
    assert T3.factors == (A, B)
    
    # Single factor with coeff 1 unwraps
    T4 = AlgebraicPureTensor(A)
    assert T4 == A
    
    # Zero coefficient produces zero tensor
    Z = AlgebraicPureTensor(0, A, B)
    assert isinstance(Z, AlgebraicZeroTensor)
    assert Z.shape == ((3, 4), (4, 5))
    
    # Empty args should raise
    raises(ValueError, lambda: AlgebraicPureTensor())


def test_algebraic_pure_tensor_props():
    """Test properties of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(2, A, B)
    
    assert T.num_factors == 2
    assert T.shape == ((3, 4), (4, 5))
    assert T.coeff == 2
    assert T.factors == (A, B)
    
    # is_AlgebraicPureTensor
    assert T.is_AlgebraicPureTensor is True


def test_algebraic_pure_tensor_commutativity_pattern():
    """Test commutativity_pattern of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(A, B)
    assert T.commutativity_pattern == (0, 0)
    
    # Numeric matrix is commutative
    M = ImmutableDenseMatrix([[1, 2, 3, 4],
                              [5, 6, 7, 8],
                              [9, 10, 11, 12]])
    T2 = AlgebraicPureTensor(M, B)
    assert T2.commutativity_pattern == (1, 0)


def test_algebraic_pure_tensor_negation():
    """Test negation of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(A, B)
    T_neg = -T
    assert T_neg.coeff == S.NegativeOne
    assert T_neg.factors == (A, B)
    
    T2 = AlgebraicPureTensor(2, A, B)
    T2_neg = -T2
    assert T2_neg.coeff == -2
    assert T2_neg.factors == (A, B)


def test_algebraic_pure_tensor_multiplication():
    """Test multiplication of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 4, 2)
    D = MatrixSymbol('D', 5, 3)
    
    T = AlgebraicPureTensor(A, B)
    
    # Scalar multiplication
    T_scaled = T * 3
    assert T_scaled.coeff == 3
    assert T_scaled.factors == (A, B)
    
    # Zero multiplication
    T_zero = T * 0
    assert isinstance(T_zero, AlgebraicZeroTensor)
    assert T_zero.shape == ((3, 4), (4, 5))
    
    # One multiplication
    assert T * 1 is T
    
    # Composition with another pure tensor
    T2 = AlgebraicPureTensor(C, D)
    result = T * T2
    assert isinstance(result, AlgebraicPureTensor)
    
    # Left scalar multiplication
    T_scaled2 = 3 * T
    assert T_scaled2.coeff == 3
    assert T_scaled2.factors == (A, B)
    
    # Left zero multiplication
    T_zero2 = 0 * T
    assert isinstance(T_zero2, AlgebraicZeroTensor)


def test_algebraic_pure_tensor_addition():
    """Test addition of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    # Addition produces AlgebraicTensor
    result = T1 + T2
    assert isinstance(result, AlgebraicTensor)
    
    # Addition with zero tensor of same shape
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result2 = T1 + Z
    assert result2 == T1
    
    # Right addition
    result3 = T2 + T1
    assert isinstance(result3, AlgebraicTensor)


def test_algebraic_pure_tensor_subtraction():
    """Test subtraction of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    result = T1 - T2
    assert isinstance(result, AlgebraicTensor)
    
    result2 = T2 - T1
    assert isinstance(result2, AlgebraicTensor)


def test_algebraic_pure_tensor_transpose():
    """Test transpose of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(A, B)
    TT = T.T
    
    assert TT.shape == ((4, 3), (5, 4))
    assert TT.coeff == S.One
    
    # With coefficient
    T2 = AlgebraicPureTensor(2, A, B)
    T2T = T2.T
    assert T2T.coeff == 2
    assert T2T.shape == ((4, 3), (5, 4))


def test_algebraic_pure_tensor_conjugate():
    """Test conjugate of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(1 + I, A, B)
    TC = T.conjugate()
    
    assert TC.coeff == 1 - I
    
    # Numeric matrix conjugate
    M = ImmutableDenseMatrix([[1 + I, 2], [3, 4 - I]])
    N = MatrixSymbol('N', 2, 3)
    T2 = AlgebraicPureTensor(M, N)
    T2C = T2.conjugate()
    assert T2C.factors[0][0, 0] == 1 - I


def test_algebraic_pure_tensor_diff():
    """Test differentiation of AlgebraicPureTensor."""
    from sympy.abc import x
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    # Differentiate coefficient
    T = AlgebraicPureTensor(x**2, A, B)
    T_diff = T.diff(x)
    assert T_diff.coeff == 2*x
    
    # Symbol not present
    T2 = AlgebraicPureTensor(A, B)
    T2_diff = T2.diff(x)
    assert isinstance(T2_diff, AlgebraicZeroTensor)


def test_algebraic_pure_tensor_doit():
    """Test doit of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(A, B)
    assert T.doit() == T
    
    # With symbolic coefficient
    from sympy.abc import x
    T2 = AlgebraicPureTensor(x, A, B)
    T2_doit = T2.doit()
    assert T2_doit.coeff == x


def test_algebraic_pure_tensor_expand():
    """Test expand of AlgebraicPureTensor."""
    from sympy.matrices.expressions import MatAdd
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 4, 5)
    
    T = AlgebraicPureTensor(A, MatAdd(B, C))
    T_exp = T.expand()
    
    assert isinstance(T_exp, AlgebraicTensor)


def test_algebraic_pure_tensor_has_zero_term():
    """Test has_zero_term of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(A, B)
    assert T.has_zero_term() is False


def test_algebraic_pure_tensor_simplify():
    """Test simplify of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(2, A, B)
    T_simp = T.simplify()
    # Should return a simplified version
    assert isinstance(T_simp, AlgebraicPureTensor)


def test_algebraic_tensor_product():
    """Test algebraic_tensor_product function."""
    A = MatrixSymbol('A', 2, 3)
    v = MatrixSymbol('v', 3, 1)
    
    result = algebraic_tensor_product(A, v)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((2, 3), (3, 1))
    
    # With scalar
    result2 = algebraic_tensor_product(2, A, v)
    assert result2.coeff == 2


def test_algebraic_pure_tensor_pickle():
    """Test pickling of AlgebraicPureTensor."""
    import pickle
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(2, A, B)
    T_pickled = pickle.loads(pickle.dumps(T))
    
    assert T_pickled == T
    assert T_pickled.shape == T.shape
    assert T_pickled.coeff == T.coeff


def test_algebraic_pure_tensor_1x1_matrix():
    """Test AlgebraicPureTensor with 1x1 matrix."""
    M = ImmutableDenseMatrix([[5]])
    A = MatrixSymbol('A', 3, 4)

    T = AlgebraicPureTensor(M, A)
    assert T.shape == ((1, 1), (3, 4))
    assert T.commutativity_pattern == (1, 0)


def test_pure_tensor_add_with_dispatched_zeromatrix():
    """PureTensor cannot be added to a ZeroMatrix (dispatched single-factor zero)."""
    from sympy.tensor.algebraic.algebraic_tensor import ShapeMismatchError

    zm = AlgebraicZeroTensor((3, 4))  # dispatches to ZeroMatrix
    A = MatrixSymbol("A", 3, 4)
    B = MatrixSymbol("B", 4, 5)
    pt = AlgebraicPureTensor(A, B)

    # ZeroMatrix.__add__ routes through MatAdd, rejects non-MatrixExpr
    raises(TypeError, lambda: zm + pt)
    raises(TypeError, lambda: zm - pt)

    # Reverse: AlgebraicPureTensor.__radd__ validates shape mismatch
    raises(ShapeMismatchError, lambda: pt + zm)
    raises(ShapeMismatchError, lambda: pt - zm)


def test_pure_tensor_compose_with_dispatched_zeromatrix():
    """PureTensor * ZeroMatrix raises TypeError; ZeroMatrix * PureTensor through MatMul."""
    zm = AlgebraicZeroTensor((3, 4))  # dispatches to ZeroMatrix
    A = MatrixSymbol("A", 3, 2)
    B = MatrixSymbol("B", 4, 5)
    pt = AlgebraicPureTensor(A, B)

    # ZeroMatrix.__mul__ routes through MatMul, not compose_algebraic_tensors.
    # MatMul.as_coeff_matrices does not handle AlgebraicPureTensor (is_Matrix=False,
    # is_commutative=False) and raises NotImplementedError. This is a known Matrix
    # bug: __mul__ lacks sufficient instance checks for tensor types.
    raises(NotImplementedError, lambda: zm * pt)

    # Reverse: AlgebraicPureTensor.__mul__ rejects ZeroMatrix with TypeError
    raises(TypeError, lambda: pt * zm)


# ---------------------------------------------------------------------------
# subs tests
# ---------------------------------------------------------------------------

def test_pure_tensor_subs_coefficient():
    """Substitute a symbol in the coefficient."""
    x = Symbol('x')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T = AlgebraicPureTensor(x**2, A, B)
    result = T.subs(x, 3)
    assert result.coeff == 9
    assert result.factors == (A, B)


def test_pure_tensor_subs_matrix_factor():
    """Substitute a MatrixSymbol factor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    T = AlgebraicPureTensor(A, B)
    result = T.subs(A, C)
    assert result.factors[0] == C
    assert result.factors[1] == B


def test_pure_tensor_subs_concrete_matrix_entries():
    """Substitute a symbol inside concrete matrix entries."""
    x = Symbol('x')
    M = ImmutableDenseMatrix([[x, 1], [2, x]])
    N = MatrixSymbol('N', 2, 3)
    T = AlgebraicPureTensor(M, N)
    result = T.subs(x, 5)
    assert result.factors[0][0, 0] == 5
    assert result.factors[0][1, 1] == 5


def test_pure_tensor_subs_multiple():
    """Substitute multiple symbols at once."""
    x, y = Symbol('x'), Symbol('y')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    T = AlgebraicPureTensor(x, A, B)
    result = T.subs({x: y, A: C})
    assert result.coeff == y
    assert result.factors[0] == C


def test_pure_tensor_subs_no_match():
    """Substitution that doesn't match anything returns unchanged tensor."""
    x, y = Symbol('x'), Symbol('y')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T = AlgebraicPureTensor(x, A, B)
    result = T.subs(y, 10)
    assert result == T


def test_pure_tensor_subs_to_zero():
    """Substituting coefficient to zero yields AlgebraicZeroTensor."""
    x = Symbol('x')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T = AlgebraicPureTensor(x, A, B)
    result = T.subs(x, 0)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


# ---------------------------------------------------------------------------
# free_symbols tests
# ---------------------------------------------------------------------------

def test_pure_tensor_free_symbols_coefficient():
    """free_symbols collects symbols from the coefficient."""
    x, y = Symbol('x'), Symbol('y')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T = AlgebraicPureTensor(x**2, A, B)
    assert x in T.free_symbols
    assert y not in T.free_symbols


def test_pure_tensor_free_symbols_no_coefficient():
    """Tensor with no symbolic coefficient has only factor symbols."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T = AlgebraicPureTensor(A, B)
    assert T.free_symbols == {A, B}


def test_pure_tensor_free_symbols_concrete_matrix():
    """free_symbols collects symbols from concrete matrix entries."""
    x, y = Symbol('x'), Symbol('y')
    M = ImmutableDenseMatrix([[x, 1], [2, y]])
    N = MatrixSymbol('N', 2, 3)
    T = AlgebraicPureTensor(M, N)
    assert x in T.free_symbols
    assert y in T.free_symbols
    assert N in T.free_symbols


def test_pure_tensor_free_symbols_mixed():
    """free_symbols unions symbols from coefficient and all factors."""
    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    A = MatrixSymbol('A', 3, 4)
    M = ImmutableDenseMatrix([[x, 1], [2, x]])
    N = MatrixSymbol('N', 2, 3)
    T = AlgebraicPureTensor(y, M, N)
    assert T.free_symbols == {y, x, N}
    assert z not in T.free_symbols


def test_pure_tensor_free_symbols_numeric_only():
    """Tensor with purely numeric content has no free symbols."""
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    N = ImmutableDenseMatrix([[5, 6, 7], [8, 9, 10]])
    T = AlgebraicPureTensor(3, M, N)
    assert T.free_symbols == set()


def test_pure_tensor_free_symbols_three_factors():
    """free_symbols works for tensors with three or more factors."""
    x = Symbol('x')
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    T = AlgebraicPureTensor(x, A, B, C)
    assert T.free_symbols == {x, A, B, C}


# ---------------------------------------------------------------------------
# hash tests
# ---------------------------------------------------------------------------

def test_pure_tensor_hash_consistency():
    """Two equal AlgebraicPureTensor objects have the same hash."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(A, B)
    assert hash(T1) == hash(T2)
    assert T1 == T2


def test_pure_tensor_hash_with_coefficient():
    """Hash is consistent when coefficient is present."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(2, A, B)
    assert hash(T1) == hash(T2)
    assert T1 == T2


def test_pure_tensor_hash_different_factors():
    """Tensors with different factors have different content."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    assert T1 != T2


def test_pure_tensor_hash_noncommutative_order():
    """Order matters: A⊗B != B⊗A (different shapes)."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(B, A)
    assert T1 != T2


# ---------------------------------------------------------------------------
# 3+ factor tensor tests
# ---------------------------------------------------------------------------

def test_three_factor_pure_tensor():
    """Test construction and properties of a 3-factor pure tensor."""
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    T = AlgebraicPureTensor(A, B, C)
    assert T.num_factors == 3
    assert T.shape == ((2, 3), (3, 4), (4, 5))
    assert T.factors == (A, B, C)
    assert T.coeff == S.One


def test_three_factor_with_coefficient():
    """Test 3-factor pure tensor with symbolic coefficient."""
    from sympy.abc import x
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    T = AlgebraicPureTensor(x, A, B, C)
    assert T.num_factors == 3
    assert T.coeff == x
    assert T.factors == (A, B, C)


def test_compose_three_factor_tensors():
    """Test composition of two 3-factor pure tensors."""
    from sympy.tensor.algebraic import compose_algebraic_pure_tensors
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    D = MatrixSymbol('D', 3, 2)
    E = MatrixSymbol('E', 4, 3)
    F = MatrixSymbol('F', 5, 4)
    T1 = AlgebraicPureTensor(A, B, C)
    T2 = AlgebraicPureTensor(D, E, F)
    result = compose_algebraic_pure_tensors(T1, T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 3
    assert result.shape == ((2, 2), (3, 3), (4, 4))


def test_compose_three_factor_via_mul():
    """Test composition of 3-factor tensors via __mul__."""
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    D = MatrixSymbol('D', 3, 2)
    E = MatrixSymbol('E', 4, 3)
    F = MatrixSymbol('F', 5, 4)
    T1 = AlgebraicPureTensor(A, B, C)
    T2 = AlgebraicPureTensor(D, E, F)
    result = T1 * T2
    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 3


def test_expand_three_factors():
    """Test expand with MatAdd in a 3-factor tensor."""
    from sympy.matrices.expressions import MatAdd
    A = MatrixSymbol('A', 2, 3)
    B1 = MatrixSymbol('B1', 3, 4)
    B2 = MatrixSymbol('B2', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    T = AlgebraicPureTensor(A, MatAdd(B1, B2), C)
    expanded = T.expand()
    assert isinstance(expanded, AlgebraicTensor)
    assert len(expanded.args) == 2


def test_transpose_three_factors():
    """Test transpose of a 3-factor pure tensor."""
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    T = AlgebraicPureTensor(A, B, C)
    TT = T.T
    assert TT.shape == ((3, 2), (4, 3), (5, 4))
    assert TT.num_factors == 3


def test_conjugate_three_factors():
    """Test conjugate of a 3-factor pure tensor with complex coefficient."""
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    T = AlgebraicPureTensor(1 + I, A, B, C)
    TC = T.conjugate()
    assert TC.coeff == 1 - I
    assert TC.num_factors == 3


def test_diff_three_factors():
    """Test differentiation of a 3-factor pure tensor."""
    from sympy.abc import x
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    T = AlgebraicPureTensor(x**2, A, B, C)
    T_diff = T.diff(x)
    assert T_diff.coeff == 2 * x
    assert T_diff.num_factors == 3


def test_subs_three_factors():
    """Test substitution in a 3-factor pure tensor."""
    x = Symbol('x')
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    D = MatrixSymbol('D', 2, 3)
    T = AlgebraicPureTensor(x, A, B, C)
    result = T.subs(x, 5)
    assert result.coeff == 5
    result2 = T.subs(A, D)
    assert result2.factors[0] == D
    assert result2.factors[1] == B
    assert result2.factors[2] == C


def test_four_factor_pure_tensor():
    """Test construction of a 4-factor pure tensor."""
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    D = MatrixSymbol('D', 5, 6)
    T = AlgebraicPureTensor(A, B, C, D)
    assert T.num_factors == 4
    assert T.shape == ((2, 3), (3, 4), (4, 5), (5, 6))


def test_compose_four_factor_tensors():
    """Test composition of two 4-factor pure tensors."""
    from sympy.tensor.algebraic import compose_algebraic_pure_tensors
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 3, 4)
    C = MatrixSymbol('C', 4, 5)
    D = MatrixSymbol('D', 5, 6)
    E = MatrixSymbol('E', 3, 2)
    F = MatrixSymbol('F', 4, 3)
    G = MatrixSymbol('G', 5, 4)
    H = MatrixSymbol('H', 6, 5)
    T1 = AlgebraicPureTensor(A, B, C, D)
    T2 = AlgebraicPureTensor(E, F, G, H)
    result = compose_algebraic_pure_tensors(T1, T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 4
    assert result.shape == ((2, 2), (3, 3), (4, 4), (5, 5))
