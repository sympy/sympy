"""Tests for algebraic_tensor.py"""
from sympy.core.numbers import I
from sympy.core.symbol import Symbol
from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices import ImmutableDenseMatrix
from sympy.testing.pytest import raises
from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor
from sympy.tensor.algebraic.algebraic_tensor import ShapeMismatchError, compose_algebraic_tensors


def test_algebraic_tensor_constructor():
    """Test construction of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    # Two terms
    S = AlgebraicTensor(T1, T2)
    assert S.shape == ((3, 4), (4, 5))
    assert len(S.args) == 2
    
    # Single term unwraps
    S2 = AlgebraicTensor(T1)
    assert S2 == T1
    
    # From addition
    S3 = T1 + T2
    assert isinstance(S3, AlgebraicTensor)
    
    # Cancellation produces zero tensor
    S4 = AlgebraicTensor(T1, -T1)
    assert isinstance(S4, AlgebraicZeroTensor)
    
    # Empty args should raise
    raises(ValueError, lambda: AlgebraicTensor())


def test_algebraic_tensor_shape():
    """Test shape property of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    assert S.shape == ((3, 4), (4, 5))


def test_algebraic_tensor_commutativity_pattern():
    """Test commutativity_pattern of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    assert S.commutativity_pattern == (0, 0)
    
    # With numeric matrix
    M = ImmutableDenseMatrix([[1, 2, 3, 4],
                              [5, 6, 7, 8],
                              [9, 10, 11, 12]])
    T3 = AlgebraicPureTensor(M, B)
    S2 = AlgebraicTensor(T1, T3)
    assert S2.commutativity_pattern == (0, 0)


def test_algebraic_tensor_terms():
    """Test terms property of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    assert len(S.terms) == 2


def test_algebraic_tensor_has_zero_term():
    """Test has_zero_term of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    
    S2 = AlgebraicTensor(T1)
    assert S2.has_zero_term() is False


def test_algebraic_tensor_addition():
    """Test addition of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    E = MatrixSymbol('E', 3, 4)
    F = MatrixSymbol('F', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(E, F)
    
    S1 = AlgebraicTensor(T1, T2)
    
    # Add pure tensor
    S2 = S1 + T3
    assert isinstance(S2, AlgebraicTensor)
    assert len(S2.args) == 3
    
    # Add another AlgebraicTensor
    S3 = AlgebraicTensor(T3)
    S4 = S1 + S3
    assert isinstance(S4, AlgebraicTensor)
    
    # Right addition
    S5 = T3 + S1
    assert isinstance(S5, AlgebraicTensor)
    
    # Shape mismatch should raise
    G = MatrixSymbol('G', 2, 3)
    H = MatrixSymbol('H', 3, 4)
    T4 = AlgebraicPureTensor(G, H)
    raises(ShapeMismatchError, lambda: S1 + T4)


def test_algebraic_tensor_subtraction():
    """Test subtraction of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    S = AlgebraicTensor(T1, T2)
    
    result = S - T1
    # Result may be unwrapped to AlgebraicPureTensor if single term
    assert isinstance(result, (AlgebraicTensor, AlgebraicPureTensor))
    
    # Subtract same tensor
    result2 = S - S
    assert isinstance(result2, AlgebraicZeroTensor)
    
    # Right subtraction
    result3 = T1 - S
    assert isinstance(result3, AlgebraicTensor)


def test_algebraic_tensor_multiplication():
    """Test multiplication of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    # Scalar multiplication
    S_scaled = S * 2
    assert isinstance(S_scaled, AlgebraicTensor)
    
    # Left scalar multiplication
    S_scaled2 = 2 * S
    assert isinstance(S_scaled2, AlgebraicTensor)
    
    # Composition with pure tensor
    E = MatrixSymbol('E', 4, 2)
    F = MatrixSymbol('F', 5, 3)
    T3 = AlgebraicPureTensor(E, F)
    result = S * T3
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_negation():
    """Test negation of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    S_neg = -S
    assert isinstance(S_neg, AlgebraicTensor)


def test_algebraic_tensor_transpose():
    """Test transpose of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    ST = S.T
    assert ST.shape == ((4, 3), (5, 4))


def test_algebraic_tensor_conjugate():
    """Test conjugate of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(1 + I, A, B)
    T2 = AlgebraicPureTensor(2 - I, C, D)
    S = AlgebraicTensor(T1, T2)
    
    SC = S.conjugate()
    assert isinstance(SC, AlgebraicTensor)
    
    # Double conjugate
    assert S.conjugate().conjugate() == S


def test_algebraic_tensor_diff():
    """Test differentiation of AlgebraicTensor."""
    from sympy.abc import x, y
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(x**2, A, B)
    T2 = AlgebraicPureTensor(x, C, D)
    S = AlgebraicTensor(T1, T2)
    
    S_diff = S.diff(x)
    assert isinstance(S_diff, AlgebraicTensor)
    
    # Differentiate with symbol not present
    S_diff2 = S.diff(y)
    assert isinstance(S_diff2, AlgebraicZeroTensor)


def test_algebraic_tensor_doit():
    """Test doit of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    assert S.doit() == S
    
    # With symbolic coefficients
    from sympy.abc import x, y
    T3 = AlgebraicPureTensor(x, A, B)
    T4 = AlgebraicPureTensor(y, A, B)
    S2 = AlgebraicTensor(T3, T4)
    S2_doit = S2.doit()
    # Coefficients should be combined
    assert S2_doit == AlgebraicTensor(T3, T4)


def test_algebraic_tensor_expand():
    """Test expand of AlgebraicTensor."""
    from sympy.matrices.expressions import MatAdd
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 4, 5)
    
    T = AlgebraicPureTensor(A, MatAdd(B, C))
    S = AlgebraicTensor(T)
    
    S_exp = S.expand()
    assert isinstance(S_exp, AlgebraicTensor)


def test_algebraic_tensor_simplify():
    """Test simplify of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    S_simp = S.simplify()
    assert isinstance(S_simp, (AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor))


def test_algebraic_tensor_pickle():
    """Test pickling of AlgebraicTensor."""
    import pickle
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    
    S_pickled = pickle.loads(pickle.dumps(S))
    assert S_pickled.shape == S.shape


def test_compose_algebraic_tensors():
    """Test compose_algebraic_tensors function."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 4, 2)
    D = MatrixSymbol('D', 5, 3)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    result = compose_algebraic_tensors(T1, T2)
    assert isinstance(result, AlgebraicPureTensor)
    
    # Compose two sums
    E = MatrixSymbol('E', 3, 4)
    F = MatrixSymbol('F', 4, 5)
    G = MatrixSymbol('G', 4, 2)
    H = MatrixSymbol('H', 5, 3)
    
    S1 = AlgebraicTensor(T1, AlgebraicPureTensor(E, F))
    S2 = AlgebraicTensor(T2, AlgebraicPureTensor(G, H))
    
    result2 = compose_algebraic_tensors(S1, S2)
    assert isinstance(result2, AlgebraicTensor)


def test_shape_mismatch_error():
    """Test ShapeMismatchError."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 2, 3)
    D = MatrixSymbol('D', 3, 4)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    raises(ShapeMismatchError, lambda: AlgebraicTensor(T1, T2))


def test_compose_shape_mismatch():
    """Test ShapeMismatchError on composition with incompatible inner dimensions."""
    from sympy.tensor.algebraic import compose_algebraic_pure_tensors
    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor

    # T1 has shape ((3, 4), (4, 5)), T2 has shape ((5, 2), (3, 3))
    # Factor 0: cols1=4 != rows2=5 -> mismatch
    # Factor 1: cols1=5 != rows2=3 -> mismatch
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 5, 2)
    D = MatrixSymbol('D', 3, 3)

    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)

    # Direct composition via function
    raises(ShapeMismatchError, lambda: compose_algebraic_pure_tensors(T1, T2))

    # Via __mul__
    raises(ShapeMismatchError, lambda: T1 * T2)

    # Via compose_algebraic_tensors
    raises(ShapeMismatchError, lambda: compose_algebraic_tensors(T1, T2))

    # With zero tensors: left zero, incompatible inner dims
    Z1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z2 = AlgebraicZeroTensor(((5, 2), (3, 3)))
    raises(ShapeMismatchError, lambda: compose_algebraic_tensors(Z1, Z2))
    raises(ShapeMismatchError, lambda: compose_algebraic_pure_tensors(Z1, Z2))

    # Zero left, non-zero right with incompatible inner dims
    raises(ShapeMismatchError, lambda: compose_algebraic_pure_tensors(Z1, T2))

    # Non-zero left, zero right with incompatible inner dims
    raises(ShapeMismatchError, lambda: compose_algebraic_pure_tensors(T1, Z2))

    # Different number of factors
    E = MatrixSymbol('E', 4, 2)
    T3 = AlgebraicPureTensor(A, E)  # shape ((3, 4), (4, 2))
    raises(ShapeMismatchError, lambda: compose_algebraic_pure_tensors(T1, T3))


def test_algebraic_tensor_coefficient_collection():
    """Test coefficient collection in AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, A, B)
    
    S = AlgebraicTensor(T1, T2)
    # Coefficients should be collected: 2 + 3 = 5
    assert S.coeff == 5


def test_algebraic_tensor_with_zero_tensor():
    """Test AlgebraicTensor with AlgebraicZeroTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    
    # Adding zero tensor to pure tensor returns the pure tensor
    S = T1 + Z
    assert S == T1


def test_algebraic_tensor_equality_commutative():
    """Test that AlgebraicTensor equality is commutative: T1+T2 == T2+T1."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    S1 = T1 + T2
    S2 = T2 + T1
    
    assert S1 == S2
    assert S2 == S1


def test_algebraic_tensor_equality_with_coefficients():
    """Test commutative equality with coefficients."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, C, D)
    
    S1 = T1 + T2
    S2 = T2 + T1
    
    assert S1 == S2


def test_algebraic_tensor_equality_different_terms():
    """Different terms are not equal."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    E = MatrixSymbol('E', 3, 4)
    F = MatrixSymbol('F', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(E, F)
    
    S1 = T1 + T2
    S2 = T1 + T3
    
    assert S1 != S2


# ---------------------------------------------------------------------------
# subs tests
# ---------------------------------------------------------------------------

def test_algebraic_tensor_subs_coefficient():
    """Substitute a symbol in coefficients of terms."""
    x = Symbol('x')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(x, A, B), AlgebraicPureTensor(2*x, C, D))
    result = S.subs(x, 3)
    # Should produce 3*A⊗B + 6*C⊗D
    assert any(t.coeff == 3 for t in result.args if hasattr(t, 'coeff'))
    assert any(t.coeff == 6 for t in result.args if hasattr(t, 'coeff'))


def test_algebraic_tensor_subs_factor():
    """Substitute a MatrixSymbol factor across all terms."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    E = MatrixSymbol('E', 3, 4)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = S.subs(A, E)
    assert AlgebraicPureTensor(E, B) in result.args
    assert AlgebraicPureTensor(C, D) in result.args
    assert A not in result.args


def test_algebraic_tensor_subs_multiple():
    """Substitute multiple symbols across terms."""
    x, y = Symbol('x'), Symbol('y')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(x, A, B), AlgebraicPureTensor(y, C, D))
    result = S.subs({x: 1, y: 2})
    assert any(t.coeff == 1 for t in result.args if hasattr(t, 'coeff'))
    assert any(t.coeff == 2 for t in result.args if hasattr(t, 'coeff'))


def test_algebraic_tensor_subs_no_match():
    """Substitution that doesn't match returns unchanged tensor."""
    x, y = Symbol('x'), Symbol('y')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(x, A, B), AlgebraicPureTensor(x, C, D))
    result = S.subs(y, 10)
    assert result == S


def test_algebraic_tensor_subs_concrete_matrix():
    """Substitute symbol inside concrete matrix entries in a sum."""
    x = Symbol('x')
    M = ImmutableDenseMatrix([[x, 1], [2, x]])
    N = MatrixSymbol('N', 2, 3)
    P = MatrixSymbol('P', 2, 2)
    Q = MatrixSymbol('Q', 2, 3)
    S = AlgebraicTensor(AlgebraicPureTensor(M, N), AlgebraicPureTensor(P, Q))
    result = S.subs(x, 7)
    for t in result.args:
        if hasattr(t, 'factors') and N in t.factors:
            assert t.factors[0][0, 0] == 7
            assert t.factors[0][1, 1] == 7


# ---------------------------------------------------------------------------
# free_symbols tests
# ---------------------------------------------------------------------------

def test_algebraic_tensor_free_symbols_basic():
    """free_symbols unions symbols from all terms."""
    x, y = Symbol('x'), Symbol('y')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(x, A, B), AlgebraicPureTensor(y, C, D))
    assert x in S.free_symbols
    assert y in S.free_symbols
    assert A in S.free_symbols
    assert D in S.free_symbols


def test_algebraic_tensor_free_symbols_shared_symbol():
    """free_symbols contains shared symbols only once."""
    x = Symbol('x')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(x, A, B), AlgebraicPureTensor(x, C, D))
    assert S.free_symbols == {x, A, B, C, D}


def test_algebraic_tensor_free_symbols_concrete_matrix():
    """free_symbols collects symbols from concrete matrices in terms."""
    x, y = Symbol('x'), Symbol('y')
    M = ImmutableDenseMatrix([[x, 1], [2, x]])
    N = MatrixSymbol('N', 2, 3)
    P = MatrixSymbol('P', 2, 2)
    Q = ImmutableDenseMatrix([[y, 0, 1], [0, 3, 2]])
    S = AlgebraicTensor(AlgebraicPureTensor(M, N), AlgebraicPureTensor(P, Q))
    assert x in S.free_symbols
    assert y in S.free_symbols


def test_algebraic_tensor_free_symbols_no_symbols():
    """Tensor with purely numeric content has no free symbols."""
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    N = ImmutableDenseMatrix([[5, 6, 7], [8, 9, 10]])
    P = ImmutableDenseMatrix([[2, 0], [0, 2]])
    Q = ImmutableDenseMatrix([[1, 1, 1], [1, 1, 1]])
    S = AlgebraicTensor(AlgebraicPureTensor(2, M, N), AlgebraicPureTensor(3, P, Q))
    assert S.free_symbols == set()


def test_algebraic_tensor_free_symbols_after_subs():
    """free_symbols updates correctly after substitution."""
    x, y = Symbol('x'), Symbol('y')
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(x, A, B), AlgebraicPureTensor(y, C, D))
    result = S.subs(x, 1)
    assert x not in result.free_symbols
    assert y in result.free_symbols


# ---------------------------------------------------------------------------
# _compose_with_term tests
# ---------------------------------------------------------------------------

def test_compose_with_term_basic():
    """Test _compose_with_term produces same result as compose_algebraic_tensors."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    E = MatrixSymbol('E', 4, 2)
    F = MatrixSymbol('F', 5, 3)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    T = AlgebraicPureTensor(E, F)
    result = S._compose_with_term(T)
    expected = compose_algebraic_tensors(S, T)
    assert result == expected


def test_compose_with_term_coefficient_merge():
    """Test _compose_with_term with terms that have same factors but different coefficients."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    E = MatrixSymbol('E', 4, 2)
    F = MatrixSymbol('F', 5, 3)
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, A, B)
    T3 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2, T3)
    T = AlgebraicPureTensor(E, F)
    result = S._compose_with_term(T)
    expected = compose_algebraic_tensors(S, T)
    assert result == expected


def test_compose_with_term_with_coefficients():
    """Test _compose_with_term with coefficients on terms."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    E = MatrixSymbol('E', 4, 2)
    F = MatrixSymbol('F', 5, 3)
    S = AlgebraicTensor(AlgebraicPureTensor(2, A, B), AlgebraicPureTensor(3, C, D))
    T = AlgebraicPureTensor(E, F)
    result = S._compose_with_term(T)
    assert isinstance(result, AlgebraicTensor)


def test_compose_with_term_with_zero_tensor():
    """Test _compose_with_term handles zero tensor terms in the sum."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    E = MatrixSymbol('E', 4, 2)
    F = MatrixSymbol('F', 5, 3)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D), Z)
    T = AlgebraicPureTensor(E, F)
    result = S._compose_with_term(T)
    assert isinstance(result, (AlgebraicTensor, AlgebraicPureTensor))


# ---------------------------------------------------------------------------
# is_Add flag tests
# ---------------------------------------------------------------------------

def test_is_add_flag():
    """Test that AlgebraicTensor has is_Add = True."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    assert S.is_Add is True


def test_is_add_flag_basic_inheritance():
    """Test that AlgebraicTensor with is_Add works as Basic subclass."""
    from sympy.core.basic import Basic
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    S = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    assert isinstance(S, Basic)
    assert S.is_Add is True


def test_is_not_add_for_pure_tensor():
    """Test that AlgebraicPureTensor does not have is_Add = True."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    T = AlgebraicPureTensor(A, B)
    assert getattr(T, 'is_Add', None) is not True


def test_is_not_add_for_zero_tensor():
    """Test that AlgebraicZeroTensor does not have is_Add = True."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert getattr(Z, 'is_Add', None) is not True
