"""Tests for algebraic_tensor.py"""
from sympy.core.numbers import I
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
    assert isinstance(result3, (AlgebraicTensor, AlgebraicPureTensor))


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
