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
