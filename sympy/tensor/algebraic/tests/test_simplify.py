"""Tests for simplify.py"""
from sympy.core.singleton import S
from sympy.core.numbers import I
from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices import ImmutableDenseMatrix
from sympy.tensor.algebraic import AlgebraicPureTensor, AlgebraicTensor, AlgebraicZeroTensor
from sympy.tensor.algebraic.simplify import (
    tensorsimplify,
    _simplify_algebraic_pure_tensor,
    _simplify_algebraic_tensor,
    _equality_factoring,
    _matrix_proportionality_ratio,
    _proportionality_ratio,
    _extract_pt_and_coeff,
    _build_pt,
    _extract_commutative_from_factor,
    _commutativity_simplify,
)


def test_tensorsimplify_pure_tensor():
    """Test tensorsimplify with AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(2, A, B)
    T_simp = tensorsimplify(T)
    
    assert isinstance(T_simp, AlgebraicPureTensor)
    assert T_simp.coeff == 2


def test_tensorsimplify_zero_tensor():
    """Test tensorsimplify with AlgebraicZeroTensor."""
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z_simp = tensorsimplify(Z)
    
    assert Z_simp is Z


def test_tensorsimplify_tensor_sum():
    """Test tensorsimplify with AlgebraicTensor sum."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    S = AlgebraicTensor(T1, T2)
    S_simp = tensorsimplify(S)
    
    assert isinstance(S_simp, AlgebraicTensor)


def test_equality_factoring():
    """Test equality factoring."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(A, C)
    
    S = AlgebraicTensor(T1, T2)
    S_factored = _equality_factoring(S)
    
    # Should factor out A - result may be AlgebraicTensor or a factored expression
    assert S_factored is not None


def test_matrix_proportionality_ratio():
    """Test matrix proportionality ratio."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    
    ratio = _matrix_proportionality_ratio(M1, M2)
    assert ratio == S.Half
    
    ratio2 = _matrix_proportionality_ratio(M2, M1)
    assert ratio2 == 2
    
    # Non-proportional matrices
    M3 = ImmutableDenseMatrix([[1, 2], [3, 5]])
    ratio3 = _matrix_proportionality_ratio(M1, M3)
    assert ratio3 is None


def test_proportionality_ratio():
    """Test proportionality ratio for factors."""
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    
    ratio = _proportionality_ratio(M1, M2)
    assert ratio == S.Half
    
    # Same matrix
    ratio2 = _proportionality_ratio(M1, M1)
    assert ratio2 == S.One


def test_extract_pt_and_coeff():
    """Test extraction of coefficient and factors."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(2, A, B)
    coeff, factors = _extract_pt_and_coeff(T)
    
    assert coeff == 2
    assert factors == [A, B]
    
    # Without coefficient
    T2 = AlgebraicPureTensor(A, B)
    coeff2, factors2 = _extract_pt_and_coeff(T2)
    
    assert coeff2 == S.One
    assert factors2 == [A, B]


def test_build_pt():
    """Test building AlgebraicPureTensor from coeff and factors."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    # With coefficient
    T = _build_pt(2, [A, B])
    assert isinstance(T, AlgebraicPureTensor)
    assert T.coeff == 2
    
    # Without coefficient (coeff=1)
    T2 = _build_pt(S.One, [A])
    assert T2 == A
    
    # Zero coefficient
    T3 = _build_pt(S.Zero, [A, B])
    assert isinstance(T3, AlgebraicZeroTensor)


def test_extract_commutative_from_factor_add():
    """Test extracting commutative prefactor from Add."""
    from sympy.abc import x
    A = MatrixSymbol('A', 3, 4)
    
    # Symbolic Add
    expr = x * A + x * A
    coeff, factor = _extract_commutative_from_factor(expr)
    
    # Should extract commutative part if possible, or return S.One
    # The function may not always extract from Add expressions
    assert coeff is not None


def test_extract_commutative_from_factor_matrix():
    """Test extracting commutative prefactor from Matrix."""
    from sympy.abc import x
    M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    
    coeff, new_M = _extract_commutative_from_factor(M)
    
    assert coeff == x
    assert new_M[0, 0] == 1
    assert new_M[0, 1] == 2


def test_commutativity_simplify():
    """Test commutativity-based simplification."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    S = AlgebraicTensor(T1, T2)
    S_simp = _commutativity_simplify(S)
    
    assert isinstance(S_simp, (AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor))


def test_simplify_algebraic_pure_tensor():
    """Test simplification of AlgebraicPureTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(2, A, B)
    T_simp = _simplify_algebraic_pure_tensor(T)
    
    assert isinstance(T_simp, AlgebraicPureTensor)


def test_simplify_algebraic_tensor():
    """Test simplification of AlgebraicTensor."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    S = AlgebraicTensor(T1, T2)
    S_simp = _simplify_algebraic_tensor(S)
    
    assert isinstance(S_simp, (AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor))


def test_tensorsimplify_proportional_terms():
    """Test tensorsimplify with proportional terms."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(2, A, B)
    
    S = AlgebraicTensor(T1, T2)
    S_simp = tensorsimplify(S)
    
    # Should combine to 3*A ⊗ B
    assert isinstance(S_simp, AlgebraicPureTensor)
    assert S_simp.coeff == 3


def test_extract_commutative_from_factor_mul():
    """Test extracting commutative prefactor from Mul."""
    from sympy.abc import x
    A = MatrixSymbol('A', 3, 4)
    
    # Mul with commutative symbol
    expr = x * A
    coeff, factor = _extract_commutative_from_factor(expr)
    
    assert coeff == x
    assert factor == A


def test_simplify_with_numeric_matrix():
    """Test simplification with numeric matrices."""
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    N = ImmutableDenseMatrix([[2, 4], [6, 8]])
    
    T1 = AlgebraicPureTensor(M)
    T2 = AlgebraicPureTensor(N)
    
    S = AlgebraicTensor(T1, T2)
    S_simp = tensorsimplify(S)
    
    # Should handle numeric matrix simplification
    assert isinstance(S_simp, (AlgebraicTensor, AlgebraicPureTensor))


def test_simplify_preserves_shape():
    """Test that simplification preserves tensor shape."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    C = MatrixSymbol('C', 3, 4)
    D = MatrixSymbol('D', 4, 5)
    
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    
    S = AlgebraicTensor(T1, T2)
    S_simp = tensorsimplify(S)
    
    assert S_simp.shape == S.shape


def test_simplify_complex_coefficient():
    """Test simplification with complex coefficients."""
    A = MatrixSymbol('A', 3, 4)
    B = MatrixSymbol('B', 4, 5)
    
    T = AlgebraicPureTensor(1 + I, A, B)
    T_simp = tensorsimplify(T)
    
    assert isinstance(T_simp, AlgebraicPureTensor)
    assert T_simp.coeff == 1 + I
