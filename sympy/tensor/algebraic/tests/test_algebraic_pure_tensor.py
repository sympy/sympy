from __future__ import annotations

import pickle

from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.matrices.expressions import MatrixSymbol, MatAdd, MatMul
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    compose_algebraic_pure_tensors,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 4, 5)
C = MatrixSymbol("C", 3, 4)
D = MatrixSymbol("D", 4, 5)
E = MatrixSymbol("E", 4, 2)
F = MatrixSymbol("F", 5, 3)
G = MatrixSymbol("G", 5, 3)
H = MatrixSymbol("H", 3, 3)

x = Symbol("x")
y = Symbol("y")


# ---------------------------------------------------------------------------
# Type flags and attributes
# ---------------------------------------------------------------------------

def test_type_flag_is_AlgebraicPureTensor():
    T = AlgebraicPureTensor(A, B)
    assert T.is_AlgebraicPureTensor is True


def test_type_flag_is_Mul():
    T = AlgebraicPureTensor(A, B)
    assert T.is_Mul is False


def test_type_flag_is_commutative():
    T = AlgebraicPureTensor(A, B)
    assert T.is_commutative is False


def test_op_priority():
    T = AlgebraicPureTensor(A, B)
    assert T._op_priority == 11


# ---------------------------------------------------------------------------
# Constructor: basic
# ---------------------------------------------------------------------------

def test_constructor_two_factors():
    T = AlgebraicPureTensor(A, B)
    assert isinstance(T, AlgebraicPureTensor)
    assert T.shape == ((3, 4), (4, 5))


def test_constructor_three_factors():
    T = AlgebraicPureTensor(A, B, G)
    assert isinstance(T, AlgebraicPureTensor)
    assert T.shape == ((3, 4), (4, 5), (5, 3))


def test_constructor_numeric_coeff():
    T = AlgebraicPureTensor(2, A, B)
    assert isinstance(T, AlgebraicPureTensor)
    assert T.coeff == 2
    assert T.factors == (A, B)


def test_constructor_symbolic_coeff():
    T = AlgebraicPureTensor(x, A, B)
    assert isinstance(T, AlgebraicPureTensor)
    assert T.coeff == x
    assert T.factors == (A, B)


def test_constructor_negative_coeff():
    T = AlgebraicPureTensor(-1, A, B)
    assert isinstance(T, AlgebraicPureTensor)
    assert T.coeff == -1
    assert T.factors == (A, B)


def test_constructor_zero_coeff():
    T = AlgebraicPureTensor(0, A, B)
    assert isinstance(T, AlgebraicZeroTensor)
    assert T.shape == ((3, 4), (4, 5))


def test_constructor_zero_coeff_symbolic():
    T = AlgebraicPureTensor(S.Zero, A, B)
    assert isinstance(T, AlgebraicZeroTensor)
    assert T.shape == ((3, 4), (4, 5))


def test_constructor_single_factor_unwrap():
    result = AlgebraicPureTensor(A)
    assert result is A


def test_constructor_single_factor_with_coeff():
    """Two args with commutative first arg returns coeff * factor (MatMul)."""
    from sympy.matrices.expressions.matexpr import MatMul
    T = AlgebraicPureTensor(2, A)
    assert isinstance(T, MatMul)
    assert T == 2 * A


def test_constructor_single_factor_coeff_one():
    result = AlgebraicPureTensor(1, A)
    assert result is A


def test_constructor_coeff_S_One():
    result = AlgebraicPureTensor(S.One, A)
    assert result is A


def test_constructor_coeff_S_One_multi():
    T = AlgebraicPureTensor(S.One, A, B)
    assert isinstance(T, AlgebraicPureTensor)
    assert T.coeff == S.One
    assert T.factors == (A, B)


# ---------------------------------------------------------------------------
# Constructor: errors
# ---------------------------------------------------------------------------

def test_constructor_no_args():
    raises(ValueError, lambda: AlgebraicPureTensor())


def test_constructor_only_coeff():
    """Single arg returns the arg directly."""
    assert AlgebraicPureTensor(2) == 2


def test_constructor_only_zero_coeff():
    """Single arg 0 returns 0 directly."""
    assert AlgebraicPureTensor(0) == 0


def test_constructor_scalar_as_factor():
    raises(TypeError, lambda: AlgebraicPureTensor(A, 3))


def test_constructor_no_shape():
    """Single arg without shape returns the arg directly."""
    z = Symbol("z")
    assert AlgebraicPureTensor(z) is z


# ---------------------------------------------------------------------------
# coeff property
# ---------------------------------------------------------------------------

def test_coeff_no_coeff():
    T = AlgebraicPureTensor(A, B)
    assert T.coeff == S.One


def test_coeff_numeric():
    T = AlgebraicPureTensor(2, A, B)
    assert T.coeff == 2


def test_coeff_symbolic():
    T = AlgebraicPureTensor(x, A, B)
    assert T.coeff == x


def test_coeff_negative():
    T = AlgebraicPureTensor(-3, A, B)
    assert T.coeff == -3


# ---------------------------------------------------------------------------
# factors property
# ---------------------------------------------------------------------------

def test_factors_no_coeff():
    T = AlgebraicPureTensor(A, B)
    assert T.factors == (A, B)


def test_factors_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    assert T.factors == (A, B)


def test_factors_symbolic_coeff():
    T = AlgebraicPureTensor(x, A, B)
    assert T.factors == (A, B)


def test_factors_single():
    """Two args with commutative first returns MatMul, not PureTensor."""
    from sympy.matrices.expressions.matexpr import MatMul
    T = AlgebraicPureTensor(2, A)
    assert isinstance(T, MatMul)


def test_factors_three():
    T = AlgebraicPureTensor(x, A, B, G)
    assert T.factors == (A, B, G)


# ---------------------------------------------------------------------------
# num_factors property
# ---------------------------------------------------------------------------

def test_num_factors_two():
    T = AlgebraicPureTensor(A, B)
    assert T.num_factors == 2


def test_num_factors_one():
    """Two args with commutative first returns MatMul, not PureTensor."""
    from sympy.matrices.expressions.matexpr import MatMul
    T = AlgebraicPureTensor(2, A)
    assert isinstance(T, MatMul)


def test_num_factors_three():
    T = AlgebraicPureTensor(A, B, G)
    assert T.num_factors == 3


# ---------------------------------------------------------------------------
# shape property
# ---------------------------------------------------------------------------

def test_shape_two_factors():
    T = AlgebraicPureTensor(A, B)
    assert T.shape == ((3, 4), (4, 5))


def test_shape_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    assert T.shape == ((3, 4), (4, 5))


def test_shape_single_factor():
    T = AlgebraicPureTensor(2, A)
    # Two args with commutative first returns MatMul, shape is bare (m, n)
    assert T.shape == (3, 4)


def test_shape_three_factors():
    T = AlgebraicPureTensor(A, B, G)
    assert T.shape == ((3, 4), (4, 5), (5, 3))


# ---------------------------------------------------------------------------
# commutativity_pattern
# ---------------------------------------------------------------------------

def test_commutativity_pattern_all_symbolic():
    T = AlgebraicPureTensor(A, B)
    assert T.commutativity_pattern == (0, 0)


def test_commutativity_pattern_numeric_matrix():
    M = ImmutableDenseMatrix([[1, 2, 3, 4],
                               [5, 6, 7, 8],
                               [9, 10, 11, 12]])
    T = AlgebraicPureTensor(M, B)
    assert T.commutativity_pattern == (1, 0)


def test_commutativity_pattern_all_numeric():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[5, 6], [7, 8]])
    T = AlgebraicPureTensor(M1, M2)
    assert T.commutativity_pattern == (1, 1)


def test_commutativity_pattern_single():
    """Two args with commutative first returns MatMul, not PureTensor."""
    from sympy.matrices.expressions.matexpr import MatMul
    T = AlgebraicPureTensor(2, A)
    assert isinstance(T, MatMul)


def test_commutativity_pattern_three_factors():
    M = ImmutableDenseMatrix([[1, 2, 3, 4],
                               [5, 6, 7, 8],
                               [9, 10, 11, 12]])
    T = AlgebraicPureTensor(M, B, G)
    assert T.commutativity_pattern == (1, 0, 0)


# ---------------------------------------------------------------------------
# _get_coeff
# ---------------------------------------------------------------------------

def test_get_coeff():
    T = AlgebraicPureTensor(2, A, B)
    assert T._get_coeff() == 2


def test_get_coeff_default():
    T = AlgebraicPureTensor(A, B)
    assert T._get_coeff() == S.One


# ---------------------------------------------------------------------------
# Negation
# ---------------------------------------------------------------------------

def test_neg_no_coeff():
    T = AlgebraicPureTensor(A, B)
    nT = -T
    assert isinstance(nT, AlgebraicPureTensor)
    assert nT.coeff == S.NegativeOne
    assert nT.factors == (A, B)


def test_neg_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    nT = -T
    assert isinstance(nT, AlgebraicPureTensor)
    assert nT.coeff == -2
    assert nT.factors == (A, B)


def test_neg_negative_coeff():
    T = AlgebraicPureTensor(-3, A, B)
    nT = -T
    assert isinstance(nT, AlgebraicPureTensor)
    assert nT.coeff == 3
    assert nT.factors == (A, B)


def test_double_neg():
    T = AlgebraicPureTensor(2, A, B)
    assert -(-T) == T


# ---------------------------------------------------------------------------
# Scalar multiplication (__mul__)
# ---------------------------------------------------------------------------

def test_mul_numeric():
    T = AlgebraicPureTensor(A, B)
    result = T * 3
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 3
    assert result.factors == (A, B)


def test_mul_with_existing_coeff():
    T = AlgebraicPureTensor(2, A, B)
    result = T * 3
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6
    assert result.factors == (A, B)


def test_mul_one():
    T = AlgebraicPureTensor(A, B)
    result = T * 1
    assert result is T


def test_mul_S_One():
    T = AlgebraicPureTensor(A, B)
    result = T * S.One
    assert result is T


def test_mul_zero():
    T = AlgebraicPureTensor(A, B)
    result = T * 0
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_mul_S_Zero():
    T = AlgebraicPureTensor(A, B)
    result = T * S.Zero
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_mul_symbol():
    T = AlgebraicPureTensor(A, B)
    result = T * x
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x
    assert result.factors == (A, B)


def test_mul_symbol_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    result = T * x
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2 * x
    assert result.factors == (A, B)


# ---------------------------------------------------------------------------
# Scalar right-multiplication (__rmul__)
# ---------------------------------------------------------------------------

def test_rmul_numeric():
    T = AlgebraicPureTensor(A, B)
    result = 3 * T
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 3
    assert result.factors == (A, B)


def test_rmul_with_existing_coeff():
    T = AlgebraicPureTensor(2, A, B)
    result = 3 * T
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6
    assert result.factors == (A, B)


def test_rmul_one():
    T = AlgebraicPureTensor(A, B)
    result = 1 * T
    assert result is T


def test_rmul_S_One():
    T = AlgebraicPureTensor(A, B)
    result = S.One * T
    assert result is T


def test_rmul_zero():
    T = AlgebraicPureTensor(A, B)
    result = 0 * T
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_rmul_S_Zero():
    T = AlgebraicPureTensor(A, B)
    result = S.Zero * T
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_rmul_symbol():
    T = AlgebraicPureTensor(A, B)
    result = x * T
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x
    assert result.factors == (A, B)


def test_rmul_symbol_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    result = x * T
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2 * x
    assert result.factors == (A, B)


# ---------------------------------------------------------------------------
# Composition via __mul__ (non-commutative)
# ---------------------------------------------------------------------------

def test_mul_pure_tensor_compose():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(E, F)
    result = T1 * T2
    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 2


def test_mul_compose_with_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, E, F)
    result = T1 * T2
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6


def test_mul_compose_result_shapes():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(E, F)
    result = T1 * T2
    assert result.shape == ((3, 2), (4, 3))


def test_mul_bare_matrix_compose():
    T = AlgebraicPureTensor(A, B)
    # T has 2 factors, E is a bare matrix (1 factor) -> factor count mismatch
    raises(ValueError, lambda: T * E)


def test_rmul_pure_tensor_compose():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(E, F)
    # T2*T1: E(4,2)*A(3,4) has incompatible inner dims (2 vs 3)
    raises(ValueError, lambda: T2 * T1)


def test_rmul_bare_matrix_compose():
    T = AlgebraicPureTensor(A, B)
    # A * T: MatrixSymbol.__mul__ creates MatMul, which rejects noncommutative AlgebraicPureTensor
    raises(NotImplementedError, lambda: A * T)


# ---------------------------------------------------------------------------
# Composition with AlgebraicZeroTensor
# ---------------------------------------------------------------------------

def test_mul_zero_tensor():
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = T * Z
    assert isinstance(result, AlgebraicZeroTensor)


def test_rmul_zero_tensor():
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z * T
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# Addition
# ---------------------------------------------------------------------------

def test_add_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    result = T1 + T2
    assert isinstance(result, AlgebraicTensor)


def test_add_with_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(C, D)
    result = T1 + T2
    assert isinstance(result, AlgebraicTensor)


def test_add_zero_tensor_same_shape():
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = T + Z
    assert result is T


def test_radd_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    result = T2 + T1
    assert isinstance(result, AlgebraicTensor)


def test_radd_zero_tensor_same_shape():
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = Z + T
    assert result is T


# ---------------------------------------------------------------------------
# Subtraction
# ---------------------------------------------------------------------------

def test_sub_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    result = T1 - T2
    assert isinstance(result, AlgebraicTensor)


def test_sub_self():
    T = AlgebraicPureTensor(A, B)
    result = T - T
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_rsub_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    result = T2 - T1
    assert isinstance(result, AlgebraicTensor)


def test_rsub_self():
    T = AlgebraicPureTensor(A, B)
    result = T - T
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# has_zero_term
# ---------------------------------------------------------------------------

def test_has_zero_term_false():
    T = AlgebraicPureTensor(A, B)
    assert T.has_zero_term() is False


def test_has_zero_term_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    assert T.has_zero_term() is False


# ---------------------------------------------------------------------------
# Transpose
# ---------------------------------------------------------------------------

def test_transpose_basic():
    T = AlgebraicPureTensor(A, B)
    TT = T.T
    assert isinstance(TT, AlgebraicPureTensor)
    assert TT.shape == ((4, 3), (5, 4))
    assert TT.coeff == S.One


def test_transpose_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    TT = T.T
    assert isinstance(TT, AlgebraicPureTensor)
    assert TT.shape == ((4, 3), (5, 4))
    assert TT.coeff == 2


def test_transpose_single_factor():
    """Two args with commutative first returns MatMul, .T gives 3*A.T."""
    T = AlgebraicPureTensor(3, A)
    TT = T.T
    assert TT == 3 * A.T


def test_transpose_double():
    T = AlgebraicPureTensor(A, B)
    assert T.T.T.shape == T.shape


def test_transpose_preserves_factors():
    T = AlgebraicPureTensor(A, B)
    TT = T.T
    assert TT.factors[0] == A.T
    assert TT.factors[1] == B.T


# ---------------------------------------------------------------------------
# Expand: _eval_expand_mul
# ---------------------------------------------------------------------------

def test_expand_no_add():
    T = AlgebraicPureTensor(A, B)
    result = T.expand()
    assert result is T


def test_expand_factor_with_matadd():
    T = AlgebraicPureTensor(A, MatAdd(B, D))
    result = T.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_coeff_with_add():
    T = AlgebraicPureTensor(x + y, A, B)
    result = T.expand()
    # Terms share the same factors, so they're combined back into AlgebraicPureTensor
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x + y


def test_expand_multiple_add_factors():
    T = AlgebraicPureTensor(MatAdd(A, C), MatAdd(B, D))
    result = T.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_preserves_coeff():
    T = AlgebraicPureTensor(2, A, MatAdd(B, D))
    result = T.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_single_factor_add():
    """Two args with commutative first returns MatMul, expand distributes."""
    T = AlgebraicPureTensor(2, MatAdd(A, C))
    result = T.expand()
    assert result == 2 * A + 2 * C


def test_expand_no_change_returns_self():
    T = AlgebraicPureTensor(2, A, B)
    result = T.expand()
    assert result is T


def test_expand_deep_hint():
    T = AlgebraicPureTensor(A, MatAdd(B, D))
    result = T.expand(deep=True)
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# compose_algebraic_pure_tensors
# ---------------------------------------------------------------------------

def test_compose_basic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_pure_tensors(T1, T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 2
    assert result.shape == ((3, 2), (4, 3))


def test_compose_with_coefficients():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, E, F)
    result = compose_algebraic_pure_tensors(T1, T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6


def test_compose_coeff_from_one_side():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_pure_tensors(T1, T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2


def test_compose_bare_matrices():
    result = compose_algebraic_pure_tensors(A, E)
    assert isinstance(result, MatMul)
    assert result.shape == (3, 2)


def test_compose_pure_with_bare():
    T = AlgebraicPureTensor(A, B)
    # E is 4x2 (1 factor), T is 2 factors -> mismatch
    raises(ValueError, lambda: compose_algebraic_pure_tensors(T, E))


def test_compose_single_factor_unwrap():
    """Two args with commutative first returns MatMul; compose still works."""
    T1 = AlgebraicPureTensor(2, A)
    T2 = AlgebraicPureTensor(E)
    result = compose_algebraic_pure_tensors(T1, T2)
    # T1 is MatMul(2, A), T2 is E (bare matrix)
    # compose gives MatMul(MatMul(2, A), E)
    assert result == MatMul(T1, E, evaluate=False)


def test_compose_single_factor_coeff_one_unwrap():
    T1 = AlgebraicPureTensor(A)
    T2 = AlgebraicPureTensor(E)
    result = compose_algebraic_pure_tensors(T1, T2)
    # Single factor, coeff 1 -> unwraps to bare MatMul
    assert isinstance(result, MatMul)


def test_compose_factor_count_mismatch():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(A)
    raises(ValueError, lambda: compose_algebraic_pure_tensors(T1, T2))


def test_compose_inner_dim_mismatch():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(A, B)
    # A is 3x4, A is 3x4 -> inner dims 4 vs 3 don't match for slot 0
    raises(ValueError, lambda: compose_algebraic_pure_tensors(T1, T2))


def test_compose_invalid_left_type():
    raises(TypeError, lambda: compose_algebraic_pure_tensors(x, A))


def test_compose_invalid_right_type():
    T = AlgebraicPureTensor(A, B)
    raises(TypeError, lambda: compose_algebraic_pure_tensors(T, x))


# ---------------------------------------------------------------------------
# compose_algebraic_pure_tensors with AlgebraicZeroTensor
# ---------------------------------------------------------------------------

def test_compose_both_zero():
    Z1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z2 = AlgebraicZeroTensor(((4, 2), (5, 3)))
    result = compose_algebraic_pure_tensors(Z1, Z2)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 2), (4, 3))


def test_compose_left_zero():
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    T = AlgebraicPureTensor(E, F)
    result = compose_algebraic_pure_tensors(Z, T)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 2), (4, 3))


def test_compose_right_zero():
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((4, 2), (5, 3)))
    result = compose_algebraic_pure_tensors(T, Z)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 2), (4, 3))


def test_compose_zero_with_bare_matrix():
    Z = AlgebraicZeroTensor(((3, 4),))
    result = compose_algebraic_pure_tensors(Z, E)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 2),)


def test_compose_bare_matrix_with_zero():
    Z = AlgebraicZeroTensor(((4, 2),))
    result = compose_algebraic_pure_tensors(A, Z)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 2),)


def test_compose_both_zero_factor_mismatch():
    Z1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z2 = AlgebraicZeroTensor(((3, 4),))
    raises(ValueError, lambda: compose_algebraic_pure_tensors(Z1, Z2))


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_constructor_sympified_coeff():
    T = AlgebraicPureTensor("2", A, B)
    assert isinstance(T, AlgebraicPureTensor)
    assert T.coeff == 2


def test_mul_commutative_dispatch():
    """Ensure commutative * PureTensor goes through __rmul__, not composition."""
    T = AlgebraicPureTensor(A, B)
    result = x * T
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x


def test_mul_noncommutative_dispatch():
    """Ensure PureTensor * PureTensor goes through composition."""
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(E, F)
    result = T1 * T2
    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 2


def test_identity_preserved_through_neg():
    T = AlgebraicPureTensor(2, A, B)
    nT = -T
    assert isinstance(nT, AlgebraicPureTensor)


def test_identity_preserved_through_mul():
    T = AlgebraicPureTensor(A, B)
    result = T * 3
    assert isinstance(result, AlgebraicPureTensor)


def test_identity_preserved_through_transpose():
    T = AlgebraicPureTensor(A, B)
    result = T.T
    assert isinstance(result, AlgebraicPureTensor)


def test_pickle():
    T = AlgebraicPureTensor(2, A, B)
    data = pickle.dumps(T)
    T2 = pickle.loads(data)
    assert T2 == T
    assert T2.coeff == 2
    assert T2.factors == (A, B)


def test_hash_consistency():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(A, B)
    assert hash(T1) == hash(T2)
    assert T1 == T2


def test_set_membership():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(A, B)
    T3 = AlgebraicPureTensor(C, D)
    s = {T1, T2, T3}
    assert len(s) == 2


def test_args_structure():
    T = AlgebraicPureTensor(2, A, B)
    assert len(T.args) == 3
    assert T.args[0] == 2
    assert T.args[1] == A
    assert T.args[2] == B


def test_args_no_coeff():
    T = AlgebraicPureTensor(A, B)
    assert len(T.args) == 2
    assert T.args[0] == A
    assert T.args[1] == B


# ---------------------------------------------------------------------------
# simplify
# ---------------------------------------------------------------------------

def test_simplify_basic():
    T = AlgebraicPureTensor(2, A, B)
    result = T.simplify()
    assert isinstance(result, AlgebraicPureTensor)


def test_simplify_symbolic_coeff():
    T = AlgebraicPureTensor(x + x, A, B)
    result = T.simplify()
    assert isinstance(result, AlgebraicPureTensor)


# ---------------------------------------------------------------------------
# display (smoke test)
# ---------------------------------------------------------------------------

def test_display_text_mode():
    T = AlgebraicPureTensor(A, B)
    # Should not raise
    T.display(mode="text")


def test_display_latex_mode():
    T = AlgebraicPureTensor(A, B)
    # Should not raise
    T.display(mode="latex")
