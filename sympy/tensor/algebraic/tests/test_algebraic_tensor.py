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
    compose_algebraic_tensors,
)
from sympy.tensor.algebraic.algebraic_tensor import (
    ShapeMismatchError,
    _commutativity_pattern_of,
    _shape_of,
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
G = MatrixSymbol("G", 3, 4)
H = MatrixSymbol("H", 4, 5)
I = MatrixSymbol("I", 5, 3)

# Matrices with compatible dimensions for composition:
# A(3,4) can compose with E(4,2) -> (3,2)
# B(4,5) can compose with F(5,3) -> (4,3)
# J(4,2) can compose with K(2,1) -> (4,1)
# L(5,3) can compose with M(3,2) -> (5,2)
J = MatrixSymbol("J", 4, 2)
K = MatrixSymbol("K", 2, 1)
L = MatrixSymbol("L", 5, 3)
M = MatrixSymbol("M", 3, 2)

x = Symbol("x")
y = Symbol("y")


# ---------------------------------------------------------------------------
# _shape_of helper
# ---------------------------------------------------------------------------

def test_shape_of_pure_tensor():
    T = AlgebraicPureTensor(A, B)
    assert _shape_of(T) == ((3, 4), (4, 5))


def test_shape_of_algebraic_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert _shape_of(S) == ((3, 4), (4, 5))


def test_shape_of_zero_tensor():
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert _shape_of(Z) == ((3, 4), (4, 5))


def test_shape_of_bare_matrix():
    assert _shape_of(A) == ((3, 4),)


def test_shape_of_mul_coeff_pure_tensor():
    T = AlgebraicPureTensor(2, A, B)
    mul_expr = x * T
    assert _shape_of(mul_expr) == ((3, 4), (4, 5))


def test_shape_of_unknown():
    assert _shape_of(x) is None


# ---------------------------------------------------------------------------
# _commutativity_pattern_of helper
# ---------------------------------------------------------------------------

def test_commutativity_pattern_of_pure_tensor():
    T = AlgebraicPureTensor(A, B)
    assert _commutativity_pattern_of(T) == (0, 0)


def test_commutativity_pattern_of_zero_tensor():
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    assert _commutativity_pattern_of(Z) == (1, 1)


def test_commutativity_pattern_of_algebraic_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert _commutativity_pattern_of(S) == (0, 0)


def test_commutativity_pattern_of_bare_matrix():
    assert _commutativity_pattern_of(A) == (0,)


def test_commutativity_pattern_of_bare_numeric_matrix():
    M = ImmutableDenseMatrix([[1, 2, 3, 4],
                               [5, 6, 7, 8],
                               [9, 10, 11, 12]])
    assert _commutativity_pattern_of(M) == (1,)


def test_commutativity_pattern_of_mul_coeff_pure_tensor():
    T = AlgebraicPureTensor(2, A, B)
    mul_expr = x * T
    assert _commutativity_pattern_of(mul_expr) == (0, 0)


def test_commutativity_pattern_of_unknown():
    assert _commutativity_pattern_of(x) is None


# ---------------------------------------------------------------------------
# ShapeMismatchError
# ---------------------------------------------------------------------------

def test_shape_mismatch_error_is_type_error():
    assert issubclass(ShapeMismatchError, TypeError)


def test_shape_mismatch_raised():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(MatrixSymbol("X", 2, 3),
                              MatrixSymbol("Y", 3, 4))
    raises(ShapeMismatchError, lambda: AlgebraicTensor(T1, T2))


# ---------------------------------------------------------------------------
# AlgebraicTensor: type flags
# ---------------------------------------------------------------------------

def test_type_flag_is_AlgebraicTensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.is_AlgebraicTensor is True


def test_type_flag_is_Add():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.is_Add is True


def test_identity_none():
    assert AlgebraicTensor.identity is None


# ---------------------------------------------------------------------------
# AlgebraicTensor: constructor
# ---------------------------------------------------------------------------

def test_constructor_basic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert isinstance(S, AlgebraicTensor)
    assert len(S.args) == 2


def test_constructor_single_pure_tensor_unwrap():
    T = AlgebraicPureTensor(A, B)
    result = AlgebraicTensor(T)
    assert result is T


def test_constructor_single_zero_tensor_unwrap():
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = AlgebraicTensor(Z)
    assert result is Z


def test_constructor_single_algebraic_tensor_unwrap():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = AlgebraicTensor(S)
    assert result is S


def test_constructor_no_args():
    raises(ValueError, lambda: AlgebraicTensor())


def test_constructor_cancellation_to_zero():
    T = AlgebraicPureTensor(A, B)
    result = AlgebraicTensor(T, -T)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_constructor_flat_nested():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(S1, T3)
    assert isinstance(S2, AlgebraicTensor)
    assert len(S2.args) == 3


def test_constructor_with_zero_tensor_anchor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    assert isinstance(S, AlgebraicTensor)
    assert S.has_zero_term()


def test_constructor_coefficient_collection():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, A, B)
    S = AlgebraicTensor(T1, T2)
    assert isinstance(S, AlgebraicPureTensor)
    assert S.coeff == 5
    assert S.factors == (A, B)


def test_constructor_coefficient_collection_partial():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, A, B)
    T3 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2, T3)
    assert isinstance(S, AlgebraicTensor)


def test_constructor_with_mul_coeff_terms():
    """Mul(coeff, PureTensor) terms are handled correctly."""
    T = AlgebraicPureTensor(A, B)
    mul_term = 2 * T
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(mul_term, T2)
    assert isinstance(S, AlgebraicTensor)


def test_constructor_single_term_unwrap_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    result = AlgebraicTensor(T)
    assert result is T


def test_constructor_with_zero_tensor_keeps_anchor():
    """When user explicitly provides a zero tensor, it stays as anchor."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T, Z)
    assert isinstance(S, AlgebraicTensor)
    assert S.has_zero_term()


def test_constructor_cancellation_keeps_anchor_when_user_provided():
    """When zero tensor was user-provided, cancellation keeps it."""
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T, -T, Z)
    assert isinstance(S, AlgebraicZeroTensor)
    assert S.shape == ((3, 4), (4, 5))


def test_constructor_list_single_element():
    T = AlgebraicPureTensor(A, B)
    result = AlgebraicTensor([T])
    assert result is T


# ---------------------------------------------------------------------------
# AlgebraicTensor: shape property
# ---------------------------------------------------------------------------

def test_shape_basic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.shape == ((3, 4), (4, 5))


def test_shape_with_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.shape == ((3, 4), (4, 5))


def test_shape_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    assert S.shape == ((3, 4), (4, 5))


def test_shape_three_factors():
    T1 = AlgebraicPureTensor(A, B, I)
    T2 = AlgebraicPureTensor(C, D, I)
    S = AlgebraicTensor(T1, T2)
    assert S.shape == ((3, 4), (4, 5), (5, 3))


# ---------------------------------------------------------------------------
# AlgebraicTensor: commutativity_pattern property
# ---------------------------------------------------------------------------

def test_commutativity_pattern_all_symbolic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.commutativity_pattern == (0, 0)


def test_commutativity_pattern_mixed():
    M = ImmutableDenseMatrix([[1, 2, 3, 4],
                               [5, 6, 7, 8],
                               [9, 10, 11, 12]])
    T1 = AlgebraicPureTensor(M, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.commutativity_pattern == (0, 0)


def test_commutativity_pattern_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    assert S.commutativity_pattern == (0, 0)


def test_commutativity_pattern_all_numeric():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4], [5, 6]])
    M2 = ImmutableDenseMatrix([[1, 2, 3, 4, 5],
                                [6, 7, 8, 9, 10]])
    T1 = AlgebraicPureTensor(M1, M2)
    T2 = AlgebraicPureTensor(M1, M2)
    S = AlgebraicTensor(T1, T2)
    assert S.commutativity_pattern == (1, 1)


# ---------------------------------------------------------------------------
# AlgebraicTensor: terms property
# ---------------------------------------------------------------------------

def test_terms_basic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert len(S.terms) == 2


def test_terms_excludes_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    assert len(S.terms) == 1


def test_terms_excludes_numbers():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert len(S.terms) == 2


# ---------------------------------------------------------------------------
# AlgebraicTensor: has_zero_term
# ---------------------------------------------------------------------------

def test_has_zero_term_true():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    assert S.has_zero_term() is True


def test_has_zero_term_false():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.has_zero_term() is False


def test_has_zero_term_all_cancelled():
    T = AlgebraicPureTensor(A, B)
    result = AlgebraicTensor(T, -T)
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: transpose
# ---------------------------------------------------------------------------

def test_transpose_basic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    ST = S.T
    assert ST.shape == ((4, 3), (5, 4))


def test_transpose_preserves_type():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    ST = S.T
    assert isinstance(ST, AlgebraicTensor)


def test_transpose_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    ST = S.T
    assert ST.shape == ((4, 3), (5, 4))


def test_transpose_double():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert S.T.T.shape == S.shape


def test_transpose_preserves_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    ST = S.T
    assert isinstance(ST, AlgebraicTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: negation
# ---------------------------------------------------------------------------

def test_neg_basic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    nS = -S
    assert isinstance(nS, AlgebraicTensor)


def test_neg_preserves_shape():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert (-S).shape == S.shape


def test_neg_with_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, C, D)
    S = AlgebraicTensor(T1, T2)
    nS = -S
    assert isinstance(nS, AlgebraicTensor)


def test_double_neg():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert -(-S) == S


def test_neg_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    nS = -S
    assert isinstance(nS, AlgebraicTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: addition
# ---------------------------------------------------------------------------

def test_add_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S = AlgebraicTensor(T1, T2)
    result = S + T3
    assert isinstance(result, AlgebraicTensor)


def test_add_algebraic_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    T4 = AlgebraicPureTensor(G, H)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T3, T4)
    result = S1 + S2
    assert isinstance(result, AlgebraicTensor)


def test_add_self():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S + S
    assert isinstance(result, AlgebraicTensor)


def test_radd_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S = AlgebraicTensor(T1, T2)
    result = T3 + S
    assert isinstance(result, AlgebraicTensor)


def test_add_shape_mismatch():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(MatrixSymbol("X", 2, 3),
                              MatrixSymbol("Y", 3, 4))
    S = AlgebraicTensor(T1, T2)
    raises(ShapeMismatchError, lambda: S + T3)


def test_add_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, T2)
    result = S + Z
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: subtraction
# ---------------------------------------------------------------------------

def test_sub_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S = AlgebraicTensor(T1, T2, T3)
    result = S - T1
    assert isinstance(result, AlgebraicTensor)


def test_sub_self():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S - S
    assert isinstance(result, AlgebraicZeroTensor)


def test_sub_algebraic_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T2, T3)
    result = S1 - S2
    assert isinstance(result, AlgebraicTensor)


def test_rsub_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S = AlgebraicTensor(T1, T2, T3)
    result = T1 - S
    assert isinstance(result, AlgebraicTensor)


def test_sub_cancellation():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S - T1 - T2
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: scalar multiplication
# ---------------------------------------------------------------------------

def test_mul_numeric():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S * 2
    assert isinstance(result, AlgebraicTensor)


def test_mul_symbol():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S * x
    assert isinstance(result, AlgebraicTensor)


def test_mul_S_One():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    at = AlgebraicTensor(T1, T2)
    result = at * S.One
    assert isinstance(result, AlgebraicTensor)


def test_mul_zero():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S * 0
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_rmul_numeric():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = 3 * S
    assert isinstance(result, AlgebraicTensor)


def test_rmul_symbol():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = x * S
    assert isinstance(result, AlgebraicTensor)


def test_rmul_S_One():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    at = AlgebraicTensor(T1, T2)
    result = S.One * at
    assert isinstance(result, AlgebraicTensor)


def test_rmul_zero():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = 0 * S
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: composition via __mul__
# ---------------------------------------------------------------------------

def test_mul_pure_tensor_compose():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = S * T3
    assert isinstance(result, AlgebraicTensor)


def test_rmul_pure_tensor_compose():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    at = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = at * T3
    assert isinstance(result, AlgebraicTensor)


def test_mul_algebraic_tensor_compose():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(J, L)
    T4 = AlgebraicPureTensor(E, F)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T3, T4)
    result = S1 * S2
    assert isinstance(result, AlgebraicTensor)


def test_rmul_algebraic_tensor_compose():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(J, L)
    T4 = AlgebraicPureTensor(E, F)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T3, T4)
    result = S1 * S2
    assert isinstance(result, AlgebraicTensor)


def test_mul_compose_result_shape():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = S * T3
    assert result.shape == ((3, 2), (4, 3))


def test_mul_compose_preserves_shape():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = S * T3
    assert result.shape == ((3, 2), (4, 3))


# ---------------------------------------------------------------------------
# AlgebraicTensor: _compose_with_term
# ---------------------------------------------------------------------------

def test_compose_with_term_pure_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = S._compose_with_term(T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_with_term_bare_matrix():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    at = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = at._compose_with_term(T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_with_term_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(T1, T2, Z)
    T3 = AlgebraicPureTensor(E, F)
    result = at._compose_with_term(T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_with_term_mul_coeff():
    T = AlgebraicPureTensor(A, B)
    mul_term = 2 * T
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(mul_term, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = S._compose_with_term(T3)
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# compose_algebraic_tensors function
# ---------------------------------------------------------------------------

def test_compose_both_pure_tensors():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(T1, T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((3, 2), (4, 3))


def test_compose_both_algebraic_tensors():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(J, L)
    T4 = AlgebraicPureTensor(E, F)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T3, T4)
    result = compose_algebraic_tensors(S1, S2)
    assert isinstance(result, AlgebraicTensor)


def test_compose_left_algebraic_right_pure():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(S, T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_left_pure_right_algebraic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(E, F)
    at = AlgebraicTensor(T1, T2)
    result = compose_algebraic_tensors(at, T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_left_zero_tensor():
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    T = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(Z, T)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_compose_right_zero_tensor():
    T = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((4, 2), (5, 3)))
    result = compose_algebraic_tensors(T, Z)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((4, 2), (5, 3))


def test_compose_both_zero_tensors():
    Z1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    Z2 = AlgebraicZeroTensor(((4, 2), (5, 3)))
    result = compose_algebraic_tensors(Z1, Z2)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_compose_bare_matrices():
    result = compose_algebraic_tensors(A, E)
    assert isinstance(result, MatMul)
    assert result.shape == (3, 2)


def test_compose_pure_with_bare():
    T = AlgebraicPureTensor(A)
    result = compose_algebraic_tensors(T, E)
    assert isinstance(result, MatMul)
    assert result.shape == (3, 2)


def test_compose_algebraic_with_bare():
    T1 = AlgebraicPureTensor(A)
    T2 = AlgebraicPureTensor(C)
    at = AlgebraicTensor(T1, T2)
    result = compose_algebraic_tensors(at, E)
    assert isinstance(result, AlgebraicTensor)


def test_compose_invalid_types():
    raises(TypeError, lambda: compose_algebraic_tensors(x, A))


def test_compose_algebraic_with_zero_tensor_in_args():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(T1, T2, Z)
    T3 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(at, T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_algebraic_with_zero_tensor_on_right():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    Z = AlgebraicZeroTensor(((4, 2), (5, 3)))
    result = compose_algebraic_tensors(S, Z)
    assert isinstance(result, AlgebraicZeroTensor)


def test_compose_algebraic_self():
    H1 = MatrixSymbol("H1", 3, 3)
    H2 = MatrixSymbol("H2", 3, 3)
    H3 = MatrixSymbol("H3", 3, 3)
    H4 = MatrixSymbol("H4", 3, 3)
    T1 = AlgebraicPureTensor(H1, H2)
    T2 = AlgebraicPureTensor(H3, H4)
    at = AlgebraicTensor(T1, T2)
    result = compose_algebraic_tensors(at, at)
    assert isinstance(result, AlgebraicTensor)


def test_compose_with_coefficients():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, C, D)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(S, T3)
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: simplify
# ---------------------------------------------------------------------------

def test_simplify_basic():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S.simplify()
    assert isinstance(result, AlgebraicTensor)


def test_simplify_with_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, C, D)
    S = AlgebraicTensor(T1, T2)
    result = S.simplify()
    assert isinstance(result, AlgebraicTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor: display (smoke test)
# ---------------------------------------------------------------------------

def test_display_text_mode():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    S.display(mode="text")


def test_display_latex_mode():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    S.display(mode="latex")


# ---------------------------------------------------------------------------
# AlgebraicTensor: expand
# ---------------------------------------------------------------------------

def test_expand_no_add():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S.expand()
    assert result is S


def test_expand_factor_with_matadd():
    T = AlgebraicPureTensor(A, MatAdd(B, D))
    S = AlgebraicTensor(T)
    result = S.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_multiple_terms_with_add():
    T1 = AlgebraicPureTensor(A, MatAdd(B, D))
    T2 = AlgebraicPureTensor(C, MatAdd(B, D))
    S = AlgebraicTensor(T1, T2)
    result = S.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    result = S.expand()
    assert isinstance(result, AlgebraicTensor)


def test_expand_deep_hint():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S.expand(deep=True)
    assert result is S


def test_expand_no_change_returns_self():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, C, D)
    S = AlgebraicTensor(T1, T2)
    result = S.expand()
    assert result is S


# ---------------------------------------------------------------------------
# AlgebraicTensor: args structure
# ---------------------------------------------------------------------------

def test_args_structure():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    assert len(S.args) == 2
    assert T1 in S.args
    assert T2 in S.args


def test_args_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    assert T1 in S.args
    assert Z in S.args


def test_args_with_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, C, D)
    S = AlgebraicTensor(T1, T2)
    assert len(S.args) == 2


# ---------------------------------------------------------------------------
# AlgebraicTensor: pickle
# ---------------------------------------------------------------------------

def test_pickle():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    data = pickle.dumps(S)
    S2 = pickle.loads(data)
    assert S2.shape == S.shape
    assert set(S2.args) == set(S.args)


def test_pickle_with_coeff():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(3, C, D)
    S = AlgebraicTensor(T1, T2)
    data = pickle.dumps(S)
    S2 = pickle.loads(data)
    assert S2.shape == S.shape
    assert set(S2.args) == set(S.args)


def test_pickle_with_zero_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    S = AlgebraicTensor(T1, Z)
    data = pickle.dumps(S)
    S2 = pickle.loads(data)
    assert S2 == S


# ---------------------------------------------------------------------------
# AlgebraicTensor: hash and equality
# ---------------------------------------------------------------------------

def test_hash_consistency():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T1, T2)
    assert hash(S1) == hash(S2)
    assert S1 == S2


def test_set_membership():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T1, T2)
    S3 = AlgebraicTensor(T1, T3)
    s = {S1, S2, S3}
    assert len(s) == 2


def test_equality_different_terms():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S1 = AlgebraicTensor(T1, T2)
    S2 = AlgebraicTensor(T1, T3)
    assert S1 != S2


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_constructor_sympified_coeff():
    T = AlgebraicPureTensor("2", A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T, T2)
    assert isinstance(S, AlgebraicTensor)


def test_mul_commutative_dispatch():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = x * S
    assert isinstance(result, AlgebraicTensor)


def test_mul_noncommutative_dispatch():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = S * T3
    assert isinstance(result, AlgebraicTensor)


def test_identity_preserved_through_neg():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = -S
    assert isinstance(result, AlgebraicTensor)


def test_identity_preserved_through_mul():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S * 3
    assert isinstance(result, AlgebraicTensor)


def test_identity_preserved_through_transpose():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(T1, T2)
    result = S.T
    assert isinstance(result, AlgebraicTensor)


def test_constructor_coefficient_collection_zero_result():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(-2, A, B)
    result = AlgebraicTensor(T1, T2)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_constructor_coefficient_collection_partial_zero():
    T1 = AlgebraicPureTensor(2, A, B)
    T2 = AlgebraicPureTensor(-2, A, B)
    T3 = AlgebraicPureTensor(C, D)
    result = AlgebraicTensor(T1, T2, T3)
    assert result is T3


def test_constructor_with_S_Zero_in_args():
    T = AlgebraicPureTensor(A, B)
    result = AlgebraicTensor(T, S.Zero)
    assert result is T


def test_add_dispatcher_routing():
    """Test that PureTensor + PureTensor routes through AlgebraicTensor."""
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    result = T1 + T2
    assert isinstance(result, AlgebraicTensor)


def test_add_dispatcher_algebraic_plus_pure():
    """Test that AlgebraicTensor + PureTensor routes correctly."""
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    T3 = AlgebraicPureTensor(G, H)
    S = AlgebraicTensor(T1, T2)
    result = S + T3
    assert isinstance(result, AlgebraicTensor)


def test_compose_algebraic_tensor_with_mul_in_args():
    """Composition handles Mul(coeff, PureTensor) terms in AlgebraicTensor."""
    T = AlgebraicPureTensor(A, B)
    mul_term = 2 * T
    T2 = AlgebraicPureTensor(C, D)
    S = AlgebraicTensor(mul_term, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(S, T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_reassemble_cancellation():
    """Composition that cancels out produces zero tensor."""
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(-1, A, B)
    S = AlgebraicTensor(T1, T2)
    T3 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(S, T3)
    assert isinstance(result, AlgebraicZeroTensor)


def test_compose_algebraic_tensor_zero_in_args():
    """Compose AlgebraicTensor containing zero tensor with another tensor."""
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(T1, T2, Z)
    T3 = AlgebraicPureTensor(E, F)
    result = compose_algebraic_tensors(at, T3)
    assert isinstance(result, AlgebraicTensor)


def test_compose_algebraic_tensor_both_with_zero():
    """Compose two AlgebraicTensors both containing zero tensors."""
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    Z1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
    T3 = AlgebraicPureTensor(J, L)
    T4 = AlgebraicPureTensor(E, F)
    Z2 = AlgebraicZeroTensor(((4, 2), (5, 3)))
    S1 = AlgebraicTensor(T1, T2, Z1)
    S2 = AlgebraicTensor(T3, T4, Z2)
    result = compose_algebraic_tensors(S1, S2)
    assert isinstance(result, AlgebraicTensor)
