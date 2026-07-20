from __future__ import annotations

from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.matrices.expressions import MatrixSymbol, MatAdd
from sympy.matrices.immutable import ImmutableDenseMatrix

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    tensorsimplify,
)
from sympy.tensor.algebraic.simplify import (
    _build_pt,
    _deduplicate_proportional,
    _decompose_commutative_factors,
    _extract_commutative_from_factor,
    _extract_commutative_prefactors,
    _extract_pt_and_coeff,
    _is_exactly_divisible,
    _matrix_proportionality_ratio,
    _normalize_factor_sign,
    _proportionality_ratio,
    _reconstruct_term,
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

x = Symbol("x")
y = Symbol("y")
w = Symbol("w")
z = Symbol("z")


# ---------------------------------------------------------------------------
# _matrix_proportionality_ratio
# ---------------------------------------------------------------------------

def test_matrix_proportionality_ratio_identical():
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    assert _matrix_proportionality_ratio(M, M) == S.One


def test_matrix_proportionality_ratio_scaled():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    assert _matrix_proportionality_ratio(M1, M2) == S.Half


def test_matrix_proportionality_ratio_scaled_reverse():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    assert _matrix_proportionality_ratio(M2, M1) == 2


def test_matrix_proportionality_ratio_negative_scale():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[-1, -2], [-3, -4]])
    assert _matrix_proportionality_ratio(M1, M2) == -1


def test_matrix_proportionality_ratio_not_proportional():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[1, 2], [3, 5]])
    assert _matrix_proportionality_ratio(M1, M2) is None


def test_matrix_proportionality_ratio_shape_mismatch():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[1, 2, 3], [4, 5, 6]])
    assert _matrix_proportionality_ratio(M1, M2) is None


def test_matrix_proportionality_ratio_with_zeros():
    M1 = ImmutableDenseMatrix([[1, 0], [0, 2]])
    M2 = ImmutableDenseMatrix([[2, 0], [0, 4]])
    assert _matrix_proportionality_ratio(M1, M2) == S.Half


def test_matrix_proportionality_ratio_zero_mismatch():
    M1 = ImmutableDenseMatrix([[1, 0], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 1], [6, 8]])
    assert _matrix_proportionality_ratio(M1, M2) is None


def test_matrix_proportionality_ratio_all_zero():
    M1 = ImmutableDenseMatrix([[0, 0], [0, 0]])
    M2 = ImmutableDenseMatrix([[0, 0], [0, 0]])
    assert _matrix_proportionality_ratio(M1, M2) is None


def test_matrix_proportionality_ratio_symbolic():
    M1 = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    M2 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    assert _matrix_proportionality_ratio(M1, M2) == x


# ---------------------------------------------------------------------------
# _proportionality_ratio
# ---------------------------------------------------------------------------

def test_proportionality_ratio_identical():
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    assert _proportionality_ratio(M, M) == S.One


def test_proportionality_ratio_scaled():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    assert _proportionality_ratio(M1, M2) == S.Half


def test_proportionality_ratio_not_proportional():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[1, 2], [3, 5]])
    assert _proportionality_ratio(M1, M2) is None


# ---------------------------------------------------------------------------
# _extract_pt_and_coeff
# ---------------------------------------------------------------------------

def test_extract_pt_and_coeff_pure_tensor_no_coeff():
    T = AlgebraicPureTensor(A, B)
    coeff, factors = _extract_pt_and_coeff(T)
    assert coeff == S.One
    assert factors == [A, B]


def test_extract_pt_and_coeff_pure_tensor_with_coeff():
    T = AlgebraicPureTensor(2, A, B)
    coeff, factors = _extract_pt_and_coeff(T)
    assert coeff == 2
    assert factors == [A, B]


def test_extract_pt_and_coeff_pure_tensor_symbolic_coeff():
    T = AlgebraicPureTensor(x, A, B)
    coeff, factors = _extract_pt_and_coeff(T)
    assert coeff == x
    assert factors == [A, B]


def test_extract_pt_and_coeff_mul_with_pure_tensor():
    T = AlgebraicPureTensor(2, A, B)
    mul_expr = x * T
    # x * AlgebraicPureTensor(2, A, B) -> AlgebraicPureTensor(2*x, A, B)
    assert isinstance(mul_expr, AlgebraicPureTensor)
    coeff, factors = _extract_pt_and_coeff(mul_expr)
    assert coeff == 2 * x
    assert factors == [A, B]


def test_extract_pt_and_coeff_bare_matrix():
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    coeff, factors = _extract_pt_and_coeff(M)
    assert coeff == S.One
    assert factors == [M]


def test_extract_pt_and_coeff_mul_with_matrices():
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    mul_expr = x * M
    # x * ImmutableDenseMatrix -> new ImmutableDenseMatrix with scaled entries
    assert isinstance(mul_expr, ImmutableDenseMatrix)
    coeff, factors = _extract_pt_and_coeff(mul_expr)
    assert coeff == S.One
    assert factors == [mul_expr]


def test_extract_pt_and_coeff_unknown():
    coeff, factors = _extract_pt_and_coeff(x)
    assert coeff == S.One
    assert factors == [x]


# ---------------------------------------------------------------------------
# _build_pt
# ---------------------------------------------------------------------------

def test_build_pt_zero_coeff():
    result = _build_pt(S.Zero, [A, B])
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_build_pt_zero_coeff_empty_factors():
    result = _build_pt(S.Zero, [])
    assert result == S.Zero


def test_build_pt_empty_factors():
    result = _build_pt(2, [])
    assert result == 2


def test_build_pt_coeff_one_single_factor():
    result = _build_pt(S.One, [A])
    assert result is A


def test_build_pt_coeff_one_multi_factor():
    result = _build_pt(S.One, [A, B])
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, B)


def test_build_pt_numeric_coeff():
    result = _build_pt(2, [A, B])
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2
    assert result.factors == (A, B)


def test_build_pt_symbolic_coeff():
    result = _build_pt(x, [A, B])
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x
    assert result.factors == (A, B)


def test_build_pt_negative_coeff():
    result = _build_pt(-3, [A, B])
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == -3
    assert result.factors == (A, B)


# ---------------------------------------------------------------------------
# _normalize_factor_sign
# ---------------------------------------------------------------------------

def test_normalize_factor_sign_already_positive():
    assert _normalize_factor_sign(x - y) == x - y


def test_normalize_factor_sign_negative_leading():
    assert _normalize_factor_sign(-x + y) == x - y


def test_normalize_factor_sign_numeric_mul():
    result = _normalize_factor_sign(-2 * x)
    assert result == 2 * x


def test_normalize_factor_sign_no_change():
    assert _normalize_factor_sign(x + y) == x + y


def test_normalize_factor_sign_symbol():
    assert _normalize_factor_sign(x) == x


# ---------------------------------------------------------------------------
# _is_exactly_divisible
# ---------------------------------------------------------------------------

def test_is_exactly_divisible_true():
    assert _is_exactly_divisible(6, 2) is True


def test_is_exactly_divisible_false():
    assert _is_exactly_divisible(5, 2) is False


def test_is_exactly_divisible_zero():
    assert _is_exactly_divisible(0, 3) is True


def test_is_exactly_divisible_symbolic_true():
    assert _is_exactly_divisible(x * y, x) is True


def test_is_exactly_divisible_symbolic_false():
    assert _is_exactly_divisible(x + 1, x) is False


def test_is_exactly_divisible_negative():
    assert _is_exactly_divisible(-6, 3) is True


def test_is_exactly_divisible_float_false():
    assert _is_exactly_divisible(5, 2) is False


# ---------------------------------------------------------------------------
# _deduplicate_proportional
# ---------------------------------------------------------------------------

def test_deduplicate_proportional_basic():
    result = _deduplicate_proportional([x, -x, x + 1])
    assert len(result) == 2


def test_deduplicate_proportional_empty():
    assert _deduplicate_proportional([]) == []


def test_deduplicate_proportional_no_duplicates():
    result = _deduplicate_proportional([x, y, x + y])
    assert len(result) == 3


def test_deduplicate_proportional_all_same():
    result = _deduplicate_proportional([x, -x, 2*x])
    assert len(result) == 1


def test_deduplicate_proportional_with_zero():
    result = _deduplicate_proportional([x, S.Zero, -x])
    assert len(result) == 2


# ---------------------------------------------------------------------------
# _extract_commutative_from_factor
# ---------------------------------------------------------------------------

def test_extract_commutative_from_factor_no_extraction():
    coeff, new_factor = _extract_commutative_from_factor(A)
    assert coeff == S.One
    assert new_factor is A


def test_extract_commutative_from_factor_matrix_common():
    M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    coeff, new_M = _extract_commutative_from_factor(M)
    assert coeff == x
    expected = ImmutableDenseMatrix([[1, 2], [3, 4]])
    assert new_M == expected


def test_extract_commutative_from_factor_matrix_numeric_gcd():
    M = ImmutableDenseMatrix([[2, 4], [6, 8]])
    coeff, new_M = _extract_commutative_from_factor(M)
    assert coeff == 2
    expected = ImmutableDenseMatrix([[1, 2], [3, 4]])
    assert new_M == expected


def test_extract_commutative_from_factor_matrix_no_common():
    M = ImmutableDenseMatrix([[1, 2], [3, 5]])
    coeff, new_M = _extract_commutative_from_factor(M)
    assert coeff == S.One
    assert new_M == M


def test_extract_commutative_from_factor_matmul():
    expr = x * A  # This is a MatMul, not a Mul
    coeff, new_factor = _extract_commutative_from_factor(expr)
    assert coeff == x
    assert new_factor == A


def test_extract_commutative_from_factor_add():
    expr = x * A + x * C
    coeff, new_factor = _extract_commutative_from_factor(expr)
    assert coeff == x


def test_extract_commutative_from_factor_matadd():
    expr = MatAdd(x * A, x * C)
    coeff, new_factor = _extract_commutative_from_factor(expr)
    assert coeff == x


def test_extract_commutative_from_factor_matadd_partial():
    expr = MatAdd(x * A, C)
    coeff, new_factor = _extract_commutative_from_factor(expr)
    assert coeff == S.One


def test_extract_commutative_from_factor_zero_matrix():
    M = ImmutableDenseMatrix.zeros(2, 2)
    coeff, new_M = _extract_commutative_from_factor(M)
    assert coeff == S.One


def test_extract_commutative_from_factor_symbolic_matrix():
    M = ImmutableDenseMatrix([[x*y, 2*x*y], [3*x*y, 4*x*y]])
    coeff, new_M = _extract_commutative_from_factor(M)
    assert coeff == x * y


# ---------------------------------------------------------------------------
# _extract_commutative_prefactors
# ---------------------------------------------------------------------------

def test_extract_commutative_prefactors_no_extraction():
    T = AlgebraicPureTensor(A, B)
    coeff, factors = _extract_commutative_prefactors(T)
    assert coeff == S.One
    assert factors == [A, B]


def test_extract_commutative_prefactors_with_matrix_factor():
    M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    T = AlgebraicPureTensor(M, B)
    coeff, factors = _extract_commutative_prefactors(T)
    assert coeff == x


def test_extract_commutative_prefactors_multiple_factors():
    M1 = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    M2 = ImmutableDenseMatrix([[y, 2*y], [3*y, 4*y]])
    T = AlgebraicPureTensor(M1, M2)
    coeff, factors = _extract_commutative_prefactors(T)
    assert coeff == x * y


def test_extract_commutative_prefactors_bare_matrix():
    M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    coeff, factors = _extract_commutative_prefactors(M)
    assert coeff == x


# ---------------------------------------------------------------------------
# _decompose_commutative_factors
# ---------------------------------------------------------------------------

def test_decompose_commutative_factors_empty():
    result = _decompose_commutative_factors([], S.One)
    assert result == [(S.One, [])]


def test_decompose_commutative_factors_single_matrix():
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    result = _decompose_commutative_factors([M], S.One)
    assert len(result) == 4


def test_decompose_commutative_factors_with_coeff():
    M = ImmutableDenseMatrix([[1, 0], [0, 0]])
    result = _decompose_commutative_factors([M], x)
    assert len(result) == 1
    assert result[0][0] == x


def test_decompose_commutative_factors_zero_entries_skipped():
    M = ImmutableDenseMatrix([[1, 0], [0, 0]])
    result = _decompose_commutative_factors([M], S.One)
    assert len(result) == 1


def test_decompose_commutative_factors_basis_matrices():
    M = ImmutableDenseMatrix([[1, 0], [0, 0]])
    result = _decompose_commutative_factors([M], S.One)
    assert len(result) == 1
    basis = result[0][1][0]
    assert basis == M


# ---------------------------------------------------------------------------
# _reconstruct_term
# ---------------------------------------------------------------------------

def test_reconstruct_term_all_commutative():
    E00 = ImmutableDenseMatrix([[1, 0], [0, 0]])
    result = _reconstruct_term((E00,), None, S.One, (1,), [0], [])
    assert result == E00


def test_reconstruct_term_zero_coeff():
    E00 = ImmutableDenseMatrix([[1, 0], [0, 0]])
    result = _reconstruct_term((E00,), None, S.Zero, (1,), [0], [])
    assert result is None


def test_reconstruct_term_with_coeff():
    E00 = ImmutableDenseMatrix([[1, 0], [0, 0]])
    result = _reconstruct_term((E00,), None, 2, (1,), [0], [])
    # Single factor with coeff -> _build_pt returns ImmutableDenseMatrix
    assert isinstance(result, ImmutableDenseMatrix)
    assert result == ImmutableDenseMatrix([[2, 0], [0, 0]])


# ---------------------------------------------------------------------------
# _simplify_algebraic_pure_tensor
# ---------------------------------------------------------------------------

def test_simplify_pure_tensor_basic():
    T = AlgebraicPureTensor(2, A, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2


def test_simplify_pure_tensor_coeff_simplification():
    T = AlgebraicPureTensor(x + x, A, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2 * x


def test_simplify_pure_tensor_zero_coeff():
    T = AlgebraicPureTensor(S.Zero, A, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicZeroTensor)


def test_simplify_pure_tensor_factor_extraction():
    M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    T = AlgebraicPureTensor(M, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x


def test_simplify_pure_tensor_no_change():
    T = AlgebraicPureTensor(A, B)
    result = tensorsimplify(T)
    assert result is T


# ---------------------------------------------------------------------------
# tensorsimplify: public API
# ---------------------------------------------------------------------------

def test_tensorsimplify_zero_tensor():
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = tensorsimplify(Z)
    assert result is Z


def test_tensorsimplify_pure_tensor():
    T = AlgebraicPureTensor(2, A, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)


def test_tensorsimplify_algebraic_tensor():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(2, A, B)
    at = AlgebraicTensor(T1, T2)
    result = tensorsimplify(at)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 3


def test_tensorsimplify_non_tensor():
    result = tensorsimplify(x + x)
    assert result == 2 * x


def test_tensorsimplify_proportional_merge():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(2, A, B)
    result = tensorsimplify(T1 + T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 3


def test_tensorsimplify_cancel():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(-1, A, B)
    result = tensorsimplify(T1 + T2)
    assert isinstance(result, AlgebraicZeroTensor)


def test_tensorsimplify_no_merge():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, D)
    result = tensorsimplify(T1 + T2)
    assert isinstance(result, AlgebraicTensor)


def test_tensorsimplify_with_kwargs():
    T = AlgebraicPureTensor(2, A, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)


def test_tensorsimplify_complex_coeff():
    T1 = AlgebraicPureTensor(x + y, A, B)
    T2 = AlgebraicPureTensor(x - y, A, B)
    result = tensorsimplify(T1 + T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2 * x


def test_tensorsimplify_numeric_matrix_proportional():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[2, 4], [6, 8]])
    T1 = AlgebraicPureTensor(M1, B)
    T2 = AlgebraicPureTensor(M2, B)
    result = tensorsimplify(T1 + T2)
    assert result is not None


def test_tensorsimplify_symbolic_matrix_factor_extraction():
    M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    T = AlgebraicPureTensor(M, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x


def test_tensorsimplify_matmul_term():
    M = ImmutableDenseMatrix([[1, 2], [3, 4]])
    T = AlgebraicPureTensor(M, B)
    result = tensorsimplify(T)
    assert result is not None


def test_tensorsimplify_three_factor_tensor():
    G = MatrixSymbol("G", 5, 3)
    T1 = AlgebraicPureTensor(A, B, G)
    T2 = AlgebraicPureTensor(2, A, B, G)
    result = tensorsimplify(T1 + T2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 3


def test_tensorsimplify_mixed_proportional():
    M = ImmutableDenseMatrix([[1, 2, 3, 4],
                               [5, 6, 7, 8],
                               [9, 10, 11, 12]])
    T1 = AlgebraicPureTensor(M, B)
    T2 = AlgebraicPureTensor(A, B)
    result = tensorsimplify(T1 + T2)
    assert result is not None


def test_tensorsimplify_linear_combination_slot():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(C, B)
    result = tensorsimplify(T1 + T2)
    assert result is not None


def test_tensorsimplify_all_terms_cancel():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(-1, A, B)
    T3 = AlgebraicPureTensor(C, D)
    T4 = AlgebraicPureTensor(-1, C, D)
    result = tensorsimplify(T1 + T2 + T3 + T4)
    assert isinstance(result, AlgebraicZeroTensor)


def test_tensorsimplify_preserves_shape():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(-1, A, B)
    result = tensorsimplify(T1 + T2)
    assert result.shape == ((3, 4), (4, 5))


def test_tensorsimplify_single_term_algebraic_tensor():
    T1 = AlgebraicPureTensor(A, B)
    Z = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = tensorsimplify(T1 + Z)
    assert result is T1


def test_tensorsimplify_commutative_matrix_factor():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[5, 6], [7, 8]])
    T = AlgebraicPureTensor(M1, M2)
    result = tensorsimplify(T)
    assert result is not None


def test_tensorsimplify_with_matadd_factor():
    T = AlgebraicPureTensor(MatAdd(A, C), B)
    result = tensorsimplify(T)
    assert result is not None


def test_tensorsimplify_numeric_gcd_extraction():
    M = ImmutableDenseMatrix([[2, 4], [6, 8]])
    T = AlgebraicPureTensor(M, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2


def test_tensorsimplify_multiple_factor_extraction():
    M1 = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    M2 = ImmutableDenseMatrix([[y, 2*y], [3*y, 4*y]])
    T = AlgebraicPureTensor(M1, M2)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x * y


def test_tensorsimplify_symbolic_coeff_with_extraction():
    M = ImmutableDenseMatrix([[x, 2*x], [3*x, 4*x]])
    T = AlgebraicPureTensor(2, M, B)
    result = tensorsimplify(T)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2 * x


def test_tensorsimplify_matadd_with_common_factor():
    expr = MatAdd(x * A, x * C)
    T = AlgebraicPureTensor(expr, B)
    result = tensorsimplify(T)
    assert result is not None


def test_tensorsimplify_negative_proportional():
    M1 = ImmutableDenseMatrix([[1, 2], [3, 4]])
    M2 = ImmutableDenseMatrix([[-1, -2], [-3, -4]])
    T1 = AlgebraicPureTensor(M1, B)
    T2 = AlgebraicPureTensor(M2, B)
    result = tensorsimplify(T1 + T2)
    assert isinstance(result, AlgebraicZeroTensor)


def test_tensorsimplify_preserves_non_tensor_terms():
    T1 = AlgebraicPureTensor(A, B)
    T2 = AlgebraicPureTensor(2, A, B)
    result = tensorsimplify(T1 + T2)
    assert isinstance(result, AlgebraicPureTensor)
