from __future__ import annotations

from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.core.symbol import Symbol
from sympy.matrices.expressions import MatrixSymbol, MatMul
from sympy.matrices.immutable import ImmutableDenseMatrix
from sympy.testing.pytest import raises

from sympy.tensor.algebraic import (
    AlgebraicPureTensor,
    AlgebraicTensor,
    AlgebraicZeroTensor,
    algebraic_tensor_product,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 4, 5)
C = MatrixSymbol("C", 3, 4)
D = MatrixSymbol("D", 4, 5)
E = MatrixSymbol("E", 5, 3)
F = MatrixSymbol("F", 5, 3)
G = MatrixSymbol("G", 3, 4)
H = MatrixSymbol("H", 4, 5)
I3 = MatrixSymbol("I", 3, 3)
J = MatrixSymbol("J", 2, 3)
K = MatrixSymbol("K", 3, 2)
# Same-shape matrices for AlgebraicTensor terms
P = MatrixSymbol("P", 5, 3)
Q = MatrixSymbol("Q", 5, 3)

x = Symbol("x")
y = Symbol("y")


# ---------------------------------------------------------------------------
# No arguments
# ---------------------------------------------------------------------------

def test_no_args():
    raises(ValueError, lambda: algebraic_tensor_product())


# ---------------------------------------------------------------------------
# Single argument
# ---------------------------------------------------------------------------

def test_single_matrix():
    result = algebraic_tensor_product(A)
    assert result is A


def test_single_pure_tensor():
    pt = AlgebraicPureTensor(A, B)
    result = algebraic_tensor_product(pt)
    assert result is pt


def test_single_algebraic_tensor():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(at)
    assert result is at


def test_single_zero_tensor():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = algebraic_tensor_product(zt)
    assert result is zt


def test_single_scalar():
    result = algebraic_tensor_product(2)
    assert result == 2


def test_single_symbol():
    result = algebraic_tensor_product(x)
    assert result is x


# ---------------------------------------------------------------------------
# Two arguments: scalar first (special fast path)
# ---------------------------------------------------------------------------

def test_two_args_number_first():
    result = algebraic_tensor_product(2, A)
    assert result == 2 * A


def test_two_args_number_first_with_pure_tensor():
    pt = AlgebraicPureTensor(A, B)
    result = algebraic_tensor_product(2, pt)
    assert result == 2 * pt


def test_two_args_number_first_with_algebraic_tensor():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(2, at)
    assert result == 2 * at


def test_two_args_number_first_with_zero_tensor():
    """0 * zt returns zt (AlgebraicZeroTensor.__mul__ returns self for commutative)."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = algebraic_tensor_product(0, zt)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


def test_two_args_symbol_first():
    result = algebraic_tensor_product(x, A)
    assert result == x * A


def test_two_args_symbol_first_with_pure_tensor():
    pt = AlgebraicPureTensor(A, B)
    result = algebraic_tensor_product(x, pt)
    assert result == x * pt


def test_two_args_symbol_first_with_algebraic_tensor():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(x, at)
    assert result == x * at


def test_two_args_S_One_first():
    result = algebraic_tensor_product(S.One, A)
    assert result == A


def test_two_args_S_Zero_first():
    """S.Zero * pt returns AlgebraicZeroTensor (zero * tensor = zero tensor)."""
    pt = AlgebraicPureTensor(A, B)
    result = algebraic_tensor_product(S.Zero, pt)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


# ---------------------------------------------------------------------------
# Two arguments: basic matrix product
# ---------------------------------------------------------------------------

def test_two_matrices():
    result = algebraic_tensor_product(A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((3, 4), (4, 5))
    assert result.factors == (A, B)
    assert result.coeff == S.One


def test_two_matrices_single_factor_unwrap():
    """When result is a single factor with coeff 1, it unwraps."""
    result = algebraic_tensor_product(A)
    assert result is A


def test_three_matrices():
    result = algebraic_tensor_product(A, B, E)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((3, 4), (4, 5), (5, 3))
    assert result.factors == (A, B, E)


# ---------------------------------------------------------------------------
# Scalar coefficient in product
# ---------------------------------------------------------------------------

def test_scalar_then_two_matrices():
    result = algebraic_tensor_product(2, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2
    assert result.factors == (A, B)


def test_symbolic_coeff_then_matrices():
    result = algebraic_tensor_product(x, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x
    assert result.factors == (A, B)


def test_two_symbolic_coeffs():
    result = algebraic_tensor_product(x, y, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x * y
    assert result.factors == (A, B)


def test_coeff_in_middle():
    result = algebraic_tensor_product(A, x, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x
    assert result.factors == (A, B)


def test_coeff_at_end():
    result = algebraic_tensor_product(A, B, x)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x
    assert result.factors == (A, B)


def test_multiple_coeffs_scattered():
    result = algebraic_tensor_product(x, A, y, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x * y
    assert result.factors == (A, B)


def test_negative_coeff():
    result = algebraic_tensor_product(-1, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == -1
    assert result.factors == (A, B)


def test_zero_coeff():
    result = algebraic_tensor_product(0, A, B)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


# ---------------------------------------------------------------------------
# PureTensor argument
# ---------------------------------------------------------------------------

def test_pure_tensor_with_matrix():
    pt = AlgebraicPureTensor(A, B)
    result = algebraic_tensor_product(pt, E)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((3, 4), (4, 5), (5, 3))
    assert result.factors == (A, B, E)


def test_matrix_with_pure_tensor():
    pt = AlgebraicPureTensor(B, E)
    result = algebraic_tensor_product(A, pt)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((3, 4), (4, 5), (5, 3))
    assert result.factors == (A, B, E)


def test_pure_tensor_with_coeff_flattens():
    pt = AlgebraicPureTensor(y, B, A.T)
    result = algebraic_tensor_product(x, A, pt)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x * y
    assert result.factors == (A, B, A.T)


def test_pure_tensor_coeff_accumulates():
    pt1 = AlgebraicPureTensor(2, A, B)
    pt2 = AlgebraicPureTensor(3, C, D)
    result = algebraic_tensor_product(pt1, pt2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6
    assert result.factors == (A, B, C, D)


def test_two_pure_tensors_no_coeff():
    pt1 = AlgebraicPureTensor(A, B)
    pt2 = AlgebraicPureTensor(C, D)
    result = algebraic_tensor_product(pt1, pt2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == S.One
    assert result.factors == (A, B, C, D)


# ---------------------------------------------------------------------------
# AlgebraicTensor argument (distributes over sum)
# ---------------------------------------------------------------------------

def test_algebraic_tensor_with_matrix():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(at, E)
    assert isinstance(result, AlgebraicTensor)
    assert result.shape == ((3, 4), (4, 5), (5, 3))


def test_matrix_with_algebraic_tensor():
    at = AlgebraicTensor(AlgebraicPureTensor(B, E), AlgebraicPureTensor(D, E))
    result = algebraic_tensor_product(A, at)
    assert isinstance(result, AlgebraicTensor)


def test_two_algebraic_tensors():
    at1 = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    at2 = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(at1, at2)
    assert isinstance(result, AlgebraicTensor)
    # Should have 4 terms: (A,B,E), (A,B,F), (C,D,E), (C,D,F)
    assert len(result.terms) == 4


def test_algebraic_tensor_coeff_distributes():
    at = AlgebraicTensor(
        AlgebraicPureTensor(2, A, B),
        AlgebraicPureTensor(3, C, D)
    )
    result = algebraic_tensor_product(at, E)
    assert isinstance(result, AlgebraicTensor)
    # Check that coefficients are preserved
    found_2 = False
    found_3 = False
    for term in result.terms:
        if isinstance(term, AlgebraicPureTensor):
            if term.coeff == 2:
                found_2 = True
            if term.coeff == 3:
                found_3 = True
    assert found_2
    assert found_3


def test_algebraic_tensor_with_scalar():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(2, at)
    assert result == 2 * at


def test_algebraic_tensor_with_symbolic_coeff():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(x, at)
    assert result == x * at


# ---------------------------------------------------------------------------
# Zero tensor argument
# ---------------------------------------------------------------------------

def test_zero_tensor_with_matrix():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = algebraic_tensor_product(A, zt)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (3, 4), (4, 5))


def test_matrix_with_zero_tensor():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = algebraic_tensor_product(zt, E)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5), (5, 3))


def test_zero_tensor_with_pure_tensor():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(E)
    result = algebraic_tensor_product(zt, pt)
    assert isinstance(result, AlgebraicZeroTensor)


def test_zero_tensor_with_algebraic_tensor():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(zt, at)
    assert isinstance(result, AlgebraicZeroTensor)


def test_two_zero_tensors():
        zt1 = AlgebraicZeroTensor(((3, 4), (4, 5)))
        zt2 = AlgebraicZeroTensor(((5, 3),))
        result = algebraic_tensor_product(zt1, zt2)
        assert isinstance(result, AlgebraicZeroTensor)
        assert result.shape == ((3, 4), (4, 5), (5, 3))


def test_zero_tensor_with_scalar():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = algebraic_tensor_product(2, zt)
    # Two-arg fast path: 2 * zt, which returns zt (commutative mul returns self)
    assert isinstance(result, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# AlgebraicTensor with only zero terms
# ---------------------------------------------------------------------------

def test_algebraic_tensor_only_zero_terms():
    """AlgebraicTensor containing only a zero tensor produces zero result."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), -AlgebraicPureTensor(A, B), zt)
    # at should be a zero tensor since terms cancel
    assert isinstance(at, AlgebraicZeroTensor)


def test_algebraic_tensor_empty_terms_with_zero():
    """When AlgebraicTensor has only zero terms, product is zero."""
    # Build an AlgebraicTensor that internally has no real terms
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(-1, A, B))
    # This should be a zero tensor
    assert isinstance(at, AlgebraicZeroTensor)


# ---------------------------------------------------------------------------
# Single factor result unwrap
# ---------------------------------------------------------------------------

def test_single_factor_unwrap():
    """Single factor with coeff 1 unwraps to bare matrix."""
    result = algebraic_tensor_product(A)
    assert result is A


def test_single_factor_with_coeff():
    """Single factor with coeff != 1 stays as PureTensor or Mul."""
    result = algebraic_tensor_product(2, A)
    # Two-arg fast path: 2 * A (plain Mul)
    assert result == 2 * A


def test_single_factor_with_symbolic_coeff():
    result = algebraic_tensor_product(x, A)
    # Two-arg fast path: x * A
    assert result == x * A


# ---------------------------------------------------------------------------
# All scalar arguments (no tensor factors)
# ---------------------------------------------------------------------------

def test_all_numbers():
    """Two numbers use the two-arg fast path: first * sympify(second)."""
    result = algebraic_tensor_product(2, 3)
    assert result == 6


def test_all_symbols():
    """Two symbols use the two-arg fast path: first_s * sympify(second)."""
    result = algebraic_tensor_product(x, y)
    assert result == x * y


def test_single_scalar_then_single_matrix():
    """Scalar + single matrix uses two-arg fast path."""
    result = algebraic_tensor_product(2, A)
    assert result == 2 * A


# ---------------------------------------------------------------------------
# TypeError for unrecognized types
# ---------------------------------------------------------------------------

def test_none_argument():
    raises(TypeError, lambda: algebraic_tensor_product(None, A))


# ---------------------------------------------------------------------------
# String arguments are sympified to symbols
# ---------------------------------------------------------------------------

def test_string_sympified_to_symbol():
    """String "x" is sympified to Symbol('x'), treated as commutative coeff."""
    result = algebraic_tensor_product("x", A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == Symbol("x")
    assert result.factors == (A, B)


def test_string_in_middle_becomes_coeff():
    """String in the middle is sympified to a symbol (coefficient)."""
    result = algebraic_tensor_product(A, "y", B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == Symbol("y")
    assert result.factors == (A, B)


# ---------------------------------------------------------------------------
# Combined shape computation
# ---------------------------------------------------------------------------

def test_combined_shape_two_factors():
    result = algebraic_tensor_product(A, B)
    assert result.shape == ((3, 4), (4, 5))


def test_combined_shape_three_factors():
    result = algebraic_tensor_product(A, B, E)
    assert result.shape == ((3, 4), (4, 5), (5, 3))


def test_combined_shape_with_pure_tensor():
    pt = AlgebraicPureTensor(B, E)
    result = algebraic_tensor_product(A, pt)
    assert result.shape == ((3, 4), (4, 5), (5, 3))


def test_combined_shape_with_algebraic_tensor():
    at = AlgebraicTensor(AlgebraicPureTensor(B, E), AlgebraicPureTensor(D, E))
    result = algebraic_tensor_product(A, at)
    assert result.shape == ((3, 4), (4, 5), (5, 3))


def test_combined_shape_zero_tensor():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = algebraic_tensor_product(A, zt)
    assert result.shape == ((3, 4), (3, 4), (4, 5))


# ---------------------------------------------------------------------------
# Multiple AlgebraicTensor arguments
# ---------------------------------------------------------------------------

def test_three_args_two_algebraic_tensors():
    at1 = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    at2 = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(at1, at2)
    assert isinstance(result, AlgebraicTensor)
    assert len(result.terms) == 4


def test_three_algebraic_tensors():
    at1 = AlgebraicTensor(AlgebraicPureTensor(A), AlgebraicPureTensor(C))
    at2 = AlgebraicTensor(AlgebraicPureTensor(B), AlgebraicPureTensor(D))
    at3 = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(at1, at2, at3)
    assert isinstance(result, AlgebraicTensor)
    # 2 * 2 * 2 = 8 terms
    assert len(result.terms) == 8


def test_algebraic_tensor_with_scalar_and_matrix():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(x, at, E)
    assert isinstance(result, AlgebraicTensor)
    # Both terms should have coefficient x
    for term in result.terms:
        if isinstance(term, AlgebraicPureTensor):
            assert term.coeff == x


# ---------------------------------------------------------------------------
# Numeric matrix arguments
# ---------------------------------------------------------------------------

def test_numeric_matrix():
    M = ImmutableDenseMatrix([[1, 2], [3, 4], [5, 6]])
    N = ImmutableDenseMatrix([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
    result = algebraic_tensor_product(M, N)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((3, 2), (2, 5))


def test_numeric_with_symbolic():
    M = ImmutableDenseMatrix([[1, 2], [3, 4], [5, 6]])
    result = algebraic_tensor_product(M, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((3, 2), (4, 5))


# ---------------------------------------------------------------------------
# Sympify behavior
# ---------------------------------------------------------------------------

def test_string_number_sympified():
    result = algebraic_tensor_product("2", A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 2


# ---------------------------------------------------------------------------
# Edge cases with MatMul factors
# ---------------------------------------------------------------------------

def test_matmul_factor():
    """MatMul expressions are accepted as matrix-like factors."""
    mm = MatMul(A, B, evaluate=False)
    result = algebraic_tensor_product(mm, E)
    assert isinstance(result, AlgebraicPureTensor)
    assert mm in result.factors


# ---------------------------------------------------------------------------
# Coefficient multiplication across PureTensor args
# ---------------------------------------------------------------------------

def test_coeff_multiply_across_pure_tensors():
    pt1 = AlgebraicPureTensor(2, A, B)
    pt2 = AlgebraicPureTensor(3, C, D)
    result = algebraic_tensor_product(pt1, pt2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6
    assert result.factors == (A, B, C, D)


def test_symbolic_coeff_multiply():
    pt1 = AlgebraicPureTensor(x, A, B)
    pt2 = AlgebraicPureTensor(y, C, D)
    result = algebraic_tensor_product(pt1, pt2)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x * y
    assert result.factors == (A, B, C, D)


def test_coeff_with_scalar_arg():
    pt = AlgebraicPureTensor(2, A, B)
    result = algebraic_tensor_product(3, pt)
    # Two-arg fast path: 3 * pt
    assert result == 3 * pt


def test_coeff_with_scalar_arg_long_form():
    pt = AlgebraicPureTensor(2, A, B)
    result = algebraic_tensor_product(3, pt, E)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6
    assert result.factors == (A, B, E)


# ---------------------------------------------------------------------------
# AlgebraicTensor with zero tensor in args
# ---------------------------------------------------------------------------

def test_algebraic_tensor_with_zero_in_args():
    """AlgebraicTensor containing AlgebraicZeroTensor in args.
    Zero terms are skipped in _terms_of_arg, so only real terms appear.
    With only 1 real term, the result unwraps to PureTensor."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), zt)
    result = algebraic_tensor_product(at, E)
    # Only 1 real term, so result unwraps to PureTensor
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, B, E)


# ---------------------------------------------------------------------------
# Commutative symbol treated as coefficient
# ---------------------------------------------------------------------------

def test_commutative_symbol_as_coeff():
    """Commutative symbols are treated as coefficients, not factors."""
    result = algebraic_tensor_product(x, A)
    # Two-arg fast path: x * A
    assert result == x * A


def test_commutative_expr_two_args():
    """With 2 args, commutative expr uses fast path."""
    result = algebraic_tensor_product(x + y, A)
    assert result == (x + y) * A


def test_commutative_expr_three_args():
    """With 3+ args, commutative expr becomes coefficient in product."""
    result = algebraic_tensor_product(x + y, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == x + y
    assert result.factors == (A, B)


# ---------------------------------------------------------------------------
# Result type checks
# ---------------------------------------------------------------------------

def test_result_is_pure_tensor():
    result = algebraic_tensor_product(A, B)
    assert isinstance(result, AlgebraicPureTensor)


def test_result_is_algebraic_tensor():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(at, E)
    assert isinstance(result, AlgebraicTensor)


def test_result_is_zero_tensor():
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    result = algebraic_tensor_product(zt, A)
    assert isinstance(result, AlgebraicZeroTensor)


def test_result_is_bare_matrix():
    result = algebraic_tensor_product(A)
    assert result is A


def test_result_is_mul():
    result = algebraic_tensor_product(2, A)
    # Two-arg fast path returns 2 * A which is a Mul
    assert isinstance(result, Mul)


# ---------------------------------------------------------------------------
# Complex combinations
# ---------------------------------------------------------------------------

def test_complex_chain():
    """Chain of mixed argument types."""
    pt = AlgebraicPureTensor(2, A, B)
    at = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(x, pt, at)
    assert isinstance(result, AlgebraicTensor)
    # Should have 2 terms, each with coefficient 2*x
    for term in result.terms:
        if isinstance(term, AlgebraicPureTensor):
            assert term.coeff == 2 * x


def test_pure_tensor_algebraic_tensor_scalar():
    pt = AlgebraicPureTensor(A, B)
    at = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(pt, at, x)
    assert isinstance(result, AlgebraicTensor)
    for term in result.terms:
        if isinstance(term, AlgebraicPureTensor):
            assert term.coeff == x


def test_multiple_pure_tensors_with_algebraic():
    pt1 = AlgebraicPureTensor(A, B)
    pt2 = AlgebraicPureTensor(C, D)
    at = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(pt1, pt2, at)
    assert isinstance(result, AlgebraicTensor)
    # 1 * 1 * 2 = 2 terms
    assert len(result.terms) == 2


def test_algebraic_tensor_product_preserves_factor_order():
    """Factors should appear in the order they are given."""
    result = algebraic_tensor_product(A, B, E)
    assert result.factors == (A, B, E)


def test_algebraic_tensor_product_with_transpose():
    """Transpose factors are handled correctly."""
    result = algebraic_tensor_product(A.T, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A.T, B)


# ---------------------------------------------------------------------------
# Identity and zero edge cases
# ---------------------------------------------------------------------------

def test_S_One_only_matrix():
    """S.One as first arg with single matrix."""
    result = algebraic_tensor_product(S.One, A)
    assert result == A


def test_S_One_with_two_matrices():
    """S.One with two matrices goes through main logic."""
    result = algebraic_tensor_product(S.One, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == S.One
    assert result.factors == (A, B)


def test_S_Zero_with_matrices():
    """S.Zero coefficient produces zero tensor."""
    result = algebraic_tensor_product(S.Zero, A, B)
    assert isinstance(result, AlgebraicZeroTensor)
    assert result.shape == ((3, 4), (4, 5))


# ---------------------------------------------------------------------------
# Hash and equality of results
# ---------------------------------------------------------------------------

def test_result_hash_consistency():
    r1 = algebraic_tensor_product(A, B)
    r2 = algebraic_tensor_product(A, B)
    assert hash(r1) == hash(r2)
    assert r1 == r2


def test_result_equality_different_order():
    r1 = algebraic_tensor_product(A, B)
    r2 = algebraic_tensor_product(B, A)
    assert r1 != r2


def test_result_set_membership():
    r1 = algebraic_tensor_product(A, B)
    r2 = algebraic_tensor_product(A, B)
    r3 = algebraic_tensor_product(C, D)
    s = {r1, r2, r3}
    assert len(s) == 2


# ---------------------------------------------------------------------------
# Deep nesting and edge cases
# ---------------------------------------------------------------------------

def test_nested_algebraic_tensor_product():
    """algebraic_tensor_product result used as argument to another call."""
    pt1 = algebraic_tensor_product(A, B)
    pt2 = algebraic_tensor_product(pt1, E)
    assert isinstance(pt2, AlgebraicPureTensor)
    assert pt2.factors == (A, B, E)


def test_product_with_itself():
    pt = AlgebraicPureTensor(A, B)
    result = algebraic_tensor_product(pt, pt)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, B, A, B)


def test_algebraic_tensor_product_itself():
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    result = algebraic_tensor_product(at, at)
    assert isinstance(result, AlgebraicTensor)
    # 2 * 2 = 4 terms
    assert len(result.terms) == 4


def test_product_empty_factors_from_all_coeffs():
    """When all args are coefficients but 3+, produces ValueError."""
    # 3 symbols -> no factors -> ValueError
    raises(ValueError, lambda: algebraic_tensor_product(x, y, x + y))


def test_product_with_single_factor_algebraic_tensor():
    """AlgebraicTensor with single term unwraps, but product still works."""
    at = AlgebraicTensor(AlgebraicPureTensor(A, B))
    # Single term unwraps to PureTensor
    assert isinstance(at, AlgebraicPureTensor)
    result = algebraic_tensor_product(at, E)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, B, E)


def test_product_preserves_all_factor_types():
    """Product handles MatMul, MatrixSymbol, and transpose factors."""
    mm = MatMul(A, B, evaluate=False)
    result = algebraic_tensor_product(A, A.T, mm, E)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, A.T, mm, E)


def test_product_coeff_from_multiple_sources():
    """Coefficients from PureTensor, scalar args, and symbols all multiply."""
    pt = AlgebraicPureTensor(2, A, B)
    result = algebraic_tensor_product(3, x, pt, y, E)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == 6 * x * y
    assert result.factors == (A, B, E)


def test_product_algebraic_tensor_single_term():
    """AlgebraicTensor that unwraps to single term."""
    at = AlgebraicTensor(AlgebraicPureTensor(A, B))
    result = algebraic_tensor_product(at)
    assert isinstance(result, AlgebraicPureTensor)


def test_product_zero_coeff_in_algebraic_tensor():
    """Zero coefficient in AlgebraicTensor term produces zero tensor."""
    at = AlgebraicTensor(AlgebraicPureTensor(0, A, B))
    assert isinstance(at, AlgebraicZeroTensor)
    result = algebraic_tensor_product(at, E)
    assert isinstance(result, AlgebraicZeroTensor)


def test_product_with_immutable_matrix():
    M = ImmutableDenseMatrix([[1, 2, 3], [4, 5, 6]])
    result = algebraic_tensor_product(M, A)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.shape == ((2, 3), (3, 4))


def test_product_coeff_one_unwraps_single_factor():
    """Single factor with coeff 1 unwraps."""
    result = algebraic_tensor_product(S.One, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == S.One
    assert result.factors == (A, B)


def test_product_multiple_zero_coeffs():
    """Multiple zero coefficients multiply to zero."""
    result = algebraic_tensor_product(0, 0, A, B)
    assert isinstance(result, AlgebraicZeroTensor)


def test_product_algebraic_tensor_all_terms_cancel():
    """AlgebraicTensor where all terms cancel is a zero tensor."""
    at = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(-1, A, B))
    assert isinstance(at, AlgebraicZeroTensor)
    result = algebraic_tensor_product(at, E)
    assert isinstance(result, AlgebraicZeroTensor)


def test_product_mixed_zero_and_real_terms():
    """Zero tensor mixed with real terms in product."""
    zt = AlgebraicZeroTensor(((3, 4), (4, 5)))
    pt = AlgebraicPureTensor(A, B)
    result = algebraic_tensor_product(zt, pt)
    assert isinstance(result, AlgebraicZeroTensor)


def test_product_four_args():
    """Four matrix arguments."""
    result = algebraic_tensor_product(A, B, E, P)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.factors == (A, B, E, P)


def test_product_with_negative_symbolic_coeff():
    result = algebraic_tensor_product(-x, A, B)
    assert isinstance(result, AlgebraicPureTensor)
    assert result.coeff == -x
    assert result.factors == (A, B)


def test_product_algebraic_tensor_preserves_shape():
    at1 = AlgebraicTensor(AlgebraicPureTensor(A, B), AlgebraicPureTensor(C, D))
    at2 = AlgebraicTensor(AlgebraicPureTensor(E), AlgebraicPureTensor(F))
    result = algebraic_tensor_product(at1, at2)
    assert result.shape == ((3, 4), (4, 5), (5, 3))


def test_product_pure_tensor_shape():
    pt = AlgebraicPureTensor(2, A, B)
    result = algebraic_tensor_product(pt, E)
    assert result.shape == ((3, 4), (4, 5), (5, 3))
