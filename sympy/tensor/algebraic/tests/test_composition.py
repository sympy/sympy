from __future__ import annotations

from contextlib import contextmanager

@contextmanager
def raises(expected_exc):
    """Minimal pytest.raises replacement for environments without pytest."""
    try:
        yield
    except expected_exc:
        pass
    except Exception as e:
        raise AssertionError(f"Expected {expected_exc.__name__} but got {type(e).__name__}: {e}")
    else:
        raise AssertionError(f"Expected {expected_exc.__name__} but no exception was raised")

from sympy.matrices.expressions import MatrixSymbol
from sympy.matrices.expressions.matexpr import MatMul
from sympy.tensor.algebraic.algebraic_pure_tensor import (
    AlgebraicPureTensor,
    algebraic_tensor_product,
    compose_algebraic_pure_tensors,
)


# --- Fixtures ---

# Left tensor factors: A(3,4), C(4,5)
A = MatrixSymbol("A", 3, 4)
C = MatrixSymbol("C", 4, 5)

# Right tensor factors: B(4,6), D(5,7)
B = MatrixSymbol("B", 4, 6)
D = MatrixSymbol("D", 5, 7)

# Single factor test
G = MatrixSymbol("G", 2, 3)
H = MatrixSymbol("H", 3, 4)

# Three factor test
M1 = MatrixSymbol("M1", 2, 3)
M2 = MatrixSymbol("M2", 4, 5)
M3 = MatrixSymbol("M3", 6, 7)

N1 = MatrixSymbol("N1", 3, 8)
N2 = MatrixSymbol("N2", 5, 9)
N3 = MatrixSymbol("N3", 7, 10)


def test_compose_basic():
    """Basic composition of two 2-factor tensors."""
    left = AlgebraicPureTensor(A, C)
    right = AlgebraicPureTensor(B, D)

    result = compose_algebraic_pure_tensors(left, right)

    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 2

    factors = result.factors
    assert factors[0].shape == (3, 6)
    assert factors[1].shape == (4, 7)
    assert result.tensor_shape == ((3, 6), (4, 7))


def test_compose_single_factor():
    """Composition when both tensors have a single factor (unwrapped matrices)."""
    # AlgebraicPureTensor(G) unwraps to G directly
    left = G
    right = H

    result = compose_algebraic_pure_tensors(left, right)

    # Result is a MatMul of shape (2, 4)
    assert result.shape == (2, 4)
    assert isinstance(result, MatMul)


def test_compose_three_factors():
    """Composition with three factors per tensor."""
    left = AlgebraicPureTensor(M1, M2, M3)
    right = AlgebraicPureTensor(N1, N2, N3)

    result = compose_algebraic_pure_tensors(left, right)

    assert isinstance(result, AlgebraicPureTensor)
    assert result.num_factors == 3
    assert result.factors[0].shape == (2, 8)
    assert result.factors[1].shape == (4, 9)
    assert result.factors[2].shape == (6, 10)


def test_compose_with_coefficients():
    """Composition combines coefficients from both tensors."""
    from sympy import Symbol
    from sympy.core.numbers import Integer
    alpha = Symbol("alpha")

    left = AlgebraicPureTensor(alpha, A, C)
    right = AlgebraicPureTensor(Integer(2), B, D)

    result = compose_algebraic_pure_tensors(left, right)

    assert isinstance(result, AlgebraicPureTensor)
    assert result._get_coeff() == 2 * alpha
    assert len(result.factors) == 2
    assert result.factors[0].shape == (3, 6)
    assert result.factors[1].shape == (4, 7)


def test_compose_coeff_zero():
    """Composition fails for AlgebraicZeroTensor (no factors to compose)."""
    from sympy.core.singleton import S

    left = AlgebraicPureTensor(S.Zero, A, C)
    right = AlgebraicPureTensor(B, D)

    from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
    assert isinstance(left, AlgebraicZeroTensor)
    # AlgebraicZeroTensor has .shape so it passes the bare-matrix check,
    # but then fails the factor count check
    with raises(ValueError):
        compose_algebraic_pure_tensors(left, right)


def test_compose_mismatched_factor_count():
    """Composition fails when tensors have different numbers of factors."""
    left = AlgebraicPureTensor(A, C)
    right = AlgebraicPureTensor(B, D, MatrixSymbol("X", 7, 8))

    with raises(ValueError):
        compose_algebraic_pure_tensors(left, right)


def test_compose_incompatible_inner_dim():
    """Composition fails when inner dimensions don't match for a factor pair."""
    B2 = MatrixSymbol("B2", 5, 6)

    left = AlgebraicPureTensor(A, C)
    right = AlgebraicPureTensor(B2, D)

    with raises(ValueError):
        compose_algebraic_pure_tensors(left, right)


def test_compose_not_pure_tensor_left():
    """Composition fails when left argument is neither PureTensor nor matrix-like."""
    from sympy import Symbol
    x = Symbol("x")
    right = AlgebraicPureTensor(B, D)

    with raises(TypeError):
        compose_algebraic_pure_tensors(x, right)


def test_compose_not_pure_tensor_right():
    """Composition fails when right argument is neither PureTensor nor matrix-like."""
    from sympy import Symbol
    x = Symbol("x")
    left = AlgebraicPureTensor(A, C)

    with raises(TypeError):
        compose_algebraic_pure_tensors(left, x)


def test_compose_matmul_result():
    """Verify that composed factors are MatMul objects."""
    left = AlgebraicPureTensor(A, C)
    right = AlgebraicPureTensor(B, D)

    result = compose_algebraic_pure_tensors(left, right)

    factors = result.factors
    assert isinstance(factors[0], MatMul)
    assert isinstance(factors[1], MatMul)


def test_compose_mixed_bare_and_pure():
    """Composition of a bare matrix with a PureTensor of matching factor count."""
    # Left is bare (2,3), right is PureTensor with one factor (3,4)
    left = G
    right = H  # Also bare, but tested as PureTensor-equivalent

    result = compose_algebraic_pure_tensors(left, right)
    assert result.shape == (2, 4)


def test_compose_preserves_order():
    """Matrix multiplication order is preserved: left[j] @ right[j]."""
    left = AlgebraicPureTensor(A, C)
    right = AlgebraicPureTensor(B, D)

    result = compose_algebraic_pure_tensors(left, right)

    factors = result.factors
    # factors[0] = A @ B (not B @ A)
    assert factors[0].args[0] is A
    assert factors[0].args[1] is B
    # factors[1] = C @ D
    assert factors[1].args[0] is C
    assert factors[1].args[1] is D