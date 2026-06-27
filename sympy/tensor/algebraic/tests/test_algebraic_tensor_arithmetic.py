from __future__ import annotations

from sympy.core.singleton import S
from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
I3 = MatrixSymbol("I3", 3, 3)


def test_algebraic_tensor_mul_scalar():
    at = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    result = at * 3
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))


def test_algebraic_tensor_rmul_scalar():
    at = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    result = 3 * at
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))


def test_algebraic_tensor_mul_symbol():
    from sympy import Symbol
    x = Symbol("x")
    at = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    result = at * x
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_rmul_symbol():
    from sympy import Symbol
    x = Symbol("x")
    at = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    result = x * at
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_mul_negative_one():
    at = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    result = at * S.NegativeOne
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_mul_preserves_terms_count():
    at = AlgebraicPureTensor(A, I3) + AlgebraicPureTensor(B, I3)
    original_count = len(at.args)
    result = at * 5
    assert len(result.args) == original_count
