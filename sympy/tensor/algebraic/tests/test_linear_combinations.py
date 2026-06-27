from __future__ import annotations

from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
from sympy.tensor.algebraic.algebraic_pure_tensor import AlgebraicPureTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
I3 = MatrixSymbol("I3", 3, 3)


def test_linear_combination_puretensors():
    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, C)
    combo = 2 * pt1 + 3 * pt2
    assert isinstance(combo, AlgebraicTensor)
    assert combo.tensor_shape == ((3, 4), (4, 5))


def test_linear_combination_with_subtraction():
    pt1 = AlgebraicPureTensor(A, C)
    pt2 = AlgebraicPureTensor(B, C)
    combo = 2 * pt1 - 3 * pt2
    assert isinstance(combo, AlgebraicTensor)
    assert combo.tensor_shape == ((3, 4), (4, 5))


def test_linear_combination_chained():
    pt1 = AlgebraicPureTensor(A, I3)
    pt2 = AlgebraicPureTensor(B, I3)
    D2 = MatrixSymbol("D2", 3, 4)
    pt3 = AlgebraicPureTensor(D2, I3)
    combo = pt1 + 2 * pt2 - 3 * pt3
    assert isinstance(combo, AlgebraicTensor)


def test_scaled_puretensor_in_algebraictensor():
    pt = AlgebraicPureTensor(A, C)
    at = AlgebraicTensor(2 * pt, 3 * pt)
    assert isinstance(at, AlgebraicTensor)
    for arg in at.args:
        if isinstance(arg, AlgebraicPureTensor):
            assert arg.factors == (A, C)
