from __future__ import annotations

from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
from sympy.tensor.algebraic.pure_tensor import PureTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
E = MatrixSymbol("E", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)
I4 = MatrixSymbol("I4", 4, 4)


def test_as_common_left_basic():
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    left, rest, right = at.as_common_left()
    assert left == ()
    assert rest is at


def test_as_common_left_no_common():
    at = PureTensor(A, C) + PureTensor(B, C)
    left, rest, right = at.as_common_left()
    assert left == ()
    assert rest is at


def test_as_common_right_basic():
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    left, rest, right = at.as_common_right()
    assert len(right) >= 1
    assert isinstance(rest, (AlgebraicTensor, PureTensor, MatrixSymbol))


def test_as_common_right_no_common():
    at = PureTensor(A, I3, D) + PureTensor(B, I3, E)
    left, rest, right = at.as_common_right()
    assert left == ()
    assert right == ()
    assert rest is at


def test_as_common_factors():
    J4 = MatrixSymbol("J4", 4, 4)
    at = PureTensor(A, I4, C) + PureTensor(A, J4, C)
    assert at.tensor_shape == ((3, 4), (4, 4), (4, 5))
    left, mid, right = at.as_common_factors()
    assert len(left) >= 1
    assert left[0] == A
    assert len(right) >= 1
    assert right[-1] == C


def test_as_common_left_identity_factors():
    at = PureTensor(A, I3, D) + PureTensor(A, I3, E)
    left, rest, right = at.as_common_left()
    assert len(left) >= 2
    assert left[0] == A
    assert left[1] == I3


def test_as_common_right_identity_factors():
    at = PureTensor(A, I4, C) + PureTensor(B, I4, C)
    left, rest, right = at.as_common_right()
    assert len(right) >= 1
