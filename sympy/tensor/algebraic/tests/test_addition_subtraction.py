from __future__ import annotations

from sympy.matrices.expressions import MatrixSymbol
from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
from sympy.tensor.algebraic.pure_tensor import PureTensor
from sympy.tensor.algebraic.zero_tensor import ZeroTensor


# --- Fixtures ---

A = MatrixSymbol("A", 3, 4)
B = MatrixSymbol("B", 3, 4)
C = MatrixSymbol("C", 4, 5)
D = MatrixSymbol("D", 3, 5)
I3 = MatrixSymbol("I3", 3, 3)


def test_pure_tensor_sub_is_algebraictensor():
    pt1 = PureTensor(A, I3)
    pt2 = PureTensor(B, I3)
    result = pt1 - pt2
    assert isinstance(result, AlgebraicTensor)
    assert result.tensor_shape == ((3, 4), (3, 3))


def test_pure_tensor_sub_self():
    pt = PureTensor(A, I3)
    result = pt - pt
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_add():
    pt1 = PureTensor(A, I3)
    pt2 = PureTensor(B, I3)
    result = pt1 + pt2
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_radd():
    pt = PureTensor(A, I3)
    result = pt + PureTensor(B, I3)
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_sub():
    pt1 = PureTensor(A, I3)
    pt2 = PureTensor(B, I3)
    result = pt1 - pt2
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_rsub():
    pt = PureTensor(A, I3)
    result = pt - PureTensor(B, I3)
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_sub():
    at1 = PureTensor(A, I3) + PureTensor(B, I3)
    at2 = PureTensor(A, I3) + PureTensor(B, I3)
    result = at1 - at2
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_sub_puretensor():
    at = PureTensor(A, I3) + PureTensor(B, I3)
    pt = PureTensor(A, I3)
    result = at - pt
    assert isinstance(result, AlgebraicTensor)


def test_pure_tensor_sub_algebraictensor():
    pt = PureTensor(A, I3)
    at = PureTensor(A, I3) + PureTensor(B, I3)
    result = pt - at
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_rsub():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    result = pt_b - pt_a
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_radd():
    pt_a = PureTensor(A, I3)
    pt_b = PureTensor(B, I3)
    at = pt_a + pt_b
    result = PureTensor(A, I3) + at
    assert isinstance(result, AlgebraicTensor)


def test_algebraic_tensor_add_method():
    at1 = PureTensor(A, I3) + PureTensor(B, I3)
    at2 = PureTensor(A, I3) + PureTensor(B, I3)
    result = at1 + at2
    assert isinstance(result, AlgebraicTensor)


def test_zerotensor_sub_puretensor():
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = z - pt
    assert isinstance(result, AlgebraicTensor)


def test_puretensor_sub_zerotensor():
    z = ZeroTensor(((3, 4), (3, 3)))
    pt = PureTensor(A, I3)
    result = pt - z
    assert isinstance(result, AlgebraicTensor)
