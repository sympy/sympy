"""Tests of transforms of quantum expressions for Mul and Pow."""

from sympy.core.symbol import symbols
from sympy.core.sympify import sympify
from sympy.testing.pytest import warns_deprecated_sympy

from sympy.physics.quantum.operator import (
    Operator, OuterProduct
)
from sympy.physics.quantum.state import Ket, Bra
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.tensorproduct import TensorProduct


k1 = Ket('k1')
k2 = Ket('k2')
k3 = Ket('k3')
b1 = Bra('b1')
b2 = Bra('b2')
b3 = Bra('b3')
A = Operator('A')
B = Operator('B')
C = Operator('C')
x, y, z = symbols('x y z')


def test_bra_ket():
    assert b1*k1 == InnerProduct(b1, k1)
    assert k1*b1 == OuterProduct(k1, b1)
    # Test priority of inner product
    assert OuterProduct(k1, b1)*k2 == InnerProduct(b1, k2)*k1
    assert b1*OuterProduct(k1, b2) == InnerProduct(b1, k1)*b2


def test_tensor_product():
    assert k1*k1 == TensorProduct(k1, k1)
    assert b1*b1 == TensorProduct(b1, b1)
    assert k1*TensorProduct(k2, k3) == TensorProduct(k1, k2, k3)
    assert b1*TensorProduct(b2, b3) == TensorProduct(b1, b2, b3)
    assert TensorProduct(k2, k3)*k1 == TensorProduct(k2, k3, k1)
    assert TensorProduct(b2, b3)*b1 == TensorProduct(b2, b3, b1)

    assert TensorProduct(A, B, C)*TensorProduct(k1, k2, k3) == \
        TensorProduct(A*k1, B*k2, C*k3)
    assert TensorProduct(b1, b2, b3)*TensorProduct(A, B, C) == \
        TensorProduct(b1*A, b2*B, b3*C)
    assert TensorProduct(b1, b2, b3)*TensorProduct(k1, k2, k3) == \
        InnerProduct(b1, k1)*InnerProduct(b2, k2)*InnerProduct(b3, k3)
    assert TensorProduct(b1, b2, b3)*TensorProduct(A, B, C)*TensorProduct(k1, k2, k3) == \
        TensorProduct(b1*A*k1, b2*B*k2, b3*C*k3)


def test_outer_product():
    assert OuterProduct(k1, b1)*OuterProduct(k2, b2) == \
        InnerProduct(b1, k2)*OuterProduct(k1, b2)


def test_compound():
    e1 = b1*A*B*k1*b2*k2*b3
    assert e1 == InnerProduct(b2, k2)*b1*A*B*OuterProduct(k1, b3)

    e2 = TensorProduct(k1, k2)*TensorProduct(b1, b2)
    assert e2 == TensorProduct(
        OuterProduct(k1, b1),
        OuterProduct(k2, b2)
    )
