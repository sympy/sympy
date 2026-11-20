from __future__ import annotations
from sympy.liealgebras.cartan_type import CartanType, Standard_Cartan
from sympy.testing.pytest import raises
from sympy import eye, Matrix

def test_Standard_Cartan():
    c = CartanType("A4")
    assert c.rank() == 4
    assert c.series == "A"
    m = Standard_Cartan("A", 2)
    assert m.rank() == 2
    assert m.series == "A"
    b = CartanType("B12")
    assert b.rank() == 12
    assert b.series == "B"

    # Issue 29537
    raises(ValueError, lambda: CartanType("X5"))

def test_CartanType():
    c = CartanType("A4")

    omega = c.omega_matrix()
    cocar = c.cocartan_matrix()

    assert cocar * omega.T == eye(4)

    c2 = CartanType("A2")
    assert c2.omega_matrix() == Matrix([
                        ['2/3', '-1/3', '-1/3'],
                        ['1/3',  '1/3', '-2/3']])

    assert c2.fundamental_weight(1) == Matrix([['2/3', '-1/3', '-1/3']])
