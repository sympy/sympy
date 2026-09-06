from __future__ import annotations
from sympy.liealgebras.cartan_type import CartanType, Standard_Cartan
from sympy.testing.pytest import raises

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


def test_valid_cartan_type():
    c = CartanType("F4")
    assert c.rank() == 4
    assert c.series == "F"


def test_invalid_cartan_type_raises():
    with raises(ValueError):
        CartanType("Z6")
