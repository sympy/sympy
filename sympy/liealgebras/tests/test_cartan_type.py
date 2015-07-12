from __future__ import division, print_function

from sympy.liealgebras.cartan_type import CartanType, Standard_Cartan

def test_Standard_Cartan():
    c = CartanType("A4")
    assert c.rank() == 4
    assert c.series == "A"
    m = Standard_Cartan("A", 2)
    assert m.rank() == 2
    assert c.series == "A"
