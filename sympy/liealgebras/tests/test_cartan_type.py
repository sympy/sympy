from sympy.liealgebras.cartan_type import CartanType, Standard_Cartan
from sympy.matrices import Matrix
from sympy.core import Basic

def test_Standard_Cartan():
    c = CartanType("A4")
    assert c.rank() == 4
    assert c.series == "A"
    m = Standard_Cartan("A", 2)
    assert m.rank() == 2
    assert c.series == "A"
