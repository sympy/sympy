from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_G():
    c = CartanType("G2")
    m = Matrix(2, 2, [2, -1, -3, 2])
    assert c.cartan_matrix() == m
    assert c.simple_root(2) == [1, -2, 1]
    assert c.basis() == 14
    assert c.roots() == 12
    assert c.dimension() == 3
