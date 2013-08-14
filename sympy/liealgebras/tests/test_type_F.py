from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_F():
    c = CartanType("F4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2, -1, 0, 0, -1, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 4
    assert c.simple_root(3) == [0, 0, 0, 1]
    assert c.simple_root(4) == [-0.5, -0.5, -0.5, -0.5]
    assert c.roots() == 48
    assert c.basis() == 52
