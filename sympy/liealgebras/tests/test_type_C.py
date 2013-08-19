from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_C():
    c = CartanType("C4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -2, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 4
    assert c.simple_root(4) == [0, 0, 0, 2]
    assert c.roots() == 32
    assert c.basis() == 36
    assert c.lie_algebra() == "sp(8)"
    t = CartanType(['C', 3])
    assert t.dimension() == 3
