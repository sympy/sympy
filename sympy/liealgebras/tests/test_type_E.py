from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_E():
    c = CartanType("E6")
    m = Matrix(6, 6, [2, 0, -1, 0, 0, 0, 0, 2, 0, -1, 0, 0,
        -1, 0, 2, -1, 0, 0, 0, -1, -1, 2, -1, 0, 0, 0, 0,
        -1, 2, -1, 0, 0, 0, 0, -1, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 8
    assert c.simple_root(6) == [0, 0, 0, -1, 1, 0, 0, 0]
    assert c.roots() == 72
    assert c.basis() == 78
