from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix



def test_type_D():
    c = CartanType("D4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -1, -1, 0, -1, 2, 0, 0, -1, 0 , 2])
    assert c.cartan_matrix() == m
    assert c.basis() == 6
    assert c.lie_algebra() == "so(8)"
    assert c.roots() == 24
    assert c.simple_root(3) == [0, 0, 1, -1]
