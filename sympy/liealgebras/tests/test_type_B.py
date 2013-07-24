from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_B():
    c = CartanType("B3")
    m = Matrix(3, 3, [2, -1, 0, -1, 2, -2, 0, -1, 2])
    assert m == c.cartan_matrix()
    assert c.dimension() == 3
    assert c.roots() == 18
    assert c.simple_root(3) == [0, 0, 1]
    assert c.basis() == 3
    assert c.lie_algebra() == "so(6)"
