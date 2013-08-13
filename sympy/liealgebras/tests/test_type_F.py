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
    assert c.positive_roots() == {1: [1, -1, 0, 0], 2: [1, 1, 0, 0],
            3: [1, 0, -1, 0], 4: [1, 0, 1, 0], 5: [1, 0, 0, -1], 6: [1, 0, 0, 1],
            7: [0, 1, -1, 0], 8: [0, 1, 1, 0], 9: [0, 1, 0, -1],
            10: [0, 1, 0, 1], 11: [0, 0, 1, -1], 12: [0, 0, 1, 1],
            13: [1, 0, 0, 0], 14: [0, 1, 0, 0], 15: [0, 0, 1, 0],
            16: [0, 0, 0, 1], 17: [0, 0, 0, 0], 18: [0, -1, 0, 0],
            19: [0, 0, -1, 0], 20: [0, 0, 0, -1], 21: [0, 0, -1, -1],
            22: [0, -1, 0, -1], 23: [0, -1, -1, 0], 24: [0, -1, -1, -1]}
