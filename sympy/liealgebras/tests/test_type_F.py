from __future__ import annotations
from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix
from sympy.core.backend import S

def test_type_F():
    c = CartanType("F4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2, -1, 0, 0, -1, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 4
    assert c.simple_root(1) == Matrix([[0, 1, -1, 0]])
    assert c.simple_root(2) == Matrix([[0, 0, 1, -1]])
    assert c.simple_root(3) == Matrix([[0, 0, 0, 1]])
    assert c.simple_root(4) == Matrix([[S.Half, -S.Half, -S.Half, -S.Half]])
    assert c.roots() == 48
    assert c.basis() == 52
    diag = "0---0=>=0---0\n" + "   ".join(str(i) for i in range(1, 5))
    assert c.dynkin_diagram() == diag
    assert c.positive_roots() == {1: Matrix([[1, -1, 0, 0]]), 2: Matrix([[1, 1, 0, 0]]), 3: Matrix([[1, 0, -1, 0]]),
            4: Matrix([[1, 0, 1, 0]]), 5: Matrix([[1, 0, 0, -1]]), 6: Matrix([[1, 0, 0, 1]]), 7: Matrix([[0, 1, -1, 0]]),
            8: Matrix([[0, 1, 1, 0]]), 9: Matrix([[0, 1, 0, -1]]), 10: Matrix([[0, 1, 0, 1]]), 11: Matrix([[0, 0, 1, -1]]),
            12: Matrix([[0, 0, 1, 1]]), 13: Matrix([[1, 0, 0, 0]]), 14: Matrix([[0, 1, 0, 0]]), 15: Matrix([[0, 0, 1, 0]]),
            16: Matrix([[0, 0, 0, 1]]), 17: Matrix([[S.Half, S.Half, S.Half, S.Half]]), 18: Matrix([[S.Half, -S.Half, S.Half, S.Half]]),
            19: Matrix([[S.Half, S.Half, -S.Half, S.Half]]), 20: Matrix([[S.Half, S.Half, S.Half, -S.Half]]), 21: Matrix([[S.Half, S.Half, -S.Half, -S.Half]]),
            22: Matrix([[S.Half, -S.Half, S.Half, -S.Half]]), 23: Matrix([[S.Half, -S.Half, -S.Half, S.Half]]),
            24: Matrix([[S.Half, -S.Half, -S.Half, -S.Half]])}
