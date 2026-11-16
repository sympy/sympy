from __future__ import annotations
from sympy.liealgebras.root_system import RootSystem
from sympy.liealgebras.type_a import TypeA
from sympy.matrices import Matrix

def test_root_system():
    c = RootSystem("A3")
    assert c.cartan_type == TypeA(3)
    assert c.simple_roots() ==  {1: Matrix([[1, -1, 0, 0]]), 2: Matrix([[0, 1, -1, 0]]), 3: Matrix([[0, 0, 1, -1]])}
    assert c.root_space() == "alpha[1] + alpha[2] + alpha[3]"
    assert c.cartan_matrix() == Matrix([[ 2, -1,  0], [-1,  2, -1], [ 0, -1,  2]])
    assert c.dynkin_diagram() == "0---0---0\n1   2   3"
    assert c.add_simple_roots(1, 2) == Matrix([[1, 0, -1, 0]])
    assert c.all_roots() == {1: Matrix([[1, -1, 0, 0]]), 2: Matrix([[1, 0, -1, 0]]),
            3: Matrix([[1, 0, 0, -1]]), 4: Matrix([[0, 1, -1, 0]]), 5: Matrix([[0, 1, 0, -1]]),
            6: Matrix([[0, 0, 1, -1]]), 7: Matrix([[-1, 1, 0, 0]]), 8: Matrix([[-1, 0, 1, 0]]),
            9: Matrix([[-1, 0, 0, 1]]), 10: Matrix([[0, -1, 1, 0]]),
            11: Matrix([[0, -1, 0, 1]]), 12: Matrix([[0, 0, -1, 1]])}
    assert c.add_as_roots(Matrix([[1, 0, -1, 0]]), Matrix([[0, 0, 1, -1]])) == Matrix([[1, 0, 0, -1]])
