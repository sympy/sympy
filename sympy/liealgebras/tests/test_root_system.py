from sympy.liealgebras.root_system import RootSystem
from sympy.matrices import Matrix

def test_old():
    c = RootSystem("A3")

    assert c.dynkin_diagram() == "0---0---0\n1   2   3"
    assert c.add_simple_roots(1, 2) == Matrix([[0, 1, 0, -1]])
    assert c.add_as_roots([1, 0, -1, 0], [0, 0, 1, -1]) == Matrix([[1, 0, 0, -1]])


def test_new():
    c = RootSystem("A3")
    assert c.root_space() == "alpha[1] + alpha[2] + alpha[3]"
    assert c.cartan_matrix() == Matrix([[ 2, -1,  0], [-1,  2, -1], [ 0, -1,  2]])
    assert c.simple_roots() ==  [Matrix([[1, -1, 0, 0]]), Matrix([[0, 1, -1, 0]]), Matrix([[0, 0, 1, -1]])]
    assert c.all_roots() == [Matrix([[1, 0, 0, -1]]),
        Matrix([[0, 1, 0, -1]]),
        Matrix([[1, 0, -1, 0]]),
        Matrix([[0, 1, -1, 0]]),
        Matrix([[0, 0, 1, -1]]),
        Matrix([[1, -1, 0, 0]]),
        Matrix([[0, 0, 0, 0]]),
        Matrix([[0, 0, 0, 0]]),
        Matrix([[0, 0, 0, 0]]),
        Matrix([[-1, 1, 0, 0]]),
        Matrix([[0, 0, -1, 1]]),
        Matrix([[0, -1, 1, 0]]),
        Matrix([[-1, 0, 1, 0]]),
        Matrix([[0, -1, 0, 1]]),
        Matrix([[-1, 0, 0, 1]])]
