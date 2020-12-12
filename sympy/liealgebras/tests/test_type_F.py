from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix
from sympy.core.numbers import Rational

def test_type_F():
    c = CartanType("F4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -2, 0, 0, -1, 2, -1, 0, 0, -1, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 4
    assert c.simple_root(1) == Matrix([[0, 1, -1, 0]])
    assert c.simple_root(2) == Matrix([[0, 0, 1, -1]])
    assert c.simple_root(3) == Matrix([[0, 0, 0, 1]])
    assert c.simple_root(4) == Matrix([[Rational(1,2), -Rational(1,2), -Rational(1,2), -Rational(1,2)]])
    assert c.roots() == 48
    assert c.basis() == 52
    diag = "0---0=>=0---0\n" + "   ".join(str(i) for i in range(1, 5))
    assert c.dynkin_diagram() == diag

def test_pos_roots():
    """The order of the roots can be finicky because F4 has a
    lot of roots with the same weight level."""
    c = CartanType("F4")
    hardcoded = [Matrix([[1, 1, 0, 0]]),
        Matrix([[1, 0, 1, 0]]),
        Matrix([[1, 0, 0, 1]]),
        Matrix([[1, 0, 0, 0]]),
        Matrix([[1/2, 1/2, 1/2, 1/2]]),
        Matrix([[1, 0, 0, -1]]),
        Matrix([[1/2, 1/2, 1/2, -1/2]]),
        Matrix([[1, 0, -1, 0]]),
        Matrix([[0, 1, 1, 0]]),
        Matrix([[1, -1, 0, 0]]),
        Matrix([[1/2, 1/2, -1/2, 1/2]]),
        Matrix([[1/2, -1/2, 1/2, 1/2]]),
        Matrix([[0, 1, 0, 1]]),
        Matrix([[1/2, 1/2, -1/2, -1/2]]),
        Matrix([[0, 0, 1, 1]]),
        Matrix([[0, 1, 0, 0]]),
        Matrix([[1/2, -1/2, 1/2, -1/2]]),
        Matrix([[0, 1, 0, -1]]),
        Matrix([[1/2, -1/2, -1/2, 1/2]]),
        Matrix([[0, 0, 1, 0]]),
        Matrix([[0, 0, 1, -1]]),
        Matrix([[0, 0, 0, 1]]),
        Matrix([[0, 1, -1, 0]]),
        Matrix([[1/2, -1/2, -1/2, -1/2]])]
    p = c.positive_roots()

    for r in hardcoded:
        assert r in p
