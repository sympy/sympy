from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_A3():
    c = CartanType("A3")
    m = Matrix(3, 3, [2, -1, 0, -1, 2, -1, 0, -1, 2])
    assert m == c.cartan_matrix()
    assert c.basis() == 8
    assert c.roots() == 12
    assert c.dimension() == 4
    assert c.simple_root(1) == Matrix([[1, -1, 0, 0]])
    assert c.highest_root() == Matrix([[1, 0, 0, -1]])
    assert c.lie_algebra() == "su(4)"
    diag = "0---0---0\n1   2   3"
    assert c.dynkin_diagram() == diag
    assert c.positive_roots() == [Matrix([[1, 0, 0, -1]]),
        Matrix([[0, 1, 0, -1]]),
        Matrix([[1, 0, -1, 0]]),
        Matrix([[0, 1, -1, 0]]),
        Matrix([[0, 0, 1, -1]]),
        Matrix([[1, -1, 0, 0]])]

    flag = False
    try:
        c.simple_root(10)
    except ValueError:
        flag = True
    assert flag


def test_type_A6():
    '''Testing numpy backend'''
    c = CartanType("A6")
    m = Matrix([
        [2, -1, 0, 0, 0, 0],
        [-1, 2, -1, 0, 0, 0],
        [0, -1, 2, -1, 0, 0],
        [0, 0, -1, 2, -1, 0],
        [0, 0, 0, -1, 2, -1],
        [0, 0, 0, 0, -1, 2]
    ])
    assert m == c.cartan_matrix()

    simpleroot0 = Matrix([[1, -1, 0, 0, 0, 0, 0]])

    assert simpleroot0 == c.simple_roots()[0]

    # take sample for brevity
    orbit = c.orbit(simpleroot0)
    assert len(orbit) == 42

    assert orbit[0] == Matrix([[1, 0, 0, 0, 0, 0, -1]])
    assert orbit[6] == Matrix([[-1, 1, 0, 0, 0, 0, 0]])

    assert Matrix(c.positive_roots()) == Matrix([
        [1, 0, 0, 0, 0, 0, -1],
        [0, 1, 0, 0, 0, 0, -1],
        [1, 0, 0, 0, 0, -1, 0],
        [0, 1, 0, 0, 0, -1, 0],
        [0, 0, 1, 0, 0, 0, -1],
        [1, 0, 0, 0, -1, 0, 0],
        [0, 1, 0, 0, -1, 0, 0],
        [0, 0, 1, 0, 0, -1, 0],
        [0, 0, 0, 1, 0, 0, -1],
        [1, 0, 0, -1, 0, 0, 0],
        [0, 1, 0, -1, 0, 0, 0],
        [0, 0, 1, 0, -1, 0, 0],
        [0, 0, 0, 1, 0, -1, 0],
        [0, 0, 0, 0, 1, 0, -1],
        [1, 0, -1, 0, 0, 0, 0],
        [0, 1, -1, 0, 0, 0, 0],
        [0, 0, 1, -1, 0, 0, 0],
        [0, 0, 0, 1, -1, 0, 0],
        [0, 0, 0, 0, 1, -1, 0],
        [0, 0, 0, 0, 0, 1, -1],
        [1, -1, 0, 0, 0, 0, 0]])
