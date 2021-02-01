from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_C():
    c = CartanType("C4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -2, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 4
    assert c.simple_root(4) == Matrix([[0, 0, 0, 2]])
    assert c.roots() == 32
    assert c.basis() == 36
    assert c.lie_algebra() == "sp(8)"
    t = CartanType(['C', 3])
    assert t.dimension() == 3
    diag = "0---0---0=<=0\n1   2   3   4"
    assert c.dynkin_diagram() == diag
    assert c.positive_roots() ==[
        Matrix([[2, 0, 0, 0]]),
        Matrix([[1, 1, 0, 0]]),
        Matrix([[0, 2, 0, 0]]),
        Matrix([[1, 0, 1, 0]]),
        Matrix([[0, 1, 1, 0]]),
        Matrix([[1, 0, 0, 1]]),
        Matrix([[0, 1, 0, 1]]),
        Matrix([[0, 0, 2, 0]]),
        Matrix([[1, 0, 0, -1]]),
        Matrix([[0, 1, 0, -1]]),
        Matrix([[0, 0, 1, 1]]),
        Matrix([[1, 0, -1, 0]]),
        Matrix([[0, 1, -1, 0]]),
        Matrix([[0, 0, 1, -1]]),
        Matrix([[0, 0, 0, 2]]),
        Matrix([[1, -1, 0, 0]])]

def test_type_C6():
    '''Testing numpy backend'''
    c = CartanType("C6")
    m = Matrix([
        [ 2, -1,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0],
        [ 0,  0, -1,  2, -1,  0],
        [ 0,  0,  0, -1,  2, -1],
        [ 0,  0,  0,  0, -2,  2]])
    assert m == c.cartan_matrix()

    simpleroot0 = Matrix([[1, -1, 0, 0, 0, 0]])

    assert simpleroot0 == c.simple_roots()[0]

    # take sample for brevity
    orbit = c.orbit(simpleroot0)
    assert len(orbit) == 60

    assert orbit[0] == Matrix([[1, 1, 0, 0, 0, 0]])
    assert orbit[6] == Matrix([[-1, 0, 0, 0, 0, 1]])

    assert Matrix(c.positive_roots()) == Matrix([
        [ 2,  0,  0,  0,  0,  0],
        [ 1,  1,  0,  0,  0,  0],
        [ 0,  2,  0,  0,  0,  0],
        [ 1,  0,  1,  0,  0,  0],
        [ 0,  1,  1,  0,  0,  0],
        [ 1,  0,  0,  1,  0,  0],
        [ 0,  1,  0,  1,  0,  0],
        [ 0,  0,  2,  0,  0,  0],
        [ 1,  0,  0,  0,  1,  0],
        [ 0,  1,  0,  0,  1,  0],
        [ 0,  0,  1,  1,  0,  0],
        [ 1,  0,  0,  0,  0,  1],
        [ 0,  1,  0,  0,  0,  1],
        [ 0,  0,  1,  0,  1,  0],
        [ 0,  0,  0,  2,  0,  0],
        [ 1,  0,  0,  0,  0, -1],
        [ 0,  1,  0,  0,  0, -1],
        [ 0,  0,  1,  0,  0,  1],
        [ 0,  0,  0,  1,  1,  0],
        [ 1,  0,  0,  0, -1,  0],
        [ 0,  1,  0,  0, -1,  0],
        [ 0,  0,  1,  0,  0, -1],
        [ 0,  0,  0,  1,  0,  1],
        [ 0,  0,  0,  0,  2,  0],
        [ 1,  0,  0, -1,  0,  0],
        [ 0,  1,  0, -1,  0,  0],
        [ 0,  0,  1,  0, -1,  0],
        [ 0,  0,  0,  1,  0, -1],
        [ 0,  0,  0,  0,  1,  1],
        [ 1,  0, -1,  0,  0,  0],
        [ 0,  1, -1,  0,  0,  0],
        [ 0,  0,  1, -1,  0,  0],
        [ 0,  0,  0,  1, -1,  0],
        [ 0,  0,  0,  0,  1, -1],
        [ 0,  0,  0,  0,  0,  2],
        [ 1, -1,  0,  0,  0,  0]
    ])
