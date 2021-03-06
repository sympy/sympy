from sympy import Matrix, S
from sympy.liealgebras.cartan_type import CartanType

def test_type_E():
    c = CartanType("E6")
    assert c.dimension() == 6
    assert c.roots() == 72
    assert c.basis() == 78
    diag = " "*8 + "6\n" + " "*8 + "0\n" + " "*8 + "|\n" + " "*8 + "|\n"
    diag += "---".join("0" for i in range(1, 6))+"\n"
    diag += "1   " + "   ".join(str(i) for i in range(2, 6))
    assert c.dynkin_diagram() == diag


def test_simpleroots():
    c = CartanType("E6")
    s = c.simple_roots

    assert Matrix(s) == Matrix([
        [S.Half, -S.Half, -S.Half,
        -S.Half, -S.Half, -S.Half,
        -S.Half, S.Half],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0]])


    c = CartanType("E7")
    s = c.simple_roots

    assert Matrix(s) == Matrix([
        [S.Half, -S.Half, -S.Half,
        -S.Half, -S.Half, -S.Half,
        -S.Half, S.Half],
       [-1, 1, 0, 0, 0, 0, 0, 0],
       [0, -1, 1, 0, 0, 0, 0, 0],
       [0, 0, -1, 1, 0, 0, 0, 0],
       [0, 0, 0, -1, 1, 0, 0, 0],
       [0, 0, 0, 0, -1, 1, 0, 0],
       [1, 1, 0, 0, 0, 0, 0, 0]])

    c = CartanType("E8")
    s = c.simple_roots

    assert Matrix(s) == Matrix([
        [S.Half, -S.Half, -S.Half,
        -S.Half, -S.Half, -S.Half,
        -S.Half, S.Half],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [0, 0, 0, 0, -1, 1, 0, 0],
        [0, 0, 0, 0, 0, -1, 1, 0],
        [1, 1, 0, 0, 0, 0, 0, 0]])

def test_Cartan():
    c = CartanType("E6")
    s = c.cartan_matrix()

    assert Matrix(s) == Matrix([
        [ 2, -1,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0],
        [ 0, -1,  2, -1,  0, -1],
        [ 0,  0, -1,  2, -1,  0],
        [ 0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  2]])


    c = CartanType("E7")
    s = c.cartan_matrix()

    assert Matrix(s) == Matrix([
        [ 2, -1,  0,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0, -1],
        [ 0,  0, -1,  2, -1,  0,  0],
        [ 0,  0,  0, -1,  2, -1,  0],
        [ 0,  0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  0,  2]])

    c = CartanType("E8")
    s = c.cartan_matrix()

    assert Matrix(s) == Matrix([
        [ 2, -1,  0,  0,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0,  0, -1],
        [ 0,  0, -1,  2, -1,  0,  0,  0],
        [ 0,  0,  0, -1,  2, -1,  0,  0],
        [ 0,  0,  0,  0, -1,  2, -1,  0],
        [ 0,  0,  0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  0,  0,  2]])

def test_orbit():
    """These results are compared to outputs of
    Mathematica's implementation of lie algebras"""
    # ignoring on large tests without numpy
    try:
        import numpy # noqa
    except ImportError:
        return

    c = CartanType("E6")
    s = c.simple_roots
    orbit = c.orbit(s[0])

    assert len(orbit) == 72


    # Slow tests
    # using e6 logic to test for now
    # c = CartanType("E7")
    # s = c.simple_roots
    # orbit = c.orbit(s[0])

    # assert len(orbit) == 126


    # c = CartanType("E8")
    # s = c.simple_roots
    # orbit = c.orbit(s[0])

    # assert len(orbit) == 240
