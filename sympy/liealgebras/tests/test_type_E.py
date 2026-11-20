from __future__ import annotations
from sympy.liealgebras.cartan_type import CartanType
from sympy import Matrix, S
from sympy.core.backend import Rational

def test_type_E():
    c = CartanType("E6")
    m = Matrix(6, 6, [2, 0, -1, 0, 0, 0, 0, 2, 0, -1, 0, 0,
        -1, 0, 2, -1, 0, 0, 0, -1, -1, 2, -1, 0, 0, 0, 0,
        -1, 2, -1, 0, 0, 0, 0, -1, 2])
    assert c.cartan_matrix() == m
    assert c.dimension() == 8
    assert c.simple_root(6) == [0, 0, 0, -1, 1, 0, 0, 0]
    assert c.roots() == 72
    assert c.basis() == 78
    diag = " "*8 + "2\n" + " "*8 + "0\n" + " "*8 + "|\n" + " "*8 + "|\n"
    diag += "---".join("0" for i in range(1, 6))+"\n"
    diag += "1   " + "   ".join(str(i) for i in range(3, 7))
    assert c.dynkin_diagram() == diag
    posroots = c.positive_roots()
    assert posroots[8] == [1, 0, 0, 0, 1, 0, 0, 0]
    assert posroots[21] == [Rational(1,2),Rational(1,2),Rational(1,2),Rational(1,2),
                            Rational(1,2),Rational(-1,2),Rational(-1,2),Rational(1,2)]


def test_simpleroots():
    c = CartanType("E6")
    s = c.simple_roots()

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
    s = c.simple_roots()

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
    s = c.simple_roots()

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
