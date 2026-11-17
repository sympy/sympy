from __future__ import annotations
from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix
from sympy.core.backend import Rational

def test_type_E():
    c = CartanType("E6")
    assert c.dimension() == 8
    assert c.simple_root(6) == Matrix([[1, 1, 0, 0, 0, 0, 0, 0]])
    assert c.roots() == 72
    assert c.basis() == 78
    diag = " "*8 + "2\n" + " "*8 + "0\n" + " "*8 + "|\n" + " "*8 + "|\n"
    diag += "---".join("0" for _ in range(1, 6))+"\n"
    diag += "1   " + "   ".join(str(i) for i in range(3, 7))
    assert c.dynkin_diagram() == diag
    posroots = c.positive_roots()
    assert Matrix([[1, 0, 0, 0, 1, 0, 0, 0]]) in posroots
    assert Matrix([[Rational(1,2),Rational(1,2),Rational(1,2),Rational(1,2),
                            Rational(1,2),Rational(-1,2),Rational(-1,2),Rational(1,2)]]) in posroots

def test_Cartan():
    c = CartanType("E6")
    assert c.cartan_matrix() == Matrix([
        [ 2, -1,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0],
        [ 0, -1,  2, -1,  0, -1],
        [ 0,  0, -1,  2, -1,  0],
        [ 0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  2]])

    c = CartanType("E7")
    assert c.cartan_matrix() == Matrix([
        [ 2, -1,  0,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0, -1],
        [ 0,  0, -1,  2, -1,  0,  0],
        [ 0,  0,  0, -1,  2, -1,  0],
        [ 0,  0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  0,  2]])

    c = CartanType("E8")
    assert c.cartan_matrix() == Matrix([
        [ 2, -1,  0,  0,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0,  0, -1],
        [ 0,  0, -1,  2, -1,  0,  0,  0],
        [ 0,  0,  0, -1,  2, -1,  0,  0],
        [ 0,  0,  0,  0, -1,  2, -1,  0],
        [ 0,  0,  0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  0,  0,  2]])
