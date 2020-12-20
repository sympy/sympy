from sympy.functions.elementary.miscellaneous import sqrt
from sympy.core.numbers import Rational
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
    s = c.simple_roots()

    assert s[0].tolist() == [[1, -1, 0, 0, 0, 0]]
    assert s[3].tolist() == [[0, 0, 0, 1, 1, 0]]
    assert s[4].tolist() == [[
        -Rational(1,2),
        -Rational(1,2),
        -Rational(1,2),
        -Rational(1,2),
        -Rational(1,2),
        sqrt(3)*Rational(1,2)]]
    assert s[-1].tolist() == [[0, 0, 0, 1, -1, 0]]
    c = CartanType("E7")
    s = c.simple_roots()

    assert s[0].tolist() == [[1, -1, 0, 0, 0, 0, 0]]
    assert s[4].tolist() == [[0, 0, 0, 0, 1, 1, 0]]
    assert s[5].tolist() == [[
        -Rational(1,2),
        -Rational(1,2),
        -Rational(1,2),
        -Rational(1,2),
        -Rational(1,2),
        -Rational(1,2),
        sqrt(2)*Rational(1,2)]]
    assert s[-1].tolist() == [[0, 0, 0, 0, 1, -1, 0]]

def test_pos_roots():
    """The order of the roots can be finicky because E has a
    lot of roots with the same weight level."""

    c = CartanType("E6")
    assert c.positive_roots() == [
        Matrix([[S.Half, S.Half, S.Half, S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[S.Half, S.Half, S.Half, -S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[S.Half, S.Half, -S.Half, S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[S.Half, S.Half, -S.Half, -S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[S.Half, -S.Half, S.Half, S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[-S.Half, S.Half, S.Half, S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[1, 1, 0, 0, 0, 0]]),
        Matrix([[S.Half, -S.Half, S.Half, -S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[-S.Half, S.Half, S.Half, -S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[1, 0, 1, 0, 0, 0]]),
        Matrix([[S.Half, -S.Half, -S.Half, S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[0, 1, 1, 0, 0, 0]]),
        Matrix([[-S.Half, S.Half, -S.Half, S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[1, 0, 0, 1, 0, 0]]),
        Matrix([[S.Half, -S.Half, -S.Half, -S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[0, 1, 0, 1, 0, 0]]),
        Matrix([[-S.Half, S.Half, -S.Half, -S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[-S.Half, -S.Half, S.Half, S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[1, 0, 0, 0, -1, 0]]),
        Matrix([[1, 0, 0, 0, 1, 0]]),
        Matrix([[0, 1, 0, 0, -1, 0]]),
        Matrix([[0, 1, 0, 0, 1, 0]]),
        Matrix([[0, 0, 1, 1, 0, 0]]),
        Matrix([[-S.Half, -S.Half, S.Half, -S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[1, 0, 0, -1, 0, 0]]),
        Matrix([[0, 1, 0, -1, 0, 0]]),
        Matrix([[0, 0, 1, 0, -1, 0]]),
        Matrix([[0, 0, 1, 0, 1, 0]]),
        Matrix([[-S.Half, -S.Half, -S.Half, S.Half, S.Half, sqrt(3)/2]]),
        Matrix([[1, 0, -1, 0, 0, 0]]),
        Matrix([[0, 1, -1, 0, 0, 0]]),
        Matrix([[0, 0, 1, -1, 0, 0]]),
        Matrix([[0, 0, 0, 1, -1, 0]]),
        Matrix([[0, 0, 0, 1, 1, 0]]),
        Matrix([[-S.Half, -S.Half, -S.Half, -S.Half, -S.Half, sqrt(3)/2]]),
        Matrix([[1, -1, 0, 0, 0, 0]])]
