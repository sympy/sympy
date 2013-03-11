"""Test sparse rational functions. """

from sympy.polys.fields import field
from sympy.polys.domains import ZZ, QQ, RR, ZZ_python
from sympy.polys.monomialtools import lex, grlex

from sympy.utilities.pytest import raises
from sympy.core import Symbol, symbols

def test_FracField___repr__():
    assert repr(field("x", ZZ, lex)[0]) == "FracField((x,), ZZ, LexOrder())"
    assert repr(field("x,y", QQ, grlex)[0]) == "FracField((x, y), QQ, GradedLexOrder())"
    assert repr(field("x,y,z", ZZ["t"], lex)[0]) == "FracField((x, y, z), ZZ[t], LexOrder())"

def test_FracField___str__():
    assert str(field("x", ZZ, lex)[0]) == "Rational function field in x over ZZ with lex order"
    assert str(field("x,y", QQ, grlex)[0]) == "Rational function field in x, y over QQ with grlex order"
    assert str(field("x,y,z", ZZ["t"], lex)[0]) == "Rational function field in x, y, z over ZZ[t] with lex order"

def test_FracField___hash__():
    F, x, y, z = field("x,y,z", QQ)
    assert hash(F)

def test_FracElement___hash__():
    F, x, y, z = field("x,y,z", QQ)
    assert hash(x*y/z)

def test_FracElement___repr__():
    F, x, y = field("x,y", ZZ_python())
    assert repr((3*x**2*y + 1)/(x - y**2)) == "FracElement(FracField((x, y), ZZ, LexOrder()), [((2, 1), 3), ((0, 0), 1)], [((1, 0), 1), ((0, 2), -1)])"

def test_FracElement___str__():
    F, x, y = field("x,y", ZZ_python())
    assert str((3*x**2*y + 1)/(x - y**2)) == "(3*x**2*y + 1)/(x - y**2)"

def test_FracElement_copy():
    F, x, y, z = field("x,y,z", ZZ)

    f = x*y/3*z
    g = f.copy()

    assert f == g
    g.numer[(1, 1, 1)] = 7
    assert f != g

def test_FracElement_as_expr():
    F, x, y, z = field("x,y,z", ZZ)
    f = (3*x**2*y - x*y*z)/(7*z**3 + 1)

    X, Y, Z = F.symbols
    g = (3*X**2*Y - X*Y*Z)/(7*Z**3 + 1)

    assert f != g
    assert f.as_expr() == g

    X, Y, Z = symbols("x,y,z")
    g = (3*X**2*Y - X*Y*Z)/(7*Z**3 + 1)

    assert f != g
    assert f.as_expr(X, Y, Z) == g

    raises(ValueError, lambda: f.as_expr(X))

def test_FracElement___neg__():
    F, x,y = field("x,y", QQ)

    f = (7*x - 9)/y
    g = (-7*x + 9)/y

    assert -f == g
    assert -g == f

def test_FracElement___add__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y
    assert f + g == g + f == (x + y)/(x*y)

def test_FracElement___sub__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y
    assert f - g == (-x + y)/(x*y)

def test_FracElement___mul__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y
    assert f*g == g*f == 1/(x*y)

def test_FracElement___div__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y
    assert f/g == y/x

def test_FracElement___pow__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y

    assert f**3 == 1/x**3
    assert g**3 == 1/y**3

    assert (f*g)**3 == 1/(x**3*y**3)
    assert (f*g)**-3 == (x*y)**3
