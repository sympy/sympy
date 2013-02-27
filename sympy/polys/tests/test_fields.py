"""Test sparse rational functions. """

from sympy.polys.fields import field
from sympy.polys.domains import ZZ, QQ
from sympy.utilities.pytest import raises

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
