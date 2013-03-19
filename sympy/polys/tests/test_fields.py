"""Test sparse rational functions. """

from sympy.polys.fields import field
from sympy.polys.domains import ZZ, QQ, RR
from sympy.polys.monomialtools import lex, grlex

from sympy.utilities.pytest import raises
from sympy.core import Symbol, symbols
from sympy import sqrt, Rational

def test_FracField___hash__():
    F, x, y, z = field("x,y,z", QQ)
    assert hash(F)

def test_FracElement___hash__():
    F, x, y, z = field("x,y,z", QQ)
    assert hash(x*y/z)

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

def test_FracElement_from_expr():
    x, y, z = symbols("x,y,z")
    F, X, Y, Z = field((x, y, z), ZZ)

    f = F.from_expr(1)
    assert f == 1 and isinstance(f, F.dtype)

    f = F.from_expr(Rational(3, 7))
    assert f == F(3)/7 and isinstance(f, F.dtype)

    f = F.from_expr(x)
    assert f == X and isinstance(f, F.dtype)

    f = F.from_expr(Rational(3,7)*x)
    assert f == 3*X/7 and isinstance(f, F.dtype)

    f = F.from_expr(1/x)
    assert f == 1/X and isinstance(f, F.dtype)

    f = F.from_expr(x*y*z)
    assert f == X*Y*Z and isinstance(f, F.dtype)

    f = F.from_expr(x*y/z)
    assert f == X*Y/Z and isinstance(f, F.dtype)

    f = F.from_expr(x*y*z + x*y + x)
    assert f == X*Y*Z + X*Y + X and isinstance(f, F.dtype)

    f = F.from_expr((x*y*z + x*y + x)/(x*y + 7))
    assert f == (X*Y*Z + X*Y + X)/(X*Y + 7) and isinstance(f, F.dtype)

    f = F.from_expr(x**3*y*z + x**2*y**7 + 1)
    assert f == X**3*Y*Z + X**2*Y**7 + 1 and isinstance(f, F.dtype)

    raises(ValueError, lambda: F.from_expr(2**x))
    raises(ValueError, lambda: F.from_expr(7*x + sqrt(2)))

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

    assert x + F.ring.gens[0] == F.ring.gens[0] + x == 2*x

    F, x,y = field("x,y", ZZ)
    assert x + 3 == 3 + x
    assert x + QQ(3,7) == QQ(3,7) + x == (7*x + 3)/7

def test_FracElement___sub__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y
    assert f - g == (-x + y)/(x*y)

    assert x - F.ring.gens[0] == F.ring.gens[0] - x == 0

    F, x,y = field("x,y", ZZ)
    assert x - 3 == -(3 - x)
    assert x - QQ(3,7) == -(QQ(3,7) - x) == (7*x - 3)/7

def test_FracElement___mul__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y
    assert f*g == g*f == 1/(x*y)

    assert x*F.ring.gens[0] == F.ring.gens[0]*x == x**2

    F, x,y = field("x,y", ZZ)
    assert x*3 == 3*x
    assert x*QQ(3,7) == QQ(3,7)*x == 3*x/7

    Fuv, u,v = field("u,v", ZZ);
    Fxyzt, x,y,z,t = field("x,y,z,t", Fuv.to_domain())

    f = ((u + 1)*x*y + 1)/((v - 1)*z - t*u*v - 1)
    assert dict(f.numer) == {(1, 1, 0, 0): u + 1, (0, 0, 0, 0): 1}
    assert dict(f.denom) == {(0, 0, 1, 0): v - 1, (0, 0, 0, 1): -u*v, (0, 0, 0, 0): -1}

def test_FracElement___div__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y
    assert f/g == y/x

    assert x/F.ring.gens[0] == F.ring.gens[0]/x == 1

    F, x,y = field("x,y", ZZ)
    assert x*3 == 3*x
    assert x/QQ(3,7) == (QQ(3,7)/x)**-1 == 7*x/3

def test_FracElement___pow__():
    F, x,y = field("x,y", QQ)

    f, g = 1/x, 1/y

    assert f**3 == 1/x**3
    assert g**3 == 1/y**3

    assert (f*g)**3 == 1/(x**3*y**3)
    assert (f*g)**-3 == (x*y)**3
