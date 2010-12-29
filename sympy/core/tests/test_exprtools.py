"""Tests for tools for manipulating of large commutative expressions. """

from sympy.core.exprtools import (
    decompose_power, Factors, Term, _gcd_terms, gcd_terms)

from sympy import S, Add, sin
from sympy.abc import x, y, z, t

def test_decompose_power():
    assert decompose_power(x) == (x, 1)
    assert decompose_power(x**2) == (x, 2)
    assert decompose_power(x**(2*y)) == (x**y, 2)
    assert decompose_power(x**(2*y/3)) == (x**(y/3), 2)

def test_Factors():
    assert Factors() == Factors({})

    assert Factors().as_expr() == S.One
    assert Factors({x: 2, y: 3, sin(x): 4}).as_expr() == x**2*y**3*sin(x)**4

    a = Factors({x: 5, y: 3, z: 7})
    b = Factors({y: 4, z: 3, t: 10})

    assert a.mul(b) == a*b == Factors({x: 5, y: 7, z: 10, t: 10})

    assert a.div(b) == divmod(a, b) == (Factors({x: 5, z: 4}), Factors({y: 1, t: 10}))
    assert a.quo(b) == a/b == Factors({x: 5, z: 4})
    assert a.rem(b) == a%b == Factors({y: 1, t: 10})

    assert a.pow(3) == a**3 == Factors({x: 15, y: 9, z: 21})
    assert b.pow(3) == b**3 == Factors({y: 12, z: 9, t: 30})

    assert a.gcd(b) == Factors({y: 3, z: 3})
    assert a.lcm(b) == Factors({x: 5, y: 4, z: 7, t: 10})

    a = Factors({x: 4, y: 7, t: 7})
    b = Factors({z: 1, t: 3})

    assert a.normal(b) == (Factors({x: 4, y: 7, t: 4}), Factors({z: 1}))

def test_Term():
    a = Term(4*x*y**2/z/t**3)
    b = Term(2*x**3*y**5/t**3)

    assert a == Term(4, Factors({x: 1, y: 2}), Factors({z: 1, t: 3}))
    assert b == Term(2, Factors({x: 3, y: 5}), Factors({t: 3}))

    assert a.as_expr() == 4*x*y**2/z/t**3
    assert b.as_expr() == 2*x**3*y**5/t**3

    assert a.inv() == Term(S(1)/4, Factors({z: 1, t: 3}), Factors({x: 1, y: 2}))
    assert b.inv() == Term(S(1)/2, Factors({t: 3}), Factors({x: 3, y: 5}))

    assert a.mul(b) == a*b == Term(8, Factors({x: 4, y: 7}), Factors({z: 1, t: 6}))
    assert a.quo(b) == a/b == Term(2, Factors({}), Factors({x: 2, y: 3, z: 1}))

    assert a.pow(3) == a**3 == Term(64, Factors({x: 3, y: 6}), Factors({z: 3, t: 9}))
    assert b.pow(3) == b**3 == Term(8, Factors({x: 9, y: 15}), Factors({t: 9}))

    assert a.pow(-3) == a**(-3) == Term(S(1)/64, Factors({z: 3, t: 9}), Factors({x: 3, y: 6}))
    assert b.pow(-3) == b**(-3) == Term(S(1)/8, Factors({t: 9}), Factors({x: 9, y: 15}))

    assert a.gcd(b) == Term(2, Factors({x: 1, y: 2}), Factors({t: 3}))
    assert a.lcm(b) == Term(4, Factors({x: 3, y: 5}), Factors({z: 1, t: 3}))

    a = Term(4*x*y**2/z/t**3)
    b = Term(2*x**3*y**5*t**7)

    assert a.mul(b) == Term(8, Factors({x: 4, y: 7, t: 4}), Factors({z: 1}))

def test_gcd_terms():
    f = 2*(x + 1)*(x + 4)/(5*x**2 + 5) + (2*x + 2)*(x + 5)/(x**2 + 1)/5 + (2*x + 2)*(x + 6)/(5*x**2 + 5)

    assert _gcd_terms(f) == ((S(6)/5)*((1 + x)/(1 + x**2)), 5 + x, 1)
    assert _gcd_terms(Add.make_args(f)) == ((S(6)/5)*((1 + x)/(1 + x**2)), 5 + x, 1)

    assert gcd_terms(f) == (S(6)/5)*((1 + x)*(5 + x)/(1 + x**2))
    assert gcd_terms(Add.make_args(f)) == (S(6)/5)*((1 + x)*(5 + x)/(1 + x**2))

    assert gcd_terms(0) == 0
    assert gcd_terms(1) == 1
    assert gcd_terms(x) == x
