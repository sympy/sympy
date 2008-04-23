import py

import sympy as g
from sympy import Symbol, Wild, sin, cos, exp, sqrt, pi, Function, Derivative,\
        abc, Integer, Eq, symbols

def test_subs():
    n3=g.Rational(3)
    n2=g.Rational(2)
    n6=g.Rational(6)
    x=g.Symbol("x")
    c=g.Symbol("c")
    e=x
    assert str(e) == "x"
    e=e.subs(x,n3)
    assert str(e) == "3"

    e=2*x
    assert e == 2*x
    e=e.subs(x,n3)
    assert str(e) == "6"

    e=(g.sin(x)**2).diff(x)
    assert e == 2*g.sin(x)*g.cos(x)
    e=e.subs(x,n3)
    assert e == 2*g.cos(n3)*g.sin(n3)

    e=(g.sin(x)**2).diff(x)
    assert e == 2*g.sin(x)*g.cos(x)
    e=e.subs(g.sin(x),g.cos(x))
    assert e == 2*g.cos(x)**2

    assert exp(pi).subs(exp, sin) == 0
    assert cos(exp(pi)).subs(exp, sin) == 1

def test_logexppow():   # no eval()
    x = g.Symbol("x")
    w = g.Symbol("dummy :)")
    e = (3**(1+x)+2**(1+x))/(3**x+2**x)
    assert e.subs(2**x, w) != e
    assert e.subs(g.exp(x*g.log(g.Rational(2))),w) != e

def test_bug():
    x1=g.Symbol("x1")
    x2=g.Symbol("x2")
    y=x1*x2
    y.subs(x1,g.Real(3.0))

def test_subbug1():
    x=g.Symbol("x")
    e=(x**x).subs(x,1)
    e=(x**x).subs(x,1.0)

def test_subbug2():
    # Ensure this does not cause infinite recursion
    x = g.Symbol('x')
    assert g.Real(7.7).epsilon_eq(abs(x).subs(x, -7.7))

def test_dict():
    x = Symbol('x')
    a,b,c = map(Wild, 'abc')

    f = 3*cos(4*x)
    r = f.match(a*cos(b*x))
    assert r == {a: 3, b: 4}
    e =  a/b * sin(b*x)
    assert e._subs_dict(r) == r[a]/r[b] * sin(r[b]*x)
    assert e._subs_dict(r) == 3 * sin(4*x) / 4


def test_dict_ambigous():   # see #467
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    f = x*exp(x)
    g = z*exp(z)

    df= {x:y, exp(x): y}
    dg= {z:y, exp(z): y}

    assert f._subs_dict(df) == y**2
    assert g._subs_dict(dg) == y**2

    # and this is how order can affect the result
    assert f .subs(x,y) .subs(exp(x),y)  == y*exp(y)
    assert f .subs(exp(x),y) .subs(x,y)  == y**2


def test_deriv_sub_bug3():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    pat = Derivative(f(x), x, x)
    assert pat.subs(y, y**2) == Derivative(f(x), x, x)
    assert pat.subs(y, y**2) != Derivative(f(x), x)

def test_equality_subs1():
    f = Function("f")
    x = abc.x
    eq = Eq(f(x)**2, x)
    res = Eq(Integer(16), x)
    assert eq.subs(f(x), 4) == res

def test_equality_subs2():
    f = Function("f")
    x = abc.x
    eq = Eq(f(x)**2, 16)
    assert bool(eq.subs(f(x), 3)) == False
    assert bool(eq.subs(f(x), 4)) == True

def test_issue643():
    x = Symbol('x')
    y = Symbol('y')

    e = sqrt(x)*exp(y)
    assert e.subs(sqrt(x), 1)   == exp(y)

def test_subs_dict1():
    x, y = symbols('xy')
    assert (1+x*y).subs(x, pi) == 1 + pi*y
    assert (1+x*y).subs({x:pi, y:2}) == 1 + 2*pi

def test_subs_dict2():
    x = Symbol('x')
    a,b,c = map(Wild, 'abc')

    f = 3*cos(4*x)
    r = f.match(a*cos(b*x))
    assert r == {a: 3, b: 4}
    e =  a/b * sin(b*x)
    assert e.subs(r) == r[a]/r[b] * sin(r[b]*x)
    assert e.subs(r) == 3 * sin(4*x) / 4

def test_add():
    a, b, c, d, x = abc.a, abc.b, abc.c, abc.d, abc.x
    assert (a**2 - b - c).subs(a**2 - b, d) == d - c
    assert (a**2 - c).subs(a**2 - c, d) == d
    assert (a**2 - b - c).subs(a**2 - c, d) in [d - b, a**2 - b - c]
    assert (a**2 - x - c).subs(a**2 - c, d) == d - x
