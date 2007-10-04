import sys
sys.path.append(".")

import py

import sympy as g
from sympy import Symbol, Wild, sin, cos, exp, pi, Function, Derivative
from sympy.utilities.pytest import XFAIL

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
    assert e.subs_dict(r) == r[a]/r[b] * sin(r[b]*x)
    assert e.subs_dict(r) == 3 * sin(4*x) / 4

def test_deriv_sub_bug3():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    pat = Derivative(f(x), x, x)
    assert pat.subs(y, y**2) == Derivative(f(x), x, x)
    assert pat.subs(y, y**2) != Derivative(f(x), x)
