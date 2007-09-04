from __future__ import division

#this module tests that sympy works with true division turned on

from sympy import Rational, Symbol, Real

def test_truediv():
    assert 1/2 != 0
    assert Rational(1)/2 != 0

def dotest(s):
    x = Symbol("x")
    y = Symbol("y")
    l = [
    Rational(2),
    Real("1.3"),
    x,
    y,
    pow(x,y)*y,
    5,
    5.5
    ]
    for x in l:
        for y in l:
            s(x,y)

def test_basic():
    def s(a,b):
        x = a
        x = +a
        x = -a
        x = a+b
        x = a-b
        x = a*b
        x = a/b
        x = a**b
    dotest(s)

def test_ibasic():
    def s(a,b):
        x = a
        x += b
        x = a
        x -= b
        x = a
        x *= b
        x = a
        x /= b
    dotest(s)
