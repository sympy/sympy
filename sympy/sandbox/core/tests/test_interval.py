import py

from sympy.sandbox.core import *
from sympy.sandbox.core.interval import *

def test_membership():
    assert Interval(2) == Interval(2, 2)
    a = Interval(2, 3)
    assert 2.5 in a
    assert 2 in a
    assert 3 in a
    assert 1 not in a
    assert 4 not in a

def test_mul():
    assert Interval(2, 3) * Interval(4, 5) == Interval(8, 15)
    assert Interval(-2, 3) * Interval(4, 5) == Interval(-10, 15)
    assert Interval(-2, 3) * Interval(-4, 5) == Interval(-12, 15)
    assert Interval(-2, 3) * Interval(-5, -4) == Interval(-15, 10)

def test_div():
    assert Interval(15, 30) / Interval(3, 5) == Interval(3, 10)
    assert Interval(-15, 30) / Interval(3, 5) == Interval(-5, 10)
    assert Interval(-15, 30) / Interval(-5, -3) == Interval(-10, 5)
    # assert 3 / Interval(1, 2) == Interval(Fraction(3,2), 3)
    Float.setdps(15)
    a = Interval(1) / Float(3)
    b = Interval(-1) / Float(3)
    assert a == Interval(0.33333333333333331, 0.33333333333333337)
    assert b == Interval(-0.33333333333333337, -0.33333333333333331)
    py.test.raises(ZeroDivisionError, "3 / Interval(-2, 2)")

def test_add():
    assert Interval(4, 5) - 2 == Interval(2, 3)
    assert Interval(2, 3) - Interval(2, 3) == Interval(-1, 1)
    assert 2 - Interval(3,5) == Interval(-3, -1)
    assert Interval(4, 5) + 2 == Interval(6, 7)
    assert 2 + Interval(4, 5) == Interval(6, 7)
    assert Interval(2, 3) + Interval(2, 3) == Interval(4, 6)

def test_pi():
    def acot(x):
        x = Interval(Float(x))
        p = Float(10) ** (-Float.getdps() - 2)
        s = w = Interval(1)/x
        x = x*x
        n = 3
        while 1:
            w = w / x
            term = w / n
            if n & 2: s -= term
            else:     s += term
            if term.a <= p:
                break
            n += 2
        return s
    Float.setdps(15)
    pi = 16*acot(5) - 4*acot(239)
    a = float(pi.a)
    b = float(pi.b)
    assert a < 3.1415926535897932 < b
