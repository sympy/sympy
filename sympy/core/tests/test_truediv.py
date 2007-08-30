from __future__ import division

from sympy import Rational

def test_truediv():
    assert 1/2 != 0
    assert Rational(1)/2 != 0
