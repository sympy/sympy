import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *
from sympy.numerics.functions import *
from sympy.numerics.functions2 import *
import math

def test_gamma():
    Float.store()
    Float.setprec(53)
    assert gamma(0.25).ae('3.6256099082219083119')
    assert gamma(0.0001).ae('9999.4228832316241908')
    assert gamma(300).ae('1.0201917073881354535e612')
    assert gamma(-0.5).ae('-3.5449077018110320546')
    assert gamma(-7.43).ae('0.00026524416464197007186')
    assert gamma(Rational(1,2)) == gamma(0.5)
    assert gamma(Rational(-7,3)).ae(gamma(Float(7)/3))
    # Should always give exact factorials when they can
    # be represented as Floats under the current working precision
    fact = 1
    for i in range(1, 18):
        assert gamma(i) == fact
        fact *= i
    Float.setprec(1000)
    for i in range(18, 105):
        assert gamma(i) == fact
        fact *= i
    Float.setdps(100)
    assert gamma(0.5).ae(sqrt(pi_float()))
    Float.revert()
