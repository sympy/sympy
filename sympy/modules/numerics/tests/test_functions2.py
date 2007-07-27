import sys
sys.path.append(".")
import py
from sympy import *
from sympy.modules.numerics import *
import math

def test_gamma():
    assert Float(gamma(1)) == 1
    assert Float(gamma(2)) == 1
    assert Float(gamma(3)) == 2
    assert Float(gamma(10)) == 362880
    assert Float(gamma(0.25)).ae('3.6256099082219083119')
    assert Float(gamma(0.0001)).ae('9999.4228832316241908')
    assert Float(gamma(300)).ae('1.0201917073881354535e612')
    assert Float(gamma(-0.5)).ae('-3.5449077018110320546')
    assert Float(gamma(-7.43)).ae('0.00026524416464197007186')
    Float.setdps(100)
    assert gamma(0.5).ae(sqrt(pi_float()))
    Float.revert()
