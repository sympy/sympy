import sys
sys.path.append(".")
import py
#from sympy import *
from sympy.numerics import *
from sympy.numerics.functions import *
import math
import cmath

def test_sqrt():
    for i in range(1000):
        assert sqrt(Float(i**2)) == i
    # These should round identically
    for x in [0, 1e-7, 0.1, 0.5, 1, 2, 3, 4, 5, 0.333, 76.19]:
        assert sqrt(Float(x)) == float(x)**0.5
    assert sqrt(-1) == 1j
    assert sqrt(-2).ae(cmath.sqrt(-2))
    assert sqrt(-3).ae(cmath.sqrt(-3))
    assert sqrt(-100).ae(cmath.sqrt(-100))
    assert sqrt(1j).ae(cmath.sqrt(1j))
    assert sqrt(-1j).ae(cmath.sqrt(-1j))
    assert sqrt(math.pi + math.e*1j).ae(cmath.sqrt(math.pi + math.e*1j))
    assert sqrt(math.pi - math.e*1j).ae(cmath.sqrt(math.pi - math.e*1j))

def test_hypot():
    assert hypot(0, 0) == 0
    assert hypot(0, 0.33) == Float(0.33)
    assert hypot(0.33, 0) == Float(0.33)
    assert hypot(-0.33, 0) == Float(0.33)
    assert hypot(3, 4) == Float(5)

def test_exp():
    assert exp(0) == 1
    assert exp(10000).ae(Float('8.8068182256629215873e4342'))
    assert exp(-10000).ae(Float('1.1354838653147360985e-4343'))
    assert exp(log2_float() * Float(10)).ae(1024)
    assert exp(2+2j).ae(cmath.exp(2+2j))

def _test_log():
    assert log(1) == 0
    for x in [0.5, 1.5, 2.0, 3.0, 100, 10**50, 1e-50]:
        assert log(x) == math.log(x)
        assert log(x, x) == 1
    assert log(1024, 2) == 10
    assert log(10**1234, 10) == 1234
    assert log(2+2j).ae(cmath.log(2+2j))

def test_trig_basic():
    for x in (range(100) + range(-100,0)):
        t = x / 4.1
        assert cos(Float(t)).ae(math.cos(t))
        assert sin(Float(t)).ae(math.sin(t))
        assert tan(Float(t)).ae(math.tan(t))
    assert sin(1+1j).ae(cmath.sin(1+1j))
    assert sin(-4-3.6j).ae(cmath.sin(-4-3.6j))

def test_trig_hard():
    assert sin(Float(10**50, 150)).ae(-0.7896724934293100827)
    assert cos(Float(10**50, 150)).ae(-0.6135286082336635622)
    assert sin(1e-6).ae(9.999999999998333e-007)
    assert cos(1e-6).ae(0.9999999999995)

def test_atan():
    import math
    assert atan(-2.3).ae(math.atan(-2.3))
    assert atan2(1,1).ae(math.atan2(1,1))
    assert atan2(1,-1).ae(math.atan2(1,-1))
    assert atan2(-1,-1).ae(math.atan2(-1,-1))
    assert atan2(-1,1).ae(math.atan2(-1,1))
    assert atan2(-1,0).ae(math.atan2(-1,0))
    assert atan2(1,0).ae(math.atan2(1,0))
