import sys
sys.path.append(".")
import py
from sympy import *
from sympy.modules.numerics import *
import math

def test_sqrt():
    for i in range(1000):
        assert sqrt(Float(i**2)) == i
    # These should round identically
    for x in [0, 1e-7, 0.1, 0.5, 1, 2, 3, 4, 5, 0.333, 76.19]:
        assert sqrt(Float(x)) == float(x)**0.5

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

def test_log():
    assert log(1) == 0
    for x in [0.5, 1.5, 2.0, 3.0, 100, 10**50, 1e-50]:
        assert log(x) == math.log(x)
        assert log(x, x) == 1
    assert log(1024, 2) == 10
    assert log(10**1234, 10) == 1234

def test_atan():
    import math
    assert atan2(1,1).ae(math.atan2(1,1))
    assert atan2(1,-1).ae(math.atan2(1,-1))
    assert atan2(-1,-1).ae(math.atan2(-1,-1))
    assert atan2(-1,1).ae(math.atan2(-1,1))
    assert atan2(-1,0).ae(math.atan2(-1,0))
    assert atan2(1,0).ae(math.atan2(1,0))
