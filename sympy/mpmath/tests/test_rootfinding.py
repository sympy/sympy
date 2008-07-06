from sympy.mpmath import *

def test_secant():
    mp.dps = 15
    assert secant(lambda x: 4*x-3, 5).ae(0.75)
    assert secant(sin, 3).ae(pi)
    assert secant(sin, 3, 3.14).ae(pi)
    assert secant(lambda x: x*x+1, 2+2j).ae(1j)
