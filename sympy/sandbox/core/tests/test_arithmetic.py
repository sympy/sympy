from sympy.sandbox.core import *

def test_arithmetic():
    x = Symbol('x')
    y = Symbol('y')
    assert x+1 == 1+x
    assert x+y == y+x
    assert (x+y)+3 == x+(y+3)
    assert 2*x == x*2
    assert 0*x == 0
    assert 1*x == x
    assert 2*x - x == x
