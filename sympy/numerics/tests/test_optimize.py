import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *
from sympy.numerics.optimize import *

def _test_polyroots():
    x = Symbol('x')
    rs = polyroots(4*(x-3)*(x+I)*(x-4-5*I))
    assert rs[0][0].ae(ComplexFloat(0, -1))
    assert rs[0][1].ae(ComplexFloat(3, 0))
    assert rs[0][2].ae(ComplexFloat(4, 5))

def test_bisect():
    a, b = bisect(sin, 3, 4)
    assert sin(a).ae(0)
    assert sin(b).ae(0)
    assert sin(a)*sin(b) <= 0

def test_secant():
    assert secant(lambda x: x**2-4, 3.7).ae(2)
    Float.setdps(100)
    assert secant(sin, 3).ae(pi_float())
    Float.setdps(15)
