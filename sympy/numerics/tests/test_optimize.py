import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *
from sympy.numerics.optimize import *

def test_polyroots():
    x = Symbol('x')
    rs = polyroots(4*(x-3)*(x+I)*(x-4-5*I))
    assert rs[0][0].ae(ComplexFloat(0, -1))
    assert rs[0][1].ae(ComplexFloat(3, 0))
    assert rs[0][2].ae(ComplexFloat(4, 5))
