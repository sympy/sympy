import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import evalf

def test_simple_evalf():
    assert evalf(2) == 2
    assert evalf(Rational(1,2)) == 0.5
    assert evalf(exp(pi*I)).ae(-1)
    assert evalf(pi+E).ae(5.8598744820488378)
