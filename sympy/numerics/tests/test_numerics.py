import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *

def test_import():
    assert Float(3) + Float(6) == 9
    assert isinstance(exp(3), Float)
    assert isinstance(log(3), Float)
