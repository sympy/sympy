from sympy import *
from sympy.physics.units import *

def test_units():
    assert (3*m/s * year) / km == 94608
    assert foot / meter == Rational('0.3048')
