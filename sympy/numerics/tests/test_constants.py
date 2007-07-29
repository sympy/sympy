import sys
sys.path.append(".")
import py
from sympy import *
from sympy.numerics import *
from sympy.numerics.constants import *

pi = "3.1415926535897932384626433832795028841971693993751058209749445923078\
1640628620899862803482534211706798"

gamma = "0.5772156649015328606065120900824024310421593359399235988057672348\
84867726777664670936947063291746749516"

log2 = "0.69314718055994530941723212145817656807550013436025525412068000949\
3393621969694715605863326996418687542"

log10 = "2.3025850929940456840179914546843642076011014886287729760333279009\
6757260967735248023599720508959829834"


def test_constants():
    for p in [3, 7, 10, 15, 20, 37, 80, 100, 29]:
        Float.store()
        Float.setdps(p)
        assert pi_float() == Float(pi)
        assert gamma_float() == Float(gamma)
        assert log10_float() == Float(log10)
        assert log2_float() == Float(log2)
        Float.revert()

def test_extreme_precision():
    Float.store(); Float.setdps(1001)
    assert str(pi_float()).endswith('1989')
    Float.revert()
