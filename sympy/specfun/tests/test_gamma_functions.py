
from sympy import *

def test_gamma():
    
    assert gamma(1) == 1
    assert gamma(2) == 1
    assert gamma(3) == 2
    
    assert gamma(Rational(1,2)) == sqrt(pi)
    
    assert uppergamma(4, 0) == 6
