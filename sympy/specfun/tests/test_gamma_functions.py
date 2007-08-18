
from sympy.core import *
from sympy.specfun import gamma, lowergamma, uppergamma

def test_gamma():
    
    assert gamma(1) == 1
    assert gamma(2) == 1
    assert gamma(3) == 2
    
    assert gamma(Rational(1,2)) == sqrt(pi)
    
    assert uppergamma(4, 0) == 6
