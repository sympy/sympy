
from sympy import *

n, k = symbols('nk', integer=True)
C0, C1, C2 = symbols('C0', 'C1', 'C2')

def test_rsolve_poly():
    assert rsolve_poly([-1, k+1], k, k) == 1
    assert rsolve_poly([-k-1, k], 1, k) == C1*k - 1
    assert rsolve_poly([-4*k-2, 1], 4*k+1, k) == -1

    #assert rsolve_poly([-1, k+1], 1, k) is None

def test_rsolve_ratio():
    pass

def test_rsolve_hyper():
    assert rsolve_hyper([2*n*(n+1), -n**2-3*n+2, n-1], 0, n) \
        == C0*factorial(n) + C1*2**n