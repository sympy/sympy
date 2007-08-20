
from sympy import *

n, k = symbols('nk', integer=True)
C0, C1, C2 = symbols('C0', 'C1', 'C2')

def test_rsolve_poly():
    assert rsolve_poly([-1, -1, 1], 0, n) == 0
    assert rsolve_poly([-1, -1, 1], 1, n) == -1

    assert rsolve_poly([-1, n+1], n, n) == 1
    assert rsolve_poly([-1, 1], n, n) == C0 + (n**2 - n)/2
    assert rsolve_poly([-n-1, n], 1, n) == C1*n - 1
    assert rsolve_poly([-4*n-2, 1], 4*n+1, n) == -1

    assert rsolve_poly([-1, 1], n**5 + n**3, n) == C0 - n**3 / 2 - n**5 / 2 + n**2 / 6 + n**6 / 6 + 2*n**4 / 3

def test_rsolve_ratio():
    assert rsolve_ratio([-2*n**3+n**2+2*n-1, 2*n**3+n**2-6*n, -2*n**3-11*n**2-18*n-9, 2*n**3+13*n**2+22*n+8], 0, n) == C2*(2*n-3)/(n-1)/(n+1)/2

def test_rsolve_hyper():
    assert rsolve_hyper([-1, -1, 1], 0, n) == C0*(S.Half + S.Half*sqrt(5))**n + C1*(S.Half - S.Half*sqrt(5))**n

    assert rsolve_hyper([n**2-2, -2*n-1, 1], 0, n) == C0*rf(sqrt(2), n) + C1*rf(-sqrt(2), n)

    assert rsolve_hyper([n**2-k, -2*n-1, 1], 0, n) == C0*rf(sqrt(k), n) + C1*rf(-sqrt(k), n)

    assert rsolve_hyper([2*n*(n+1), -n**2-3*n+2, n-1], 0, n) == C0*factorial(n) + C1*2**n

    assert rsolve_hyper([n + 2, -(2*n + 3)*(17*n**2 + 51*n + 39), n + 1], 0, n) == 0

    assert rsolve_hyper([-n-1, -1, 1], 0, n) == 0
