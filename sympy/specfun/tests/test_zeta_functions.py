import sys
sys.path.append(".")

import py

from sympy import *
from sympy.specfun.zeta_functions import *

def _test_bernoulli():
    assert bernoulli(0) == 1
    assert bernoulli(1) == Rational(-1,2)
    assert bernoulli(2) == Rational(1,6)
    assert bernoulli(3) == 0
    assert bernoulli(4) == Rational(-1,30)
    assert bernoulli(5) == 0
    assert bernoulli(6) == Rational(1,42)
    assert bernoulli(7) == 0
    assert bernoulli(8) == Rational(-1,30)
    assert bernoulli(10) == Rational(5,66)
    assert bernoulli(1000001) == 0

def _test_bernoulli_poly():
    x = Symbol('x')
    assert bernoulli_poly(0, x) == 1
    assert bernoulli_poly(1, x) == x-Rational(1,2)
    assert bernoulli_poly(2, x) == x**2-x+Rational(1,6)
    assert bernoulli_poly(3, x) == x**3 - (3*x**2)/2 + x/2

def _test_zeta():
    assert zeta(0) == Rational(-1,2)
    assert zeta(1) == oo
    assert zeta(2) == pi**2/6
    assert zeta(4) == pi**4/90
    assert zeta(6) == pi**6/945
    assert zeta(-1) == Rational(-1,12)
    assert zeta(-2) == 0
    assert zeta(-3) == Rational(1,120)
    assert zeta(-4) == 0
    assert zeta(-5) == Rational(-1,252)

def _test_dirichlet_eta():
    assert dirichlet_eta(0) == Rational(1,2)
    assert dirichlet_eta(-1) == Rational(1,4)
    assert dirichlet_eta(1) == log(2)
    assert dirichlet_eta(2) == pi**2/12
    assert dirichlet_eta(4) == pi**4*Rational(7,720)

def _test_harmonic():
    assert harmonic(1) == 1
    assert harmonic(2) == Rational(3,2)
    assert harmonic(3) == Rational(11,6)
    assert harmonic(4) == Rational(25,12)
    assert harmonic(3,1) == harmonic(3)
    assert harmonic(3,5) == 1 + Rational(1,2**5) + Rational(1,3**5)
    assert harmonic(10,0) == 10
    assert harmonic(oo) == oo
    assert harmonic(oo,2) == (pi**2)/6

def _test_polygamma():
    assert polygamma(0, 1) == -euler_gamma
    assert polygamma(0, 7) == Rational(49,20)-euler_gamma
    assert polygamma(0, 0) == oo
    assert polygamma(1, 1) == pi**2/6
    assert polygamma(1, 2) == pi**2/6 - 1
    assert polygamma(1, 3) == pi**2/6 - Rational(5,4)
    assert polygamma(3, 1) == pi**4 / 15
    assert polygamma(3, 5) == 6*(Rational(-22369,20736) + pi**4/90)
    assert polygamma(5, 1) == 8 * pi**6 / 63
    assert polygamma(3, Rational(1)/2) == pi**4
    x = Symbol('x')
    assert diff(polygamma(3, 7*x), x) == 7*polygamma(4, 7*x)
