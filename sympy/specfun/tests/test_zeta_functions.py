import sys
sys.path.append(".")

import py

from sympy import *
from sympy.specfun.zeta_functions import *

def test_zeta():
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

def test_dirichlet_eta():
    assert dirichlet_eta(0) == Rational(1,2)
    assert dirichlet_eta(-1) == Rational(1,4)
    assert dirichlet_eta(1) == log(2)
    assert dirichlet_eta(2) == pi**2/12
    assert dirichlet_eta(4) == pi**4*Rational(7,720)

def test_harmonic():
    assert harmonic(1,1) == 1
    assert harmonic(2,1) == Rational(3,2)
    assert harmonic(3,1) == Rational(11,6)
    assert harmonic(4,1) == Rational(25,12)
    # assert harmonic(3,1) == harmonic(3)
    assert harmonic(3,5) == 1 + Rational(1,2**5) + Rational(1,3**5)
    assert harmonic(10,0) == 10
    assert harmonic(oo,1) == oo
    assert harmonic(oo,2) == (pi**2)/6

def test_polygamma():
    assert polygamma(0, 1) == -EulerGamma
    assert polygamma(0, 7) == Rational(49,20)-EulerGamma
    assert polygamma(0, 0) == oo
    assert polygamma(1, 1) == pi**2/6
    assert polygamma(1, 2) == pi**2/6 - 1
    assert polygamma(1, 3) == pi**2/6 - Rational(5,4)
    assert polygamma(3, 1) == pi**4 / 15
    assert polygamma(3, 5) == 6*(Rational(-22369,20736) + pi**4/90)
    assert polygamma(5, 1) == 8 * pi**6 / 63
    assert polygamma(3, Rational(1)/2) == pi**4
    x = Symbol('x')
    assert polygamma(3, 7*x).diff(x) == 7*polygamma(4, 7*x)
