from __future__ import division
from sympy import S, sqrt, symbols
from sympy.physics.quantum.cg import Wigner3j, CG, cg_simp
from sympy.physics.quantum.kronecker import KroneckerDelta


def test_cg_simp():
    j, m1, m1p, m2, m2p = symbols('j m1 m1p m2 m2p')
    # Test Varshalovich 8.7.1 Eq 1
    a = CG(S(1)/2,S(1)/2,0,0,S(1)/2,S(1)/2)
    b = CG(S(1)/2,-S(1)/2,0,0,S(1)/2,-S(1)/2)
    c = CG(1,1,0,0,1,1)
    d = CG(1,0,0,0,1,0)
    e = CG(1,-1,0,0,1,-1)
    assert cg_simp(a+b) == 2
    assert cg_simp(c+d+e) == 3
    assert cg_simp(a+b+c+d+e) == 5
    assert cg_simp(a+b+c) == 2+c
    assert cg_simp(2*a+b) == 2+a
    assert cg_simp(2*c+d+e) == 3+c
    assert cg_simp(5*a+5*b) == 10
    assert cg_simp(5*c+5*d+5*e) == 15
    assert cg_simp(-a-b) == -2
    assert cg_simp(-c-d-e) == -3
    assert cg_simp(-6*a-6*b) == -12
    assert cg_simp(-4*c-4*d-4*e) == -12
    a = CG(S(1)/2,S(1)/2,j,0,S(1)/2,S(1)/2)
    b = CG(S(1)/2,-S(1)/2,j,0,S(1)/2,-S(1)/2)
    c = CG(1,1,j,0,1,1)
    d = CG(1,0,j,0,1,0)
    e = CG(1,-1,j,0,1,-1)
    assert cg_simp(a+b) == 2*KroneckerDelta(j,0)
    assert cg_simp(c+d+e) == 3*KroneckerDelta(j,0)
    assert cg_simp(a+b+c+d+e) == 5*KroneckerDelta(j,0)
    assert cg_simp(a+b+c) == 2*KroneckerDelta(j,0)+c
    assert cg_simp(2*a+b) == 2*KroneckerDelta(j,0)+a
    assert cg_simp(2*c+d+e) == 3*KroneckerDelta(j,0)+c
    assert cg_simp(5*a+5*b) == 10*KroneckerDelta(j,0)
    assert cg_simp(5*c+5*d+5*e) == 15*KroneckerDelta(j,0)
    assert cg_simp(-a-b) == -2*KroneckerDelta(j,0)
    assert cg_simp(-c-d-e) == -3*KroneckerDelta(j,0)
    assert cg_simp(-6*a-6*b) == -12*KroneckerDelta(j,0)
    assert cg_simp(-4*c-4*d-4*e) == -12*KroneckerDelta(j,0)
    # Test Varshalovich 8.7.1 Eq 2
    a = CG(S(1)/2,S(1)/2,S(1)/2,-S(1)/2,0,0)
    b = CG(S(1)/2,-S(1)/2,S(1)/2,S(1)/2,0,0)
    c = CG(1,1,1,-1,0,0)
    d = CG(1,0,1,0,0,0)
    e = CG(1,-1,1,1,0,0)
    assert cg_simp(a-b) == sqrt(2)
    assert cg_simp(c-d+e) == sqrt(3)
    assert cg_simp(a-b+c-d+e) == sqrt(2)+sqrt(3)
    assert cg_simp(a-b+c) == sqrt(2)+c
    assert cg_simp(2*a-b) == sqrt(2)+a
    assert cg_simp(2*c-d+e) == sqrt(3)+c
    assert cg_simp(5*a-5*b) == 5*sqrt(2)
    assert cg_simp(5*c-5*d+5*e) == 5*sqrt(3)
    assert cg_simp(-a+b) == -sqrt(2)
    assert cg_simp(-c+d-e) == -sqrt(3)
    assert cg_simp(-6*a+6*b) == -6*sqrt(2)
    assert cg_simp(-4*c+4*d-4*e) == -4*sqrt(3)
    a = CG(S(1)/2,S(1)/2,S(1)/2,-S(1)/2,j,0)
    b = CG(S(1)/2,-S(1)/2,S(1)/2,S(1)/2,j,0)
    c = CG(1,1,1,-1,j,0)
    d = CG(1,0,1,0,j,0)
    e = CG(1,-1,1,1,j,0)
    assert cg_simp(a-b) == sqrt(2)*KroneckerDelta(j,0)
    assert cg_simp(c-d+e) == sqrt(3)*KroneckerDelta(j,0)
    assert cg_simp(a-b+c-d+e) == sqrt(2)*KroneckerDelta(j,0)+sqrt(3)*KroneckerDelta(j,0)
    assert cg_simp(a-b+c) == sqrt(2)*KroneckerDelta(j,0)+c
    assert cg_simp(2*a-b) == sqrt(2)*KroneckerDelta(j,0)+a
    assert cg_simp(2*c-d+e) == sqrt(3)*KroneckerDelta(j,0)+c
    assert cg_simp(5*a-5*b) == 5*sqrt(2)*KroneckerDelta(j,0)
    assert cg_simp(5*c-5*d+5*e) == 5*sqrt(3)*KroneckerDelta(j,0)
    assert cg_simp(-a+b) == -sqrt(2)*KroneckerDelta(j,0)
    assert cg_simp(-c+d-e) == -sqrt(3)*KroneckerDelta(j,0)
    assert cg_simp(-6*a+6*b) == -6*sqrt(2)*KroneckerDelta(j,0)
    assert cg_simp(-4*c+4*d-4*e) == -4*sqrt(3)*KroneckerDelta(j,0)
    # Test Varshalovich 8.7.2 Eq 9
    a = CG(S(1)/2,S(1)/2,S(1)/2,-S(1)/2,1,0)*CG(S(1)/2,S(1)/2,S(1)/2,-S(1)/2,1,0)
    b = CG(S(1)/2,S(1)/2,S(1)/2,-S(1)/2,0,0)*CG(S(1)/2,S(1)/2,S(1)/2,-S(1)/2,0,0)
    c = CG(1,0,1,1,1,1)*CG(1,0,1,1,1,1)
    d = CG(1,0,1,1,2,1)*CG(1,0,1,1,2,1)
    assert cg_simp(a+b) == 1
    assert cg_simp(c+d) == 1
    assert cg_simp(a+b+c+d) == 2
    assert cg_simp(4*a+4*b) == 4
    assert cg_simp(4*c+4*d) == 4
    assert cg_simp(5*a+3*b) == 3+2*a
    assert cg_simp(5*c+3*d) == 3+2*c
    assert cg_simp(-a-b) == -1
    assert cg_simp(-c-d) == -1
    """a = CG(S(1)/2,m1,S(1)/2,m2,1,1)*CG(S(1)/2,m1p,S(1)/2,m2p,1,1)
    b = CG(S(1)/2,m1,S(1)/2,m2,1,0)*CG(S(1)/2,m1p,S(1)/2,m2p,1,0)
    c = CG(S(1)/2,m1,S(1)/2,m2,1,-1)*CG(S(1)/2,m1p,S(1)/2,m2p,1,-1)
    d = CG(S(1)/2,m1,S(1)/2,m2,0,0)*CG(S(1)/2,m1p,S(1)/2,m2p,0,0)
    assert cg_simp(a+b+c+d) == KroneckerDelta(m1,m1p)*KroneckerDelta(m2,m2p)"""

def test_doit():
    assert Wigner3j(1/2,-1/2,1/2,1/2,0,0).doit() == -sqrt(2)/2
    assert CG(1/2,1/2,1/2,-1/2,1,0).doit() == sqrt(2)/2
