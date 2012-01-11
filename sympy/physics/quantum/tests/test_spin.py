from __future__ import division
from sympy import cos, exp, I, Matrix, pi, S, sin, sqrt, Sum, symbols
from sympy.abc import alpha, beta, gamma, j, m

from sympy.physics.quantum import hbar, represent, Commutator, InnerProduct
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.cg import CG
from sympy.physics.quantum.spin import (
    Jx, Jy, Jz, Jplus, Jminus, J2,
    JxBra, JyBra, JzBra,
    JxKet, JyKet, JzKet,
    JxKetCoupled, JyKetCoupled, JzKetCoupled,
    couple, uncouple,
    Rotation, WignerD
)

from sympy.utilities.pytest import XFAIL

j1, j2, m1, m2, mi, mp = symbols('j1 j2 m1 m2 mi mp')

def test_represent():
    # Spin operators
    assert represent(Jx) == hbar*Matrix([[0,1],[1,0]])/2
    assert represent(Jx, j=1) == hbar*sqrt(2)*Matrix([[0,1,0],[1,0,1],[0,1,0]])/2
    assert represent(Jy) == hbar*I*Matrix([[0,-1],[1,0]])/2
    assert represent(Jy, j=1) == hbar*I*sqrt(2)*Matrix([[0,-1,0],[1,0,-1],[0,1,0]])/2
    assert represent(Jz) == hbar*Matrix([[1,0],[0,-1]])/2
    assert represent(Jz, j=1) == hbar*Matrix([[1,0,0],[0,0,0],[0,0,-1]])
    # Spin states
    # Jx basis
    assert represent(JxKet(S(1)/2,S(1)/2), basis=Jx) == Matrix([1,0])
    assert represent(JxKet(S(1)/2,-S(1)/2), basis=Jx) == Matrix([0,1])
    assert represent(JxKet(1,1), basis=Jx) == Matrix([1,0,0])
    assert represent(JxKet(1,0), basis=Jx) == Matrix([0,1,0])
    assert represent(JxKet(1,-1), basis=Jx) == Matrix([0,0,1])
    assert represent(JyKet(S(1)/2,S(1)/2), basis=Jx) == Matrix([exp(-I*pi/4),0])
    assert represent(JyKet(S(1)/2,-S(1)/2), basis=Jx) == Matrix([0,exp(I*pi/4)])
    assert represent(JyKet(1,1), basis=Jx) == Matrix([-I,0,0])
    assert represent(JyKet(1,0), basis=Jx) == Matrix([0,1,0])
    assert represent(JyKet(1,-1), basis=Jx) == Matrix([0,0,I])
    assert represent(JzKet(S(1)/2,S(1)/2), basis=Jx) == sqrt(2)*Matrix([-1,1])/2
    assert represent(JzKet(S(1)/2,-S(1)/2), basis=Jx) == sqrt(2)*Matrix([-1,-1])/2
    assert represent(JzKet(1,1), basis=Jx) == Matrix([1,-sqrt(2),1])/2
    assert represent(JzKet(1,0), basis=Jx) == sqrt(2)*Matrix([1,0,-1])/2
    assert represent(JzKet(1,-1), basis=Jx) == Matrix([1,sqrt(2),1])/2
    # Jy basis
    assert represent(JxKet(S(1)/2,S(1)/2), basis=Jy) == Matrix([exp(-3*I*pi/4),0])
    assert represent(JxKet(S(1)/2,-S(1)/2), basis=Jy) == Matrix([0,exp(3*I*pi/4)])
    assert represent(JxKet(1,1), basis=Jy) == Matrix([I,0,0])
    assert represent(JxKet(1,0), basis=Jy) == Matrix([0,1,0])
    assert represent(JxKet(1,-1), basis=Jy) == Matrix([0,0,-I])
    assert represent(JyKet(S(1)/2,S(1)/2), basis=Jy) == Matrix([1,0])
    assert represent(JyKet(S(1)/2,-S(1)/2), basis=Jy) == Matrix([0,1])
    assert represent(JyKet(1,1), basis=Jy) == Matrix([1,0,0])
    assert represent(JyKet(1,0), basis=Jy) == Matrix([0,1,0])
    assert represent(JyKet(1,-1), basis=Jy) == Matrix([0,0,1])
    assert represent(JzKet(S(1)/2,S(1)/2), basis=Jy) == sqrt(2)*Matrix([-1,I])/2
    assert represent(JzKet(S(1)/2,-S(1)/2), basis=Jy) == sqrt(2)*Matrix([I,-1])/2
    assert represent(JzKet(1,1), basis=Jy) == Matrix([1,-I*sqrt(2),-1])/2
    assert represent(JzKet(1,0), basis=Jy) == Matrix([-sqrt(2)*I,0,-sqrt(2)*I])/2
    assert represent(JzKet(1,-1), basis=Jy) == Matrix([-1,-sqrt(2)*I,1])/2
    # Jz basis
    assert represent(JxKet(S(1)/2,S(1)/2), basis=Jz) == sqrt(2)*Matrix([1,1])/2
    assert represent(JxKet(S(1)/2,-S(1)/2), basis=Jz) == sqrt(2)*Matrix([-1,1])/2
    assert represent(JxKet(1,1), basis=Jz) == Matrix([1,sqrt(2),1])/2
    assert represent(JxKet(1,0), basis=Jz) == sqrt(2)*Matrix([-1,0,1])/2
    assert represent(JxKet(1,-1), basis=Jz) == Matrix([1,-sqrt(2),1])/2
    assert represent(JyKet(S(1)/2,S(1)/2), basis=Jz) == sqrt(2)*Matrix([-1,-I])/2
    assert represent(JyKet(S(1)/2,-S(1)/2), basis=Jz) == sqrt(2)*Matrix([-I,-1])/2
    assert represent(JyKet(1,1), basis=Jz) == Matrix([1,sqrt(2)*I,-1])/2
    assert represent(JyKet(1,0), basis=Jz) == sqrt(2)*Matrix([I,0,I])/2
    assert represent(JyKet(1,-1), basis=Jz) == Matrix([-1,sqrt(2)*I,1])/2
    assert represent(JzKet(S(1)/2,S(1)/2), basis=Jz) == Matrix([1,0])
    assert represent(JzKet(S(1)/2,-S(1)/2), basis=Jz) == Matrix([0,1])
    assert represent(JzKet(1,1), basis=Jz) == Matrix([1,0,0])
    assert represent(JzKet(1,0), basis=Jz) == Matrix([0,1,0])
    assert represent(JzKet(1,-1), basis=Jz) == Matrix([0,0,1])
    # Uncoupled states
    # Jx basis
    assert represent(TensorProduct(JxKet(S(1)/2,S(1)/2),JxKet(S(1)/2,S(1)/2)), basis=Jx) == \
        Matrix([1,0,0,0])
    assert represent(TensorProduct(JxKet(S(1)/2,S(1)/2),JxKet(S(1)/2,-S(1)/2)), basis=Jx) == \
        Matrix([0,1,0,0])
    assert represent(TensorProduct(JxKet(S(1)/2,-S(1)/2),JxKet(S(1)/2,S(1)/2)), basis=Jx) == \
        Matrix([0,0,1,0])
    assert represent(TensorProduct(JxKet(S(1)/2,-S(1)/2),JxKet(S(1)/2,-S(1)/2)), basis=Jx) == \
        Matrix([0,0,0,1])
    assert represent(TensorProduct(JyKet(S(1)/2,S(1)/2),JyKet(S(1)/2,S(1)/2)), basis=Jx) == \
        Matrix([-I,0,0,0])
    assert represent(TensorProduct(JyKet(S(1)/2,S(1)/2),JyKet(S(1)/2,-S(1)/2)), basis=Jx) == \
        Matrix([0,1,0,0])
    assert represent(TensorProduct(JyKet(S(1)/2,-S(1)/2),JyKet(S(1)/2,S(1)/2)), basis=Jx) == \
        Matrix([0,0,1,0])
    assert represent(TensorProduct(JyKet(S(1)/2,-S(1)/2),JyKet(S(1)/2,-S(1)/2)), basis=Jx) == \
        Matrix([0,0,0,I])
    assert represent(TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,S(1)/2)), basis=Jx) == \
        Matrix([S(1)/2,-S(1)/2,-S(1)/2,S(1)/2])
    assert represent(TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,-S(1)/2)), basis=Jx) == \
        Matrix([S(1)/2,S(1)/2,-S(1)/2,-S(1)/2])
    assert represent(TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,S(1)/2)), basis=Jx) == \
        Matrix([S(1)/2,-S(1)/2,S(1)/2,-S(1)/2])
    assert represent(TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,-S(1)/2)), basis=Jx) == \
        Matrix([S(1)/2,S(1)/2,S(1)/2,S(1)/2])
    # Jy basis
    assert represent(TensorProduct(JxKet(S(1)/2,S(1)/2),JxKet(S(1)/2,S(1)/2)), basis=Jy) == \
        Matrix([I,0,0,0])
    assert represent(TensorProduct(JxKet(S(1)/2,S(1)/2),JxKet(S(1)/2,-S(1)/2)), basis=Jy) == \
        Matrix([0,1,0,0])
    assert represent(TensorProduct(JxKet(S(1)/2,-S(1)/2),JxKet(S(1)/2,S(1)/2)), basis=Jy) == \
        Matrix([0,0,1,0])
    assert represent(TensorProduct(JxKet(S(1)/2,-S(1)/2),JxKet(S(1)/2,-S(1)/2)), basis=Jy) == \
        Matrix([0,0,0,-I])
    assert represent(TensorProduct(JyKet(S(1)/2,S(1)/2),JyKet(S(1)/2,S(1)/2)), basis=Jy) == \
        Matrix([1,0,0,0])
    assert represent(TensorProduct(JyKet(S(1)/2,S(1)/2),JyKet(S(1)/2,-S(1)/2)), basis=Jy) == \
        Matrix([0,1,0,0])
    assert represent(TensorProduct(JyKet(S(1)/2,-S(1)/2),JyKet(S(1)/2,S(1)/2)), basis=Jy) == \
        Matrix([0,0,1,0])
    assert represent(TensorProduct(JyKet(S(1)/2,-S(1)/2),JyKet(S(1)/2,-S(1)/2)), basis=Jy) == \
        Matrix([0,0,0,1])
    assert represent(TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,S(1)/2)), basis=Jy) == \
        Matrix([S(1)/2,-I/2,-I/2,-S(1)/2])
    assert represent(TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,-S(1)/2)), basis=Jy) == \
        Matrix([-I/2,S(1)/2,-S(1)/2,-I/2])
    assert represent(TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,S(1)/2)), basis=Jy) == \
        Matrix([-I/2,-S(1)/2,S(1)/2,-I/2])
    assert represent(TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,-S(1)/2)), basis=Jy) == \
        Matrix([-S(1)/2,-I/2,-I/2,S(1)/2])
    # Jz basis
    assert represent(TensorProduct(JxKet(S(1)/2,S(1)/2),JxKet(S(1)/2,S(1)/2)), basis=Jz) == \
        Matrix([S(1)/2,S(1)/2,S(1)/2,S(1)/2])
    assert represent(TensorProduct(JxKet(S(1)/2,S(1)/2),JxKet(S(1)/2,-S(1)/2)), basis=Jz) == \
        Matrix([-S(1)/2,S(1)/2,-S(1)/2,S(1)/2])
    assert represent(TensorProduct(JxKet(S(1)/2,-S(1)/2),JxKet(S(1)/2,S(1)/2)), basis=Jz) == \
        Matrix([-S(1)/2,-S(1)/2,S(1)/2,S(1)/2])
    assert represent(TensorProduct(JxKet(S(1)/2,-S(1)/2),JxKet(S(1)/2,-S(1)/2)), basis=Jz) == \
        Matrix([S(1)/2,-S(1)/2,-S(1)/2,S(1)/2])
    assert represent(TensorProduct(JyKet(S(1)/2,S(1)/2),JyKet(S(1)/2,S(1)/2)), basis=Jz) == \
        Matrix([S(1)/2,I/2,I/2,-S(1)/2])
    assert represent(TensorProduct(JyKet(S(1)/2,S(1)/2),JyKet(S(1)/2,-S(1)/2)), basis=Jz) == \
        Matrix([I/2,S(1)/2,-S(1)/2,I/2])
    assert represent(TensorProduct(JyKet(S(1)/2,-S(1)/2),JyKet(S(1)/2,S(1)/2)), basis=Jz) == \
        Matrix([I/2,-S(1)/2,S(1)/2,I/2])
    assert represent(TensorProduct(JyKet(S(1)/2,-S(1)/2),JyKet(S(1)/2,-S(1)/2)), basis=Jz) == \
        Matrix([-S(1)/2,I/2,I/2,S(1)/2])
    assert represent(TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,S(1)/2)), basis=Jz) == \
        Matrix([1,0,0,0])
    assert represent(TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,-S(1)/2)), basis=Jz) == \
        Matrix([0,1,0,0])
    assert represent(TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,S(1)/2)), basis=Jz) == \
        Matrix([0,0,1,0])
    assert represent(TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,-S(1)/2)), basis=Jz) == \
        Matrix([0,0,0,1])
    # Coupled spin states
    # Jx basis
    assert represent(JxKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([1,0,0,0])
    assert represent(JxKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,1,0,0])
    assert represent(JxKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,0,1,0])
    assert represent(JxKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,0,0,1])
    assert represent(JyKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([1,0,0,0])
    assert represent(JyKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,-I,0,0])
    assert represent(JyKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,0,1,0])
    assert represent(JyKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,0,0,I])
    assert represent(JzKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([1,0,0,0])
    assert represent(JzKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,S(1)/2,-sqrt(2)/2,S(1)/2])
    assert represent(JzKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,sqrt(2)/2,0,-sqrt(2)/2])
    assert represent(JzKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jx) == \
        Matrix([0,S(1)/2,sqrt(2)/2,S(1)/2])
    # Jy basis
    assert represent(JxKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([1,0,0,0])
    assert represent(JxKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,I,0,0])
    assert represent(JxKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,0,1,0])
    assert represent(JxKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,0,0,-I])
    assert represent(JyKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([1,0,0,0])
    assert represent(JyKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,1,0,0])
    assert represent(JyKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,0,1,0])
    assert represent(JyKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,0,0,1])
    assert represent(JzKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([1,0,0,0])
    assert represent(JzKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,S(1)/2,-I*sqrt(2)/2,-S(1)/2])
    assert represent(JzKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,-I*sqrt(2)/2,0,-I*sqrt(2)/2])
    assert represent(JzKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jy) == \
        Matrix([0,-S(1)/2,-I*sqrt(2)/2,S(1)/2])
    # Jz basis
    assert represent(JxKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([1,0,0,0])
    assert represent(JxKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,S(1)/2,sqrt(2)/2,S(1)/2])
    assert represent(JxKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,-sqrt(2)/2,0,sqrt(2)/2])
    assert represent(JxKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,S(1)/2,-sqrt(2)/2,S(1)/2])
    assert represent(JyKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([1,0,0,0])
    assert represent(JyKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,S(1)/2,I*sqrt(2)/2,-S(1)/2])
    assert represent(JyKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,I*sqrt(2)/2,0,I*sqrt(2)/2])
    assert represent(JyKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,-S(1)/2,I*sqrt(2)/2,S(1)/2])
    assert represent(JzKetCoupled(0,0,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([1,0,0,0])
    assert represent(JzKetCoupled(1,1,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,1,0,0])
    assert represent(JzKetCoupled(1,0,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,0,1,0])
    assert represent(JzKetCoupled(1,-1,S(1)/2,S(1)/2), basis=Jz) == \
        Matrix([0,0,0,1])

def test_rewrite():
    # Rewrite to same basis
    assert JxBra(1,1).rewrite('Jx') == JxBra(1,1)
    assert JxBra(j,m).rewrite('Jx') == JxBra(j,m)
    assert JxKet(1,1).rewrite('Jx') == JxKet(1,1)
    assert JxKet(j,m).rewrite('Jx') == JxKet(j,m)
    # Rewriting a normal state
    # Bra
    # Numerical
    assert JxBra(1,1).rewrite('Jy') == -I*JyBra(1,1)
    assert JxBra(1,0).rewrite('Jy') == JyBra(1,0)
    assert JxBra(1,-1).rewrite('Jy') == I*JyBra(1,-1)
    assert JxBra(1,1).rewrite('Jz') == JzBra(1,1)/2+JzBra(1,0)/sqrt(2)+JzBra(1,-1)/2
    assert JxBra(1,0).rewrite('Jz') == -sqrt(2)*JzBra(1,1)/2+sqrt(2)*JzBra(1,-1)/2
    assert JxBra(1,-1).rewrite('Jz') == JzBra(1,1)/2-JzBra(1,0)/sqrt(2)+JzBra(1,-1)/2
    assert JyBra(1,1).rewrite('Jx') == I*JxBra(1,1)
    assert JyBra(1,0).rewrite('Jx') == JxBra(1,0)
    assert JyBra(1,-1).rewrite('Jx') == -I*JxBra(1,-1)
    assert JyBra(1,1).rewrite('Jz') == JzBra(1,1)/2-sqrt(2)*I*JzBra(1,0)/2-JzBra(1,-1)/2
    assert JyBra(1,0).rewrite('Jz') == -sqrt(2)*I*JzBra(1,1)/2-sqrt(2)*I*JzBra(1,-1)/2
    assert JyBra(1,-1).rewrite('Jz') == -JzBra(1,1)/2-sqrt(2)*I*JzBra(1,0)/2+JzBra(1,-1)/2
    assert JzBra(1,1).rewrite('Jx') == JxBra(1,1)/2-sqrt(2)*JxBra(1,0)/2+JxBra(1,-1)/2
    assert JzBra(1,0).rewrite('Jx') == sqrt(2)*JxBra(1,1)/2-sqrt(2)*JxBra(1,-1)/2
    assert JzBra(1,-1).rewrite('Jx') == JxBra(1,1)/2+sqrt(2)*JxBra(1,0)/2+JxBra(1,-1)/2
    assert JzBra(1,1).rewrite('Jy') == JyBra(1,1)/2+sqrt(2)*I*JyBra(1,0)/2-JyBra(1,-1)/2
    assert JzBra(1,0).rewrite('Jy') == sqrt(2)*I*JyBra(1,1)/2+sqrt(2)*I*JyBra(1,-1)/2
    assert JzBra(1,-1).rewrite('Jy') == -JyBra(1,1)/2+sqrt(2)*I*JyBra(1,0)/2+JyBra(1,-1)/2
    # Symbolic
    assert JxBra(j,m).rewrite('Jy') == Sum(WignerD(j,mi,m,3*pi/2,0,0) * JyBra(j,mi), (mi, -j, j))
    assert JxBra(j,m).rewrite('Jz') == Sum(WignerD(j,mi,m,0,pi/2,0) * JzBra(j,mi), (mi, -j, j))
    assert JyBra(j,m).rewrite('Jx') == Sum(WignerD(j,mi,m,0,0,pi/2) * JxBra(j,mi), (mi, -j, j))
    assert JyBra(j,m).rewrite('Jz') == Sum(WignerD(j,mi,m,3*pi/2,-pi/2,pi/2) * JzBra(j,mi), (mi, -j, j))
    assert JzBra(j,m).rewrite('Jx') == Sum(WignerD(j,mi,m,0,3*pi/2,0) * JxBra(j,mi), (mi, -j, j))
    assert JzBra(j,m).rewrite('Jy') == Sum(WignerD(j,mi,m,3*pi/2,pi/2,pi/2) * JyBra(j,mi), (mi, -j, j))
    # Kets
    # Numerical
    assert JxKet(1,1).rewrite('Jy') == I*JyKet(1,1)
    assert JxKet(1,0).rewrite('Jy') == JyKet(1,0)
    assert JxKet(1,-1).rewrite('Jy') == -I*JyKet(1,-1)
    assert JxKet(1,1).rewrite('Jz') == JzKet(1,1)/2+JzKet(1,0)/sqrt(2)+JzKet(1,-1)/2
    assert JxKet(1,0).rewrite('Jz') == -sqrt(2)*JzKet(1,1)/2+sqrt(2)*JzKet(1,-1)/2
    assert JxKet(1,-1).rewrite('Jz') == JzKet(1,1)/2-JzKet(1,0)/sqrt(2)+JzKet(1,-1)/2
    assert JyKet(1,1).rewrite('Jx') == -I*JxKet(1,1)
    assert JyKet(1,0).rewrite('Jx') == JxKet(1,0)
    assert JyKet(1,-1).rewrite('Jx') == I*JxKet(1,-1)
    assert JyKet(1,1).rewrite('Jz') == JzKet(1,1)/2+sqrt(2)*I*JzKet(1,0)/2-JzKet(1,-1)/2
    assert JyKet(1,0).rewrite('Jz') == sqrt(2)*I*JzKet(1,1)/2+sqrt(2)*I*JzKet(1,-1)/2
    assert JyKet(1,-1).rewrite('Jz') == -JzKet(1,1)/2+sqrt(2)*I*JzKet(1,0)/2+JzKet(1,-1)/2
    assert JzKet(1,1).rewrite('Jx') == JxKet(1,1)/2-sqrt(2)*JxKet(1,0)/2+JxKet(1,-1)/2
    assert JzKet(1,0).rewrite('Jx') == sqrt(2)*JxKet(1,1)/2-sqrt(2)*JxKet(1,-1)/2
    assert JzKet(1,-1).rewrite('Jx') == JxKet(1,1)/2+sqrt(2)*JxKet(1,0)/2+JxKet(1,-1)/2
    assert JzKet(1,1).rewrite('Jy') == JyKet(1,1)/2-sqrt(2)*I*JyKet(1,0)/2-JyKet(1,-1)/2
    assert JzKet(1,0).rewrite('Jy') == -sqrt(2)*I*JyKet(1,1)/2-sqrt(2)*I*JyKet(1,-1)/2
    assert JzKet(1,-1).rewrite('Jy') == -JyKet(1,1)/2-sqrt(2)*I*JyKet(1,0)/2+JyKet(1,-1)/2
    # Symbolic
    assert JxKet(j,m).rewrite('Jy') == Sum(WignerD(j,mi,m,3*pi/2,0,0) * JyKet(j,mi), (mi, -j, j))
    assert JxKet(j,m).rewrite('Jz') == Sum(WignerD(j,mi,m,0,pi/2,0) * JzKet(j,mi), (mi, -j, j))
    assert JyKet(j,m).rewrite('Jx') == Sum(WignerD(j,mi,m,0,0,pi/2) * JxKet(j,mi), (mi, -j, j))
    assert JyKet(j,m).rewrite('Jz') == Sum(WignerD(j,mi,m,3*pi/2,-pi/2,pi/2) * JzKet(j,mi), (mi, -j, j))
    assert JzKet(j,m).rewrite('Jx') == Sum(WignerD(j,mi,m,0,3*pi/2,0) * JxKet(j,mi), (mi, -j, j))
    assert JzKet(j,m).rewrite('Jy') == Sum(WignerD(j,mi,m,3*pi/2,pi/2,pi/2) * JyKet(j,mi), (mi, -j, j))
    # Rewrite an uncoupled state
    # Numerical
    assert TensorProduct(JyKet(1,1),JxKet(1,1)).rewrite('Jx') == -I*TensorProduct(JxKet(1,1),JxKet(1,1))
    assert TensorProduct(JyKet(1,0),JxKet(1,1)).rewrite('Jx') == TensorProduct(JxKet(1,0),JxKet(1,1))
    assert TensorProduct(JyKet(1,-1),JxKet(1,1)).rewrite('Jx') == I*TensorProduct(JxKet(1,-1),JxKet(1,1))
    assert TensorProduct(JzKet(1,1),JxKet(1,1)).rewrite('Jx') == \
            TensorProduct(JxKet(1,-1),JxKet(1,1))/2-sqrt(2)*TensorProduct(JxKet(1,0),JxKet(1,1))/2+TensorProduct(JxKet(1,1),JxKet(1,1))/2
    assert TensorProduct(JzKet(1,0),JxKet(1,1)).rewrite('Jx') == \
            -sqrt(2)*TensorProduct(JxKet(1,-1),JxKet(1,1))/2+sqrt(2)*TensorProduct(JxKet(1,1),JxKet(1,1))/2
    assert TensorProduct(JzKet(1,-1),JxKet(1,1)).rewrite('Jx') == \
            TensorProduct(JxKet(1,-1),JxKet(1,1))/2 + sqrt(2)*TensorProduct(JxKet(1,0),JxKet(1,1))/2 + TensorProduct(JxKet(1,1),JxKet(1,1))/2
    assert TensorProduct(JxKet(1,1),JyKet(1,1)).rewrite('Jy') == I*TensorProduct(JyKet(1,1),JyKet(1,1))
    assert TensorProduct(JxKet(1,0),JyKet(1,1)).rewrite('Jy') == TensorProduct(JyKet(1,0),JyKet(1,1))
    assert TensorProduct(JxKet(1,-1),JyKet(1,1)).rewrite('Jy') == -I*TensorProduct(JyKet(1,-1),JyKet(1,1))
    assert TensorProduct(JzKet(1,1),JyKet(1,1)).rewrite('Jy') == \
            -TensorProduct(JyKet(1,-1),JyKet(1,1))/2 - sqrt(2)*I*TensorProduct(JyKet(1,0),JyKet(1,1))/2 + TensorProduct(JyKet(1,1),JyKet(1,1))/2
    assert TensorProduct(JzKet(1,0),JyKet(1,1)).rewrite('Jy') == \
            -sqrt(2)*I*TensorProduct(JyKet(1,-1),JyKet(1,1))/2 - sqrt(2)*I*TensorProduct(JyKet(1,1),JyKet(1,1))/2
    assert TensorProduct(JzKet(1,-1),JyKet(1,1)).rewrite('Jy') == \
            TensorProduct(JyKet(1,-1),JyKet(1,1))/2 - sqrt(2)*I*TensorProduct(JyKet(1,0),JyKet(1,1))/2 - TensorProduct(JyKet(1,1),JyKet(1,1))/2
    assert TensorProduct(JxKet(1,1),JzKet(1,1)).rewrite('Jz') == \
            TensorProduct(JzKet(1,-1),JzKet(1,1))/2 + sqrt(2)*TensorProduct(JzKet(1,0),JzKet(1,1))/2 + TensorProduct(JzKet(1,1),JzKet(1,1))/2
    assert TensorProduct(JxKet(1,0),JzKet(1,1)).rewrite('Jz') == \
            sqrt(2)*TensorProduct(JzKet(1,-1),JzKet(1,1))/2 - sqrt(2)*TensorProduct(JzKet(1,1),JzKet(1,1))/2
    assert TensorProduct(JxKet(1,-1),JzKet(1,1)).rewrite('Jz') == \
            TensorProduct(JzKet(1,-1),JzKet(1,1))/2 - sqrt(2)*TensorProduct(JzKet(1,0),JzKet(1,1))/2 + TensorProduct(JzKet(1,1),JzKet(1,1))/2
    assert TensorProduct(JyKet(1,1),JzKet(1,1)).rewrite('Jz') == \
            -TensorProduct(JzKet(1,-1),JzKet(1,1))/2 + sqrt(2)*I*TensorProduct(JzKet(1,0),JzKet(1,1))/2 + TensorProduct(JzKet(1,1),JzKet(1,1))/2
    assert TensorProduct(JyKet(1,0),JzKet(1,1)).rewrite('Jz') == \
            sqrt(2)*I*TensorProduct(JzKet(1,-1),JzKet(1,1))/2 + sqrt(2)*I*TensorProduct(JzKet(1,1),JzKet(1,1))/2
    assert TensorProduct(JyKet(1,-1),JzKet(1,1)).rewrite('Jz') == \
            TensorProduct(JzKet(1,-1),JzKet(1,1))/2 + sqrt(2)*I*TensorProduct(JzKet(1,0),JzKet(1,1))/2 - TensorProduct(JzKet(1,1),JzKet(1,1))/2
    # Symbolic
    assert TensorProduct(JyKet(j1,m1), JxKet(j2,m2)).rewrite('Jy') == \
            TensorProduct(JyKet(j1,m1), Sum(WignerD(j2,mi,m2,3*pi/2,0,0) * JyKet(j2,mi), (mi, -j2, j2)))
    assert TensorProduct(JzKet(j1,m1), JxKet(j2,m2)).rewrite('Jz') == \
            TensorProduct(JzKet(j1,m1), Sum(WignerD(j2,mi,m2,0,pi/2,0) * JzKet(j2,mi), (mi, -j2, j2)))
    assert TensorProduct(JxKet(j1,m1), JyKet(j2,m2)).rewrite('Jx') == \
            TensorProduct(JxKet(j1,m1), Sum(WignerD(j2,mi,m2,0,0,pi/2) * JxKet(j2,mi), (mi, -j2, j2)))
    assert TensorProduct(JzKet(j1,m1), JyKet(j2,m2)).rewrite('Jz') == \
            TensorProduct(JzKet(j1,m1), Sum(WignerD(j2,mi,m2,3*pi/2,-pi/2,pi/2) * JzKet(j2,mi), (mi, -j2, j2)))
    assert TensorProduct(JxKet(j1,m1), JzKet(j2,m2)).rewrite('Jx') == \
            TensorProduct(JxKet(j1,m1), Sum(WignerD(j2,mi,m2,0,3*pi/2,0) * JxKet(j2,mi), (mi, -j2, j2)))
    assert TensorProduct(JyKet(j1,m1), JzKet(j2,m2)).rewrite('Jy') == \
            TensorProduct(JyKet(j1,m1), Sum(WignerD(j2,mi,m2,3*pi/2,pi/2,pi/2) * JyKet(j2,mi), (mi, -j2, j2)))
    # Rewrite a coupled state
    # Numerical
    assert JyKetCoupled(0,0,S(1)/2,S(1)/2).rewrite('Jx') == \
        JxKetCoupled(0,0,S(1)/2,S(1)/2)
    assert JyKetCoupled(1,1,S(1)/2,S(1)/2).rewrite('Jx') == \
        -I*JxKetCoupled(1,1,S(1)/2,S(1)/2)
    assert JyKetCoupled(1,0,S(1)/2,S(1)/2).rewrite('Jx') == \
        JxKetCoupled(1,0,S(1)/2,S(1)/2)
    assert JyKetCoupled(1,-1,S(1)/2,S(1)/2).rewrite('Jx') == \
        I*JxKetCoupled(1,-1,S(1)/2,S(1)/2)
    assert JzKetCoupled(0,0,S(1)/2,S(1)/2).rewrite('Jx') == \
        JxKetCoupled(0,0,S(1)/2,S(1)/2)
    assert JzKetCoupled(1,1,S(1)/2,S(1)/2).rewrite('Jx') == \
        JxKetCoupled(1,1,S(1)/2,S(1)/2)/2 - sqrt(2)*JxKetCoupled(1,0,S(1)/2,S(1)/2)/2 + JxKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JzKetCoupled(1,0,S(1)/2,S(1)/2).rewrite('Jx') == \
        sqrt(2)*JxKetCoupled(1,1,S(1)/2,S(1)/2)/2 - sqrt(2)*JxKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JzKetCoupled(1,-1,S(1)/2,S(1)/2).rewrite('Jx') == \
        JxKetCoupled(1,1,S(1)/2,S(1)/2)/2 + sqrt(2)*JxKetCoupled(1,0,S(1)/2,S(1)/2)/2 + JxKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JxKetCoupled(0,0,S(1)/2,S(1)/2).rewrite('Jy') == \
        JyKetCoupled(0,0,S(1)/2,S(1)/2)
    assert JxKetCoupled(1,1,S(1)/2,S(1)/2).rewrite('Jy') == \
        I*JyKetCoupled(1,1,S(1)/2,S(1)/2)
    assert JxKetCoupled(1,0,S(1)/2,S(1)/2).rewrite('Jy') == \
        JyKetCoupled(1,0,S(1)/2,S(1)/2)
    assert JxKetCoupled(1,-1,S(1)/2,S(1)/2).rewrite('Jy') == \
        -I*JyKetCoupled(1,-1,S(1)/2,S(1)/2)
    assert JzKetCoupled(0,0,S(1)/2,S(1)/2).rewrite('Jy') == \
        JyKetCoupled(0,0,S(1)/2,S(1)/2)
    assert JzKetCoupled(1,1,S(1)/2,S(1)/2).rewrite('Jy') == \
        JyKetCoupled(1,1,S(1)/2,S(1)/2)/2 - I*sqrt(2)*JyKetCoupled(1,0,S(1)/2,S(1)/2)/2 - JyKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JzKetCoupled(1,0,S(1)/2,S(1)/2).rewrite('Jy') == \
        -I*sqrt(2)*JyKetCoupled(1,1,S(1)/2,S(1)/2)/2 - I*sqrt(2)*JyKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JzKetCoupled(1,-1,S(1)/2,S(1)/2).rewrite('Jy') == \
        -JyKetCoupled(1,1,S(1)/2,S(1)/2)/2 - I*sqrt(2)*JyKetCoupled(1,0,S(1)/2,S(1)/2)/2 + JyKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JxKetCoupled(0,0,S(1)/2,S(1)/2).rewrite('Jz') == \
        JzKetCoupled(0,0,S(1)/2,S(1)/2)
    assert JxKetCoupled(1,1,S(1)/2,S(1)/2).rewrite('Jz') == \
        JzKetCoupled(1,1,S(1)/2,S(1)/2)/2 + sqrt(2)*JzKetCoupled(1,0,S(1)/2,S(1)/2)/2 + JzKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JxKetCoupled(1,0,S(1)/2,S(1)/2).rewrite('Jz') == \
        -sqrt(2)*JzKetCoupled(1,1,S(1)/2,S(1)/2)/2 + sqrt(2)*JzKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JxKetCoupled(1,-1,S(1)/2,S(1)/2).rewrite('Jz') == \
        JzKetCoupled(1,1,S(1)/2,S(1)/2)/2 - sqrt(2)*JzKetCoupled(1,0,S(1)/2,S(1)/2)/2 + JzKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JyKetCoupled(0,0,S(1)/2,S(1)/2).rewrite('Jz') == \
        JzKetCoupled(0,0,S(1)/2,S(1)/2)
    assert JyKetCoupled(1,1,S(1)/2,S(1)/2).rewrite('Jz') == \
        JzKetCoupled(1,1,S(1)/2,S(1)/2)/2 + I*sqrt(2)*JzKetCoupled(1,0,S(1)/2,S(1)/2)/2 - JzKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JyKetCoupled(1,0,S(1)/2,S(1)/2).rewrite('Jz') == \
        I*sqrt(2)*JzKetCoupled(1,1,S(1)/2,S(1)/2)/2 + I*sqrt(2)*JzKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    assert JyKetCoupled(1,-1,S(1)/2,S(1)/2).rewrite('Jz') == \
        -JzKetCoupled(1,1,S(1)/2,S(1)/2)/2 + I*sqrt(2)*JzKetCoupled(1,0,S(1)/2,S(1)/2)/2 + JzKetCoupled(1,-1,S(1)/2,S(1)/2)/2
    # Symbolic
    assert JyKetCoupled(j,m,j1,j2).rewrite('Jx') == \
        Sum(WignerD(j,mi,m,0,0,pi/2) * JxKetCoupled(j,mi,j1,j2), (mi,-j,j))
    assert JzKetCoupled(j,m,j1,j2).rewrite('Jx') == \
        Sum(WignerD(j,mi,m,0,3*pi/2,0) * JxKetCoupled(j,mi,j1,j2), (mi,-j,j))
    assert JxKetCoupled(j,m,j1,j2).rewrite('Jy') == \
        Sum(WignerD(j,mi,m,3*pi/2,0,0) * JyKetCoupled(j,mi,j1,j2), (mi,-j,j))
    assert JzKetCoupled(j,m,j1,j2).rewrite('Jy') == \
        Sum(WignerD(j,mi,m,3*pi/2,pi/2,pi/2) * JyKetCoupled(j,mi,j1,j2), (mi,-j,j))
    assert JxKetCoupled(j,m,j1,j2).rewrite('Jz') == \
        Sum(WignerD(j,mi,m,0,pi/2,0) * JzKetCoupled(j,mi,j1,j2), (mi,-j,j))
    assert JyKetCoupled(j,m,j1,j2).rewrite('Jz') == \
        Sum(WignerD(j,mi,m,3*pi/2,-pi/2,pi/2) * JzKetCoupled(j,mi,j1,j2), (mi,-j,j))
    # Innerproducts of rewritten states
    # Numerical
    assert qapply(JxBra(1,1)*JxKet(1,1).rewrite('Jy')).doit() == 1
    assert qapply(JxBra(1,0)*JxKet(1,0).rewrite('Jy')).doit() == 1
    assert qapply(JxBra(1,-1)*JxKet(1,-1).rewrite('Jy')).doit() == 1
    assert qapply(JxBra(1,1)*JxKet(1,1).rewrite('Jz')).doit() == 1
    assert qapply(JxBra(1,0)*JxKet(1,0).rewrite('Jz')).doit() == 1
    assert qapply(JxBra(1,-1)*JxKet(1,-1).rewrite('Jz')).doit() == 1
    assert qapply(JyBra(1,1)*JyKet(1,1).rewrite('Jx')).doit() == 1
    assert qapply(JyBra(1,0)*JyKet(1,0).rewrite('Jx')).doit() == 1
    assert qapply(JyBra(1,-1)*JyKet(1,-1).rewrite('Jx')).doit() == 1
    assert qapply(JyBra(1,1)*JyKet(1,1).rewrite('Jz')).doit() == 1
    assert qapply(JyBra(1,0)*JyKet(1,0).rewrite('Jz')).doit() == 1
    assert qapply(JyBra(1,-1)*JyKet(1,-1).rewrite('Jz')).doit() == 1
    assert qapply(JyBra(1,1)*JyKet(1,1).rewrite('Jz')).doit() == 1
    assert qapply(JyBra(1,0)*JyKet(1,0).rewrite('Jz')).doit() == 1
    assert qapply(JyBra(1,-1)*JyKet(1,-1).rewrite('Jz')).doit() == 1
    assert qapply(JzBra(1,1)*JzKet(1,1).rewrite('Jy')).doit() == 1
    assert qapply(JzBra(1,0)*JzKet(1,0).rewrite('Jy')).doit() == 1
    assert qapply(JzBra(1,-1)*JzKet(1,-1).rewrite('Jy')).doit() == 1
    assert qapply(JxBra(1,1)*JxKet(1,0).rewrite('Jy')).doit() == 0
    assert qapply(JxBra(1,1)*JxKet(1,-1).rewrite('Jy')) == 0
    assert qapply(JxBra(1,1)*JxKet(1,0).rewrite('Jz')).doit() == 0
    assert qapply(JxBra(1,1)*JxKet(1,-1).rewrite('Jz')) == 0
    assert qapply(JyBra(1,1)*JyKet(1,0).rewrite('Jx')).doit() == 0
    assert qapply(JyBra(1,1)*JyKet(1,-1).rewrite('Jx')) == 0
    assert qapply(JyBra(1,1)*JyKet(1,0).rewrite('Jz')).doit() == 0
    assert qapply(JyBra(1,1)*JyKet(1,-1).rewrite('Jz')) == 0
    assert qapply(JzBra(1,1)*JzKet(1,0).rewrite('Jx')).doit() == 0
    assert qapply(JzBra(1,1)*JzKet(1,-1).rewrite('Jx')) == 0
    assert qapply(JzBra(1,1)*JzKet(1,0).rewrite('Jy')).doit() == 0
    assert qapply(JzBra(1,1)*JzKet(1,-1).rewrite('Jy')) == 0
    assert qapply(JxBra(1,0)*JxKet(1,1).rewrite('Jy')) == 0
    assert qapply(JxBra(1,0)*JxKet(1,-1).rewrite('Jy')) == 0
    assert qapply(JxBra(1,0)*JxKet(1,1).rewrite('Jz')) == 0
    assert qapply(JxBra(1,0)*JxKet(1,-1).rewrite('Jz')) == 0
    assert qapply(JyBra(1,0)*JyKet(1,1).rewrite('Jx')) == 0
    assert qapply(JyBra(1,0)*JyKet(1,-1).rewrite('Jx')) == 0
    assert qapply(JyBra(1,0)*JyKet(1,1).rewrite('Jz')) == 0
    assert qapply(JyBra(1,0)*JyKet(1,-1).rewrite('Jz')) == 0
    assert qapply(JzBra(1,0)*JzKet(1,1).rewrite('Jx')) == 0
    assert qapply(JzBra(1,0)*JzKet(1,-1).rewrite('Jx')) == 0
    assert qapply(JzBra(1,0)*JzKet(1,1).rewrite('Jy')) == 0
    assert qapply(JzBra(1,0)*JzKet(1,-1).rewrite('Jy')) == 0
    assert qapply(JxBra(1,-1)*JxKet(1,1).rewrite('Jy')) == 0
    assert qapply(JxBra(1,-1)*JxKet(1,0).rewrite('Jy')).doit() == 0
    assert qapply(JxBra(1,-1)*JxKet(1,1).rewrite('Jz')) == 0
    assert qapply(JxBra(1,-1)*JxKet(1,0).rewrite('Jz')).doit() == 0
    assert qapply(JyBra(1,-1)*JyKet(1,1).rewrite('Jx')) == 0
    assert qapply(JyBra(1,-1)*JyKet(1,0).rewrite('Jx')).doit() == 0
    assert qapply(JyBra(1,-1)*JyKet(1,1).rewrite('Jz')) == 0
    assert qapply(JyBra(1,-1)*JyKet(1,0).rewrite('Jz')).doit() == 0
    assert qapply(JzBra(1,-1)*JzKet(1,1).rewrite('Jx')) == 0
    assert qapply(JzBra(1,-1)*JzKet(1,0).rewrite('Jx')).doit() == 0
    assert qapply(JzBra(1,-1)*JzKet(1,1).rewrite('Jy')) == 0
    assert qapply(JzBra(1,-1)*JzKet(1,0).rewrite('Jy')).doit() == 0

def test_uncouple():
    # Numerical
    assert uncouple(JzKetCoupled(0,0,S(1)/2,S(1)/2)) == \
        sqrt(2)*TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,-S(1)/2))/2 - sqrt(2)*TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,S(1)/2))/2
    assert uncouple(JzKetCoupled(1,1,S(1)/2,S(1)/2)) == \
        TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,S(1)/2))
    assert uncouple(JzKetCoupled(1,0,S(1)/2,S(1)/2)) == \
        sqrt(2)*TensorProduct(JzKet(S(1)/2,S(1)/2),JzKet(S(1)/2,-S(1)/2))/2 + sqrt(2)*TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,S(1)/2))/2
    assert uncouple(JzKetCoupled(1,-1,S(1)/2,S(1)/2)) == \
        1.0*TensorProduct(JzKet(S(1)/2,-S(1)/2),JzKet(S(1)/2,-S(1)/2))
    assert uncouple(JzKetCoupled(S(1)/2,S(1)/2,1,S(1)/2)) == \
        -sqrt(3)*TensorProduct(JzKet(1,0),JzKet(S(1)/2,S(1)/2))/3 + sqrt(6)*TensorProduct(JzKet(1,1),JzKet(S(1)/2,-S(1)/2))/3
    assert uncouple(JzKetCoupled(S(1)/2,-S(1)/2,1,S(1)/2)) == \
        sqrt(3)*TensorProduct(JzKet(1,0),JzKet(S(1)/2,-S(1)/2))/3 - sqrt(6)*TensorProduct(JzKet(1,-1),JzKet(S(1)/2,S(1)/2))/3
    assert uncouple(JzKetCoupled(S(3)/2,S(3)/2,1,S(1)/2)) == \
        TensorProduct(JzKet(1,1),JzKet(S(1)/2,S(1)/2))
    assert uncouple(JzKetCoupled(S(3)/2,S(1)/2,1,S(1)/2)) == \
        sqrt(3)*TensorProduct(JzKet(1,1),JzKet(S(1)/2,-S(1)/2))/3 + sqrt(6)*TensorProduct(JzKet(1,0),JzKet(S(1)/2,S(1)/2))/3
    assert uncouple(JzKetCoupled(S(3)/2,-S(1)/2,1,S(1)/2)) == \
        sqrt(6)*TensorProduct(JzKet(1,0),JzKet(S(1)/2,-S(1)/2))/3 + sqrt(3)*TensorProduct(JzKet(1,-1),JzKet(S(1)/2,S(1)/2))/3
    assert uncouple(JzKetCoupled(S(3)/2,-S(3)/2,1,S(1)/2)) == \
        1.0*TensorProduct(JzKet(1,-1),JzKet(S(1)/2,-S(1)/2))
    assert uncouple(JzKetCoupled(0,0,1,1)) == \
        sqrt(3)*TensorProduct(JzKet(1,1),JzKet(1,-1))/3 - sqrt(3)*TensorProduct(JzKet(1,0),JzKet(1,0))/3 + sqrt(3)*TensorProduct(JzKet(1,-1),JzKet(1,1))/3
    assert uncouple(JzKetCoupled(1,1,1,1)) == \
        sqrt(2)*TensorProduct(JzKet(1,1),JzKet(1,0))/2 - sqrt(2)*TensorProduct(JzKet(1,0),JzKet(1,1))/2
    assert uncouple(JzKetCoupled(1,0,1,1)) == \
        sqrt(2)*TensorProduct(JzKet(1,1),JzKet(1,-1))/2 - sqrt(2)*TensorProduct(JzKet(1,-1),JzKet(1,1))/2
    assert uncouple(JzKetCoupled(1,-1,1,1)) == \
        0.5*sqrt(2)*TensorProduct(JzKet(1,0),JzKet(1,-1)) - 0.5*sqrt(2)*TensorProduct(JzKet(1,-1),JzKet(1,0))
    assert uncouple(JzKetCoupled(2,2,1,1)) == \
        TensorProduct(JzKet(1,1),JzKet(1,1))
    assert uncouple(JzKetCoupled(2,1,1,1)) == \
        sqrt(2)*TensorProduct(JzKet(1,1),JzKet(1,0))/2 + sqrt(2)*TensorProduct(JzKet(1,0),JzKet(1,1))/2
    assert uncouple(JzKetCoupled(2,0,1,1)) == \
        sqrt(6)*TensorProduct(JzKet(1,1),JzKet(1,-1))/6 + sqrt(6)*TensorProduct(JzKet(1,0),JzKet(1,0))/3 + sqrt(6)*TensorProduct(JzKet(1,-1),JzKet(1,1))/6
    assert uncouple(JzKetCoupled(2,-1,1,1)) == \
        0.5*sqrt(2)*TensorProduct(JzKet(1,0),JzKet(1,-1)) + 0.5*sqrt(2)*TensorProduct(JzKet(1,-1),JzKet(1,0))
    assert uncouple(JzKetCoupled(2,-2,1,1)) == \
        1.0*TensorProduct(JzKet(1,-1),JzKet(1,-1))
    # Symbolic
    assert uncouple(JzKetCoupled(j,m,j1,j2)) == \
        Sum(CG(j1,m1,j2,m2,j,m) * TensorProduct(JzKet(j1,m1), JzKet(j2,m2)), (m1, -j1, j1), (m2, -j2, j2))

def test_couple():
    # Numerical
    assert couple(TensorProduct(JzKet(S(1)/2,S(1)/2), JzKet(S(1)/2,S(1)/2))) == \
        JzKetCoupled(1,1,S(1)/2,S(1)/2)
    assert couple(TensorProduct(JzKet(S(1)/2,S(1)/2), JzKet(S(1)/2,-S(1)/2))) == \
        sqrt(2)*JzKetCoupled(0,0,S(1)/2,S(1)/2)/2 + sqrt(2)*JzKetCoupled(1,0,S(1)/2,S(1)/2)/2
    assert couple(TensorProduct(JzKet(S(1)/2,-S(1)/2), JzKet(S(1)/2,S(1)/2))) == \
        -sqrt(2)*JzKetCoupled(0,0,S(1)/2,S(1)/2)/2 + sqrt(2)*JzKetCoupled(1,0,S(1)/2,S(1)/2)/2
    assert couple(TensorProduct(JzKet(S(1)/2,-S(1)/2), JzKet(S(1)/2,-S(1)/2))) == \
        1.0*JzKetCoupled(1,-1,S(1)/2,S(1)/2)
    assert couple(TensorProduct(JzKet(1,1), JzKet(S(1)/2,S(1)/2))) == \
        JzKetCoupled(S(3)/2,S(3)/2,1,S(1)/2)
    assert couple(TensorProduct(JzKet(1,1), JzKet(S(1)/2,-S(1)/2))) == \
        sqrt(6)*JzKetCoupled(S(1)/2,S(1)/2,1,S(1)/2)/3 + sqrt(3)*JzKetCoupled(S(3)/2,S(1)/2,1,S(1)/2)/3
    assert couple(TensorProduct(JzKet(1,0), JzKet(S(1)/2,S(1)/2))) == \
        -sqrt(3)*JzKetCoupled(S(1)/2,S(1)/2,1,S(1)/2)/3 + sqrt(6)*JzKetCoupled(S(3)/2,S(1)/2,1,S(1)/2)/3
    assert couple(TensorProduct(JzKet(1,0), JzKet(S(1)/2,-S(1)/2))) == \
        sqrt(3)*JzKetCoupled(S(1)/2,-S(1)/2,1,S(1)/2)/3 + sqrt(6)*JzKetCoupled(S(3)/2,-S(1)/2,1,S(1)/2)/3
    assert couple(TensorProduct(JzKet(1,-1), JzKet(S(1)/2,S(1)/2))) == \
        -sqrt(6)*JzKetCoupled(S(1)/2,-S(1)/2,1,S(1)/2)/3 + sqrt(3)*JzKetCoupled(S(3)/2,-S(1)/2,1,S(1)/2)/3
    assert couple(TensorProduct(JzKet(1,-1), JzKet(S(1)/2,-S(1)/2))) == \
        1.0*JzKetCoupled(S(3)/2,-S(3)/2,1,S(1)/2)
    assert couple(TensorProduct(JzKet(1,1), JzKet(1,1))) == \
        JzKetCoupled(2,2,1,1)
    assert couple(TensorProduct(JzKet(1,1), JzKet(1,0))) == \
        sqrt(2)*JzKetCoupled(1,1,1,1)/2 + sqrt(2)*JzKetCoupled(2,1,1,1)/2
    assert couple(TensorProduct(JzKet(1,1), JzKet(1,-1))) == \
        sqrt(3)*JzKetCoupled(0,0,1,1)/3 + sqrt(2)*JzKetCoupled(1,0,1,1)/2 + sqrt(6)*JzKetCoupled(2,0,1,1)/6
    assert couple(TensorProduct(JzKet(1,0), JzKet(1,1))) == \
        -sqrt(2)*JzKetCoupled(1,1,1,1)/2 + sqrt(2)*JzKetCoupled(2,1,1,1)/2
    assert couple(TensorProduct(JzKet(1,0), JzKet(1,0))) == \
        -sqrt(3)*JzKetCoupled(0,0,1,1)/3+sqrt(6)*JzKetCoupled(2,0,1,1)/3
    assert couple(TensorProduct(JzKet(1,0), JzKet(1,-1))) == \
        0.5*sqrt(2)*JzKetCoupled(1,-1,1,1) + 0.5*sqrt(2)*JzKetCoupled(2,-1,1,1)
    assert couple(TensorProduct(JzKet(1,-1), JzKet(1,1))) == \
        sqrt(3)*JzKetCoupled(0,0,1,1)/3 - sqrt(2)*JzKetCoupled(1,0,1,1)/2 + sqrt(6)*JzKetCoupled(2,0,1,1)/6
    assert couple(TensorProduct(JzKet(1,-1), JzKet(1,0))) == \
        -0.5*sqrt(2)*JzKetCoupled(1,-1,1,1) + 0.5*sqrt(2)*JzKetCoupled(2,-1,1,1)
    assert couple(TensorProduct(JzKet(1,-1), JzKet(1,-1))) == \
        1.0*JzKetCoupled(2,-2,1,1)
    # Numerical
    assert couple(TensorProduct(JzKet(j1,m1), JzKet(j2,m2))) == \
        Sum(CG(j1,m1,j2,m2,j,m1+m2) * JzKetCoupled(j,m1+m2), (j,0,j1+j2))

def test_innerproduct():
    assert InnerProduct(JzBra(1,1), JzKet(1,1)).doit() == 1
    assert InnerProduct(JzBra(S(1)/2,S(1)/2), JzKet(S(1)/2,-S(1)/2)).doit() == 0
    assert InnerProduct(JzBra(j,m), JzKet(j,m)).doit() == 1
    assert InnerProduct(JzBra(1,0), JyKet(1,1)).doit() == I/sqrt(2)
    assert InnerProduct(JxBra(S(1)/2,S(1)/2), JzKet(S(1)/2,S(1)/2)).doit() == -sqrt(2)/2
    assert InnerProduct(JyBra(1,1), JzKet(1,1)).doit() == S(1)/2
    assert InnerProduct(JxBra(1,-1), JyKet(1,1)).doit() == 0

def test_rotation_small_d():
    # Symbolic tests
    # j = 1/2
    assert Rotation.d(S(1)/2,S(1)/2,S(1)/2,beta).doit() == cos(beta/2)
    assert Rotation.d(S(1)/2,S(1)/2,-S(1)/2,beta).doit() == -sin(beta/2)
    assert Rotation.d(S(1)/2,-S(1)/2,S(1)/2,beta).doit() == sin(beta/2)
    assert Rotation.d(S(1)/2,-S(1)/2,-S(1)/2,beta).doit() == cos(beta/2)
    # j = 1
    assert Rotation.d(1,1,1,beta).doit() == (1+cos(beta))/2
    assert Rotation.d(1,1,0,beta).doit() == -sin(beta)/sqrt(2)
    assert Rotation.d(1,1,-1,beta).doit() == (1-cos(beta))/2
    assert Rotation.d(1,0,1,beta).doit() == sin(beta)/sqrt(2)
    assert Rotation.d(1,0,0,beta).doit() == cos(beta)
    assert Rotation.d(1,0,-1,beta).doit() == -sin(beta)/sqrt(2)
    assert Rotation.d(1,-1,1,beta).doit() == (1-cos(beta))/2
    assert Rotation.d(1,-1,0,beta).doit() == sin(beta)/sqrt(2)
    assert Rotation.d(1,-1,-1,beta).doit() == (1+cos(beta))/2
    # j = 3/2
    assert Rotation.d(S(3)/2,S(3)/2,S(3)/2,beta).doit() == (3*cos(beta/2)+cos(3*beta/2))/4
    assert Rotation.d(S(3)/2,S(3)/2,S(1)/2,beta).doit() == sqrt(3)*(-sin(beta/2)-sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,S(3)/2,-S(1)/2,beta).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4
    assert Rotation.d(S(3)/2,S(3)/2,-S(3)/2,beta).doit() == (-3*sin(beta/2)+sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,S(1)/2,S(3)/2,beta).doit() == sqrt(3)*(sin(beta/2)+sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,S(1)/2,S(1)/2,beta).doit() == (cos(beta/2)+3*cos(3*beta/2))/4
    assert Rotation.d(S(3)/2,S(1)/2,-S(1)/2,beta).doit() == (sin(beta/2)-3*sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,S(1)/2,-S(3)/2,beta).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(1)/2,S(3)/2,beta).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(1)/2,S(1)/2,beta).doit() == (-sin(beta/2)+3*sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(1)/2,-S(1)/2,beta).doit() == (cos(beta/2)+3*cos(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(1)/2,-S(3)/2,beta).doit() == sqrt(3)*(-sin(beta/2)-sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(3)/2,S(3)/2,beta).doit() == (3*sin(beta/2)-sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(3)/2,S(1)/2,beta).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(3)/2,-S(1)/2,beta).doit() == sqrt(3)*(sin(beta/2)+sin(3*beta/2))/4
    assert Rotation.d(S(3)/2,-S(3)/2,-S(3)/2,beta).doit() == (3*cos(beta/2)+cos(3*beta/2))/4
    # j = 2
    assert Rotation.d(2,2,2,beta).doit() == (3+4*cos(beta)+cos(2*beta))/8
    assert Rotation.d(2,2,1,beta).doit() == (-2*sin(beta)-sin(2*beta))/4
    assert Rotation.d(2,2,0,beta).doit() == sqrt(6)*(1-cos(2*beta))/8
    assert Rotation.d(2,2,-1,beta).doit() == (-2*sin(beta)+sin(2*beta))/4
    assert Rotation.d(2,2,-2,beta).doit() == (3-4*cos(beta)+cos(2*beta))/8
    assert Rotation.d(2,1,2,beta).doit() == (2*sin(beta)+sin(2*beta))/4
    assert Rotation.d(2,1,1,beta).doit() == (cos(beta)+cos(2*beta))/2
    assert Rotation.d(2,1,0,beta).doit() == -sqrt(6)*sin(2*beta)/4
    assert Rotation.d(2,1,-1,beta).doit() == (cos(beta)-cos(2*beta))/2
    assert Rotation.d(2,1,-2,beta).doit() == (-2*sin(beta)+sin(2*beta))/4
    assert Rotation.d(2,0,2,beta).doit() == sqrt(6)*(1-cos(2*beta))/8
    assert Rotation.d(2,0,1,beta).doit() == sqrt(6)*sin(2*beta)/4
    assert Rotation.d(2,0,0,beta).doit() == (1+3*cos(2*beta))/4
    assert Rotation.d(2,0,-1,beta).doit() == -sqrt(6)*sin(2*beta)/4
    assert Rotation.d(2,0,-2,beta).doit() == sqrt(6)*(1-cos(2*beta))/8
    assert Rotation.d(2,-1,2,beta).doit() == (2*sin(beta)-sin(2*beta))/4
    assert Rotation.d(2,-1,1,beta).doit() == (cos(beta)-cos(2*beta))/2
    assert Rotation.d(2,-1,0,beta).doit() == sqrt(6)*sin(2*beta)/4
    assert Rotation.d(2,-1,-1,beta).doit() == (cos(beta)+cos(2*beta))/2
    assert Rotation.d(2,-1,-2,beta).doit() == (-2*sin(beta)-sin(2*beta))/4
    assert Rotation.d(2,-2,2,beta).doit() == (3-4*cos(beta)+cos(2*beta))/8
    assert Rotation.d(2,-2,1,beta).doit() == (2*sin(beta)-sin(2*beta))/4
    assert Rotation.d(2,-2,0,beta).doit() == sqrt(6)*(1-cos(2*beta))/8
    assert Rotation.d(2,-2,-1,beta).doit() == (2*sin(beta)+sin(2*beta))/4
    assert Rotation.d(2,-2,-2,beta).doit() == (3+4*cos(beta)+cos(2*beta))/8
    # Numerical tests
    # j = 1/2
    assert Rotation.d(S(1)/2,S(1)/2,S(1)/2,pi/2).doit() == sqrt(2)/2
    assert Rotation.d(S(1)/2,S(1)/2,-S(1)/2,pi/2).doit() == -sqrt(2)/2
    assert Rotation.d(S(1)/2,-S(1)/2,S(1)/2,pi/2).doit() == sqrt(2)/2
    assert Rotation.d(S(1)/2,-S(1)/2,-S(1)/2,pi/2).doit() == sqrt(2)/2
    # j = 1
    assert Rotation.d(1,1,1,pi/2).doit() == 1/2
    assert Rotation.d(1,1,0,pi/2).doit() == -sqrt(2)/2
    assert Rotation.d(1,1,-1,pi/2).doit() == 1/2
    assert Rotation.d(1,0,1,pi/2).doit() == sqrt(2)/2
    assert Rotation.d(1,0,0,pi/2).doit() == 0
    assert Rotation.d(1,0,-1,pi/2).doit() == -sqrt(2)/2
    assert Rotation.d(1,-1,1,pi/2).doit() == 1/2
    assert Rotation.d(1,-1,0,pi/2).doit() == sqrt(2)/2
    assert Rotation.d(1,-1,-1,pi/2).doit() == 1/2
    # j = 3/2
    assert Rotation.d(S(3)/2,S(3)/2,S(3)/2,pi/2).doit() == sqrt(2)/4
    assert Rotation.d(S(3)/2,S(3)/2,S(1)/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.d(S(3)/2,S(3)/2,-S(1)/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(S(3)/2,S(3)/2,-S(3)/2,pi/2).doit() == -sqrt(2)/4
    assert Rotation.d(S(3)/2,S(1)/2,S(3)/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(S(3)/2,S(1)/2,S(1)/2,pi/2).doit() == -sqrt(2)/4
    assert Rotation.d(S(3)/2,S(1)/2,-S(1)/2,pi/2).doit() == -sqrt(2)/4
    assert Rotation.d(S(3)/2,S(1)/2,-S(3)/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(S(3)/2,-S(1)/2,S(3)/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(S(3)/2,-S(1)/2,S(1)/2,pi/2).doit() == sqrt(2)/4
    assert Rotation.d(S(3)/2,-S(1)/2,-S(1)/2,pi/2).doit() == -sqrt(2)/4
    assert Rotation.d(S(3)/2,-S(1)/2,-S(3)/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.d(S(3)/2,-S(3)/2,S(3)/2,pi/2).doit() == sqrt(2)/4
    assert Rotation.d(S(3)/2,-S(3)/2,S(1)/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(S(3)/2,-S(3)/2,-S(1)/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(S(3)/2,-S(3)/2,-S(3)/2,pi/2).doit() == sqrt(2)/4
    # j = 2
    assert Rotation.d(2,2,2,pi/2).doit() == 1/4
    assert Rotation.d(2,2,1,pi/2).doit() == -1/2
    assert Rotation.d(2,2,0,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(2,2,-1,pi/2).doit() == -1/2
    assert Rotation.d(2,2,-2,pi/2).doit() == 1/4
    assert Rotation.d(2,1,2,pi/2).doit() == 1/2
    assert Rotation.d(2,1,1,pi/2).doit() == -1/2
    assert Rotation.d(2,1,0,pi/2).doit() == 0
    assert Rotation.d(2,1,-1,pi/2).doit() == 1/2
    assert Rotation.d(2,1,-2,pi/2).doit() == -1/2
    assert Rotation.d(2,0,2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(2,0,1,pi/2).doit() == 0
    assert Rotation.d(2,0,0,pi/2).doit() == -1/2
    assert Rotation.d(2,0,-1,pi/2).doit() == 0
    assert Rotation.d(2,0,-2,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(2,-1,2,pi/2).doit() == 1/2
    assert Rotation.d(2,-1,1,pi/2).doit() == 1/2
    assert Rotation.d(2,-1,0,pi/2).doit() == 0
    assert Rotation.d(2,-1,-1,pi/2).doit() == -1/2
    assert Rotation.d(2,-1,-2,pi/2).doit() == -1/2
    assert Rotation.d(2,-2,2,pi/2).doit() == 1/4
    assert Rotation.d(2,-2,1,pi/2).doit() == 1/2
    assert Rotation.d(2,-2,0,pi/2).doit() == sqrt(6)/4
    assert Rotation.d(2,-2,-1,pi/2).doit() == 1/2
    assert Rotation.d(2,-2,-2,pi/2).doit() == 1/4

def test_rotation_d():
    # Symbolic tests
    # j = 1/2
    assert Rotation.D(S(1)/2,S(1)/2,S(1)/2,alpha,beta,gamma).doit() == cos(beta/2)*exp(-I*alpha/2)*exp(-I*gamma/2)
    assert Rotation.D(S(1)/2,S(1)/2,-S(1)/2,alpha,beta,gamma).doit() == -sin(beta/2)*exp(-I*alpha/2)*exp(I*gamma/2)
    assert Rotation.D(S(1)/2,-S(1)/2,S(1)/2,alpha,beta,gamma).doit() == sin(beta/2)*exp(I*alpha/2)*exp(-I*gamma/2)
    assert Rotation.D(S(1)/2,-S(1)/2,-S(1)/2,alpha,beta,gamma).doit() == cos(beta/2)*exp(I*alpha/2)*exp(I*gamma/2)
    # j = 1
    assert Rotation.D(1,1,1,alpha,beta,gamma).doit() == (1+cos(beta))/2*exp(-I*alpha)*exp(-I*gamma)
    assert Rotation.D(1,1,0,alpha,beta,gamma).doit() == -sin(beta)/sqrt(2)*exp(-I*alpha)
    assert Rotation.D(1,1,-1,alpha,beta,gamma).doit() == (1-cos(beta))/2*exp(-I*alpha)*exp(I*gamma)
    assert Rotation.D(1,0,1,alpha,beta,gamma).doit() == sin(beta)/sqrt(2)*exp(-I*gamma)
    assert Rotation.D(1,0,0,alpha,beta,gamma).doit() == cos(beta)
    assert Rotation.D(1,0,-1,alpha,beta,gamma).doit() == -sin(beta)/sqrt(2)*exp(I*gamma)
    assert Rotation.D(1,-1,1,alpha,beta,gamma).doit() == (1-cos(beta))/2*exp(I*alpha)*exp(-I*gamma)
    assert Rotation.D(1,-1,0,alpha,beta,gamma).doit() == sin(beta)/sqrt(2)*exp(I*alpha)
    assert Rotation.D(1,-1,-1,alpha,beta,gamma).doit() == (1+cos(beta))/2*exp(I*alpha)*exp(I*gamma)
    # j = 3/2
    assert Rotation.D(S(3)/2,S(3)/2,S(3)/2,alpha,beta,gamma).doit() == (3*cos(beta/2)+cos(3*beta/2))/4*exp(-3*I*alpha/2)*exp(-3*I*gamma/2)
    assert Rotation.D(S(3)/2,S(3)/2,S(1)/2,alpha,beta,gamma).doit() == sqrt(3)*(-sin(beta/2)-sin(3*beta/2))/4*exp(-3*I*alpha/2)*exp(-I*gamma/2)
    assert Rotation.D(S(3)/2,S(3)/2,-S(1)/2,alpha,beta,gamma).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4*exp(-3*I*alpha/2)*exp(I*gamma/2)
    assert Rotation.D(S(3)/2,S(3)/2,-S(3)/2,alpha,beta,gamma).doit() == (-3*sin(beta/2)+sin(3*beta/2))/4*exp(-3*I*alpha/2)*exp(3*I*gamma/2)
    assert Rotation.D(S(3)/2,S(1)/2,S(3)/2,alpha,beta,gamma).doit() == sqrt(3)*(sin(beta/2)+sin(3*beta/2))/4*exp(-I*alpha/2)*exp(-3*I*gamma/2)
    assert Rotation.D(S(3)/2,S(1)/2,S(1)/2,alpha,beta,gamma).doit() == (cos(beta/2)+3*cos(3*beta/2))/4*exp(-I*alpha/2)*exp(-I*gamma/2)
    assert Rotation.D(S(3)/2,S(1)/2,-S(1)/2,alpha,beta,gamma).doit() == (sin(beta/2)-3*sin(3*beta/2))/4*exp(-I*alpha/2)*exp(I*gamma/2)
    assert Rotation.D(S(3)/2,S(1)/2,-S(3)/2,alpha,beta,gamma).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4*exp(-I*alpha/2)*exp(3*I*gamma/2)
    assert Rotation.D(S(3)/2,-S(1)/2,S(3)/2,alpha,beta,gamma).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4*exp(I*alpha/2)*exp(-3*I*gamma/2)
    assert Rotation.D(S(3)/2,-S(1)/2,S(1)/2,alpha,beta,gamma).doit() == (-sin(beta/2)+3*sin(3*beta/2))/4*exp(I*alpha/2)*exp(-I*gamma/2)
    assert Rotation.D(S(3)/2,-S(1)/2,-S(1)/2,alpha,beta,gamma).doit() == (cos(beta/2)+3*cos(3*beta/2))/4*exp(I*alpha/2)*exp(I*gamma/2)
    assert Rotation.D(S(3)/2,-S(1)/2,-S(3)/2,alpha,beta,gamma).doit() == sqrt(3)*(-sin(beta/2)-sin(3*beta/2))/4*exp(I*alpha/2)*exp(3*I*gamma/2)
    assert Rotation.D(S(3)/2,-S(3)/2,S(3)/2,alpha,beta,gamma).doit() == (3*sin(beta/2)-sin(3*beta/2))/4*exp(3*I*alpha/2)*exp(-3*I*gamma/2)
    assert Rotation.D(S(3)/2,-S(3)/2,S(1)/2,alpha,beta,gamma).doit() == sqrt(3)*(cos(beta/2)-cos(3*beta/2))/4*exp(3*I*alpha/2)*exp(-I*gamma/2)
    assert Rotation.D(S(3)/2,-S(3)/2,-S(1)/2,alpha,beta,gamma).doit() == sqrt(3)*(sin(beta/2)+sin(3*beta/2))/4*exp(3*I*alpha/2)*exp(I*gamma/2)
    assert Rotation.D(S(3)/2,-S(3)/2,-S(3)/2,alpha,beta,gamma).doit() == (3*cos(beta/2)+cos(3*beta/2))/4*exp(3*I*alpha/2)*exp(3*I*gamma/2)
    # j = 2
    assert Rotation.D(2,2,2,alpha,beta,gamma).doit() == (3+4*cos(beta)+cos(2*beta))/8*exp(-2*I*alpha)*exp(-2*I*gamma)
    assert Rotation.D(2,2,1,alpha,beta,gamma).doit() == (-2*sin(beta)-sin(2*beta))/4*exp(-2*I*alpha)*exp(-I*gamma)
    assert Rotation.D(2,2,0,alpha,beta,gamma).doit() == sqrt(6)*(1-cos(2*beta))/8*exp(-2*I*alpha)
    assert Rotation.D(2,2,-1,alpha,beta,gamma).doit() == (-2*sin(beta)+sin(2*beta))/4*exp(-2*I*alpha)*exp(I*gamma)
    assert Rotation.D(2,2,-2,alpha,beta,gamma).doit() == (3-4*cos(beta)+cos(2*beta))/8*exp(-2*I*alpha)*exp(2*I*gamma)
    assert Rotation.D(2,1,2,alpha,beta,gamma).doit() == (2*sin(beta)+sin(2*beta))/4*exp(-I*alpha)*exp(-2*I*gamma)
    assert Rotation.D(2,1,1,alpha,beta,gamma).doit() == (cos(beta)+cos(2*beta))/2*exp(-I*alpha)*exp(-I*gamma)
    assert Rotation.D(2,1,0,alpha,beta,gamma).doit() == -sqrt(6)*sin(2*beta)/4*exp(-I*alpha)
    assert Rotation.D(2,1,-1,alpha,beta,gamma).doit() == (cos(beta)-cos(2*beta))/2*exp(-I*alpha)*exp(I*gamma)
    assert Rotation.D(2,1,-2,alpha,beta,gamma).doit() == (-2*sin(beta)+sin(2*beta))/4*exp(-I*alpha)*exp(2*I*gamma)
    assert Rotation.D(2,0,2,alpha,beta,gamma).doit() == sqrt(6)*(1-cos(2*beta))/8*exp(-2*I*gamma)
    assert Rotation.D(2,0,1,alpha,beta,gamma).doit() == sqrt(6)*sin(2*beta)/4*exp(-I*gamma)
    assert Rotation.D(2,0,0,alpha,beta,gamma).doit() == (1+3*cos(2*beta))/4
    assert Rotation.D(2,0,-1,alpha,beta,gamma).doit() == -sqrt(6)*sin(2*beta)/4*exp(I*gamma)
    assert Rotation.D(2,0,-2,alpha,beta,gamma).doit() == sqrt(6)*(1-cos(2*beta))/8*exp(2*I*gamma)
    assert Rotation.D(2,-1,2,alpha,beta,gamma).doit() == (2*sin(beta)-sin(2*beta))/4*exp(I*alpha)*exp(-2*I*gamma)
    assert Rotation.D(2,-1,1,alpha,beta,gamma).doit() == (cos(beta)-cos(2*beta))/2*exp(I*alpha)*exp(-I*gamma)
    assert Rotation.D(2,-1,0,alpha,beta,gamma).doit() == sqrt(6)*sin(2*beta)/4*exp(I*alpha)
    assert Rotation.D(2,-1,-1,alpha,beta,gamma).doit() == (cos(beta)+cos(2*beta))/2*exp(I*alpha)*exp(I*gamma)
    assert Rotation.D(2,-1,-2,alpha,beta,gamma).doit() == (-2*sin(beta)-sin(2*beta))/4*exp(I*alpha)*exp(2*I*gamma)
    assert Rotation.D(2,-2,2,alpha,beta,gamma).doit() == (3-4*cos(beta)+cos(2*beta))/8*exp(2*I*alpha)*exp(-2*I*gamma)
    assert Rotation.D(2,-2,1,alpha,beta,gamma).doit() == (2*sin(beta)-sin(2*beta))/4*exp(2*I*alpha)*exp(-I*gamma)
    assert Rotation.D(2,-2,0,alpha,beta,gamma).doit() == sqrt(6)*(1-cos(2*beta))/8*exp(2*I*alpha)
    assert Rotation.D(2,-2,-1,alpha,beta,gamma).doit() == (2*sin(beta)+sin(2*beta))/4*exp(2*I*alpha)*exp(I*gamma)
    assert Rotation.D(2,-2,-2,alpha,beta,gamma).doit() == (3+4*cos(beta)+cos(2*beta))/8*exp(2*I*alpha)*exp(2*I*gamma)
    # Numerical tests
    # j = 1/2
    assert Rotation.D(S(1)/2,S(1)/2,S(1)/2,pi/2,pi/2,pi/2).doit() == -I*sqrt(2)/2
    assert Rotation.D(S(1)/2,S(1)/2,-S(1)/2,pi/2,pi/2,pi/2).doit() == -sqrt(2)/2
    assert Rotation.D(S(1)/2,-S(1)/2,S(1)/2,pi/2,pi/2,pi/2).doit() == sqrt(2)/2
    assert Rotation.D(S(1)/2,-S(1)/2,-S(1)/2,pi/2,pi/2,pi/2).doit() == I*sqrt(2)/2
    # j = 1
    assert Rotation.D(1,1,1,pi/2,pi/2,pi/2).doit() == -1/2
    assert Rotation.D(1,1,0,pi/2,pi/2,pi/2).doit() == I*sqrt(2)/2
    assert Rotation.D(1,1,-1,pi/2,pi/2,pi/2).doit() == 1/2
    assert Rotation.D(1,0,1,pi/2,pi/2,pi/2).doit() == -I*sqrt(2)/2
    assert Rotation.D(1,0,0,pi/2,pi/2,pi/2).doit() == 0
    assert Rotation.D(1,0,-1,pi/2,pi/2,pi/2).doit() == -I*sqrt(2)/2
    assert Rotation.D(1,-1,1,pi/2,pi/2,pi/2).doit() == 1/2
    assert Rotation.D(1,-1,0,pi/2,pi/2,pi/2).doit() == I*sqrt(2)/2
    assert Rotation.D(1,-1,-1,pi/2,pi/2,pi/2).doit() == -1/2
    # j = 3/2
    assert Rotation.D(S(3)/2,S(3)/2,S(3)/2,pi/2,pi/2,pi/2).doit() == I*sqrt(2)/4
    assert Rotation.D(S(3)/2,S(3)/2,S(1)/2,pi/2,pi/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.D(S(3)/2,S(3)/2,-S(1)/2,pi/2,pi/2,pi/2).doit() == -I*sqrt(6)/4
    assert Rotation.D(S(3)/2,S(3)/2,-S(3)/2,pi/2,pi/2,pi/2).doit() == -sqrt(2)/4
    assert Rotation.D(S(3)/2,S(1)/2,S(3)/2,pi/2,pi/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.D(S(3)/2,S(1)/2,S(1)/2,pi/2,pi/2,pi/2).doit() == I*sqrt(2)/4
    assert Rotation.D(S(3)/2,S(1)/2,-S(1)/2,pi/2,pi/2,pi/2).doit() == -sqrt(2)/4
    assert Rotation.D(S(3)/2,S(1)/2,-S(3)/2,pi/2,pi/2,pi/2).doit() == I*sqrt(6)/4
    assert Rotation.D(S(3)/2,-S(1)/2,S(3)/2,pi/2,pi/2,pi/2).doit() == -I*sqrt(6)/4
    assert Rotation.D(S(3)/2,-S(1)/2,S(1)/2,pi/2,pi/2,pi/2).doit() == sqrt(2)/4
    assert Rotation.D(S(3)/2,-S(1)/2,-S(1)/2,pi/2,pi/2,pi/2).doit() == -I*sqrt(2)/4
    assert Rotation.D(S(3)/2,-S(1)/2,-S(3)/2,pi/2,pi/2,pi/2).doit() == sqrt(6)/4
    assert Rotation.D(S(3)/2,-S(3)/2,S(3)/2,pi/2,pi/2,pi/2).doit() == sqrt(2)/4
    assert Rotation.D(S(3)/2,-S(3)/2,S(1)/2,pi/2,pi/2,pi/2).doit() == I*sqrt(6)/4
    assert Rotation.D(S(3)/2,-S(3)/2,-S(1)/2,pi/2,pi/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.D(S(3)/2,-S(3)/2,-S(3)/2,pi/2,pi/2,pi/2).doit() == -I*sqrt(2)/4
    # j = 2
    assert Rotation.D(2,2,2,pi/2,pi/2,pi/2).doit() == 1/4
    assert Rotation.D(2,2,1,pi/2,pi/2,pi/2).doit() == -I/2
    assert Rotation.D(2,2,0,pi/2,pi/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.D(2,2,-1,pi/2,pi/2,pi/2).doit() == I/2
    assert Rotation.D(2,2,-2,pi/2,pi/2,pi/2).doit() == 1/4
    assert Rotation.D(2,1,2,pi/2,pi/2,pi/2).doit() == I/2
    assert Rotation.D(2,1,1,pi/2,pi/2,pi/2).doit() == 1/2
    assert Rotation.D(2,1,0,pi/2,pi/2,pi/2).doit() == 0
    assert Rotation.D(2,1,-1,pi/2,pi/2,pi/2).doit() == 1/2
    assert Rotation.D(2,1,-2,pi/2,pi/2,pi/2).doit() == -I/2
    assert Rotation.D(2,0,2,pi/2,pi/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.D(2,0,1,pi/2,pi/2,pi/2).doit() == 0
    assert Rotation.D(2,0,0,pi/2,pi/2,pi/2).doit() == -1/2
    assert Rotation.D(2,0,-1,pi/2,pi/2,pi/2).doit() == 0
    assert Rotation.D(2,0,-2,pi/2,pi/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.D(2,-1,2,pi/2,pi/2,pi/2).doit() == -I/2
    assert Rotation.D(2,-1,1,pi/2,pi/2,pi/2).doit() == 1/2
    assert Rotation.D(2,-1,0,pi/2,pi/2,pi/2).doit() == 0
    assert Rotation.D(2,-1,-1,pi/2,pi/2,pi/2).doit() == 1/2
    assert Rotation.D(2,-1,-2,pi/2,pi/2,pi/2).doit() == I/2
    assert Rotation.D(2,-2,2,pi/2,pi/2,pi/2).doit() == 1/4
    assert Rotation.D(2,-2,1,pi/2,pi/2,pi/2).doit() == I/2
    assert Rotation.D(2,-2,0,pi/2,pi/2,pi/2).doit() == -sqrt(6)/4
    assert Rotation.D(2,-2,-1,pi/2,pi/2,pi/2).doit() == -I/2
    assert Rotation.D(2,-2,-2,pi/2,pi/2,pi/2).doit() == 1/4

def test_wignerd():
    assert Rotation.D(j, m, mp, alpha, beta, gamma) == WignerD(j, m, mp, alpha, beta, gamma)
    assert Rotation.d(j, m, mp, beta) == WignerD(j, m, mp, 0, beta, 0)

def test_jplus():
    assert Commutator(Jplus, Jminus).doit() == 2*hbar*Jz
    assert qapply(Jplus*JzKet(1,1)) == 0
    assert Jplus.matrix_element(1,1,1,1) == 0
    assert Jplus.rewrite('xyz') == Jx + I*Jy

def test_jminus():
    assert qapply(Jminus*JzKet(1,-1)) == 0
    assert Jminus.matrix_element(1,0,1,1) == sqrt(2)*hbar
    assert Jminus.rewrite('xyz') == Jx - I*Jy

def test_j2():
    assert Commutator(J2, Jz).doit() == 0
    assert qapply(J2*JzKet(1,1)) == 2*hbar**2*JzKet(1,1)
    assert qapply(J2*JzKet(j,m)) == j**2*hbar**2*JzKet(j,m)+j*hbar**2*JzKet(j,m)
    assert J2.matrix_element(1,1,1,1) == 2*hbar**2

def test_jx():
    assert Commutator(Jx, Jz).doit() == -I*hbar*Jy
    assert qapply(Jx*JzKet(1,1)) == sqrt(2)*hbar*JzKet(1,0)/2
    assert Jx.rewrite('plusminus') == (Jminus + Jplus)/2
    assert represent(Jx, basis=Jz, j=1) == (represent(Jplus, basis=Jz, j=1)+represent(Jminus, basis=Jz, j=1))/2
    # Normal Operators
    # Numerical
    assert qapply(Jx*JxKet(1,1)) == hbar*JxKet(1,1)
    assert qapply(Jx*JyKet(1,1)) == hbar*JyKet(1,1)
    assert qapply(Jx*JzKet(1,1)) == sqrt(2)*hbar*JzKet(1,0)/2
    assert qapply(Jx*TensorProduct(JxKet(1,1), JxKet(1,1))) == 2*hbar*TensorProduct(JxKet(1,1), JxKet(1,1))
    assert qapply(Jx*TensorProduct(JyKet(1,1), JyKet(1,1))) == \
        hbar*TensorProduct(JyKet(1,1),JyKet(1,1))+hbar*TensorProduct(JyKet(1,1),JyKet(1,1))
    assert qapply(Jx*TensorProduct(JzKet(1,1), JzKet(1,1))) == \
        sqrt(2)*hbar*TensorProduct(JzKet(1,1),JzKet(1,0))/2+sqrt(2)*hbar*TensorProduct(JzKet(1,0),JzKet(1,1))/2
    assert qapply(Jx*TensorProduct(JxKet(1,1), JxKet(1,-1))) == 0
    # Symbolic
    assert qapply(Jx*JxKet(j,m)) == hbar*m*JxKet(j,m)
    assert qapply(Jx*JzKet(j,m)) == \
        hbar*sqrt(j**2+j-m**2-m)*JzKet(j,m+1)/2 + hbar*sqrt(j**2+j-m**2+m)*JzKet(j,m-1)/2
    assert qapply(Jx*TensorProduct(JxKet(j1,m1), JxKet(j2,m2))) == \
        hbar*m1*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))+hbar*m2*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))
    assert qapply(Jx*TensorProduct(JzKet(j1,m1), JzKet(j2,m2))) == \
        hbar*sqrt(j1**2+j1-m1**2-m1)*TensorProduct(JzKet(j1,m1+1),JzKet(j2,m2))/2 + \
        hbar*sqrt(j1**2+j1-m1**2+m1)*TensorProduct(JzKet(j1,m1-1),JzKet(j2,m2))/2 + \
        hbar*sqrt(j2**2+j2-m2**2-m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2+1))/2 + \
        hbar*sqrt(j2**2+j2-m2**2+m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2-1))/2
    # Uncoupled Operators
    # Numerical
    assert qapply(TensorProduct(Jx,1)*TensorProduct(JxKet(1,1),JxKet(1,-1))) == hbar*TensorProduct(JxKet(1,1),JxKet(1,-1))
    assert qapply(TensorProduct(1,Jx)*TensorProduct(JxKet(1,1),JxKet(1,-1))) == -hbar*TensorProduct(JxKet(1,1),JxKet(1,-1))
    assert qapply(TensorProduct(Jx,1)*TensorProduct(JyKet(1,1),JyKet(1,-1))) == hbar*TensorProduct(JyKet(1,1),JyKet(1,-1))
    assert qapply(TensorProduct(1,Jx)*TensorProduct(JyKet(1,1),JyKet(1,-1))) == -hbar*TensorProduct(JyKet(1,1),JyKet(1,-1))
    assert qapply(TensorProduct(Jx,1)*TensorProduct(JzKet(1,1),JzKet(1,-1))) == hbar*sqrt(2)*TensorProduct(JzKet(1,0),JzKet(1,-1))/2
    assert qapply(TensorProduct(1,Jx)*TensorProduct(JzKet(1,1),JzKet(1,-1))) == hbar*sqrt(2)*TensorProduct(JzKet(1,1),JzKet(1,0))/2
    # Symbolic
    assert qapply(TensorProduct(Jx,1)*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))) == \
        hbar*m1*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))
    assert qapply(TensorProduct(1,Jx)*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))) == \
        hbar*m2*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))
    assert qapply(TensorProduct(Jx,1)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))) == \
        hbar*sqrt(j1**2+j1-m1**2-m1)*TensorProduct(JzKet(j1,m1+1),JzKet(j2,m2))/2 + hbar*sqrt(j1**2+j1-m1**2+m1)*TensorProduct(JzKet(j1,m1-1),JzKet(j2,m2))/2
    assert qapply(TensorProduct(1,Jx)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))) == \
        hbar*sqrt(j2**2+j2-m2**2-m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2+1))/2 + hbar*sqrt(j2**2+j2-m2**2+m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2-1))/2

@XFAIL
def test_jx_failing():
    assert qapply(Jx*JyKet(j,m)) == Sum(hbar*mi*WignerD(j,mi,m,0,0,pi/2)*JxKet(j,mi),(mi,-j,j))
    assert qapply(Jx*TensorProduct(JyKet(j1,m1), JyKet(j2,m2))) == \
                TensorProduct(Sum(hbar*mi*WignerD(j1,mi,m1,0,0,pi/2)*JxKet(j1,mi),(mi,-j1,j1)),JyKet(j2,m2)) + \
                TensorProduct(JyKet(j1,m1),Sum(hbar*mi*WignerD(j2,mi,m2,0,0,pi/2)*JxKet(j2,mi),(mi,-j2,j2)))
    assert qapply(TensorProduct(Jx,1)*TensorProduct(JyKet(j1,m1),JyKet(j2,m2))) == \
                TensorProduct(Sum(hbar*mi*WignerD(j1,mi,m1,0,0,pi/2) * JxKet(j1,mi), (mi,-j1,j1)),JyKet(j2,m2))
    assert qapply(TensorProduct(1,Jx)*TensorProduct(JyKet(j1,m1),JyKet(j2,m2))) == \
        TensorProduct(JyKet(j1,m1),Sum(hbar*mi*WignerD(j2,mi,m2,0,0,pi/2) * JxKet(j2,mi), (mi,-j2,j2)))

def test_jy():
    assert Commutator(Jy, Jz).doit() == I*hbar*Jx
    assert Jy.rewrite('plusminus') == (Jplus - Jminus)/(2*I)
    assert represent(Jy, basis=Jz) == (represent(Jplus, basis=Jz) - represent(Jminus, basis=Jz))/(2*I)
    # Normal Operators
    # Numerical
    assert qapply(Jy*JxKet(1,1)) == hbar*JxKet(1,1)
    assert qapply(Jy*JyKet(1,1)) == hbar*JyKet(1,1)
    assert qapply(Jy*JzKet(1,1)) == sqrt(2)*hbar*I*JzKet(1,0)/2
    assert qapply(Jy*TensorProduct(JxKet(1,1), JxKet(1,1))) == \
        hbar*TensorProduct(JxKet(1,1),JxKet(1,1)) + hbar*TensorProduct(JxKet(1,1),JxKet(1,1))
    assert qapply(Jy*TensorProduct(JyKet(1,1), JyKet(1,1))) == 2*hbar*TensorProduct(JyKet(1,1), JyKet(1,1))
    assert qapply(Jy*TensorProduct(JzKet(1,1), JzKet(1,1))) == \
        sqrt(2)*hbar*I*TensorProduct(JzKet(1,1),JzKet(1,0))/2+sqrt(2)*hbar*I*TensorProduct(JzKet(1,0),JzKet(1,1))/2
    assert qapply(Jy*TensorProduct(JyKet(1,1), JyKet(1,-1))) == 0
    # Symbolic
    assert qapply(Jy*JyKet(j,m)) == hbar*m*JyKet(j,m)
    assert qapply(Jy*JzKet(j,m)) == \
        -hbar*I*sqrt(j**2+j-m**2-m)*JzKet(j,m+1)/2 + hbar*I*sqrt(j**2+j-m**2+m)*JzKet(j,m-1)/2
    assert qapply(Jy*TensorProduct(JyKet(j1,m1), JyKet(j2,m2))) == \
        hbar*m1*TensorProduct(JyKet(j1,m1),JyKet(j2,m2))+hbar*m2*TensorProduct(JyKet(j1,m1),JyKet(j2,m2))
    assert qapply(Jy*TensorProduct(JzKet(j1,m1), JzKet(j2,m2))) == \
        -hbar*I*sqrt(j1**2+j1-m1**2-m1)*TensorProduct(JzKet(j1,m1+1),JzKet(j2,m2))/2 + \
        hbar*I*sqrt(j1**2+j1-m1**2+m1)*TensorProduct(JzKet(j1,m1-1),JzKet(j2,m2))/2 + \
        -hbar*I*sqrt(j2**2+j2-m2**2-m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2+1))/2 + \
        hbar*I*sqrt(j2**2+j2-m2**2+m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2-1))/2
    # Uncoupled Operators
    # Numerical
    assert qapply(TensorProduct(Jy,1)*TensorProduct(JxKet(1,1),JxKet(1,-1))) == hbar*TensorProduct(JxKet(1,1),JxKet(1,-1))
    assert qapply(TensorProduct(1,Jy)*TensorProduct(JxKet(1,1),JxKet(1,-1))) == -hbar*TensorProduct(JxKet(1,1),JxKet(1,-1))
    assert qapply(TensorProduct(Jy,1)*TensorProduct(JyKet(1,1),JyKet(1,-1))) == hbar*TensorProduct(JyKet(1,1),JyKet(1,-1))
    assert qapply(TensorProduct(1,Jy)*TensorProduct(JyKet(1,1),JyKet(1,-1))) == -hbar*TensorProduct(JyKet(1,1),JyKet(1,-1))
    assert qapply(TensorProduct(Jy,1)*TensorProduct(JzKet(1,1),JzKet(1,-1))) == hbar*sqrt(2)*I*TensorProduct(JzKet(1,0),JzKet(1,-1))/2
    assert qapply(TensorProduct(1,Jy)*TensorProduct(JzKet(1,1),JzKet(1,-1))) == -hbar*sqrt(2)*I*TensorProduct(JzKet(1,1),JzKet(1,0))/2
    # Symbolic
    assert qapply(TensorProduct(Jy,1)*TensorProduct(JyKet(j1,m1),JyKet(j2,m2))) == \
        hbar*m1*TensorProduct(JyKet(j1,m1), JyKet(j2,m2))
    assert qapply(TensorProduct(1,Jy)*TensorProduct(JyKet(j1,m1),JyKet(j2,m2))) == \
        hbar*m2*TensorProduct(JyKet(j1,m1), JyKet(j2,m2))
    assert qapply(TensorProduct(Jy,1)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))) == \
        -hbar*I*sqrt(j1**2+j1-m1**2-m1)*TensorProduct(JzKet(j1,m1+1),JzKet(j2,m2))/2 + hbar*I*sqrt(j1**2+j1-m1**2+m1)*TensorProduct(JzKet(j1,m1-1),JzKet(j2,m2))/2
    assert qapply(TensorProduct(1,Jy)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))) == \
        -hbar*I*sqrt(j2**2+j2-m2**2-m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2+1))/2 + hbar*I*sqrt(j2**2+j2-m2**2+m2)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2-1))/2

@XFAIL
def test_jy_failing():
    assert qapply(Jy*JxKet(j,m)) == Sum(hbar*mi*WignerD(j,mi,m,3*pi/2,0,0)*JyKet(j,mi), (mi,-j,j))
    assert qapply(Jy*TensorProduct(JxKet(j1,m1), JxKet(j2,m2))) == \
                    TensorProduct(JxKet(j1,m1),Sum(hbar*mi*WignerD(j2,mi,m2,3*pi/2,0,0)*JyKet(j2,mi),(mi,-j2,j2))) + \
                    TensorProduct(Sum(hbar*mi*WignerD(j1,mi,m1,3*pi/2,0,0)*JyKet(j1,mi),(mi,-j1,j1)),JxKet(j2,m2))
    assert qapply(TensorProduct(Jy,1)*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))) == \
                    TensorProduct(Sum(hbar*mi*WignerD(j1,mi,m1,3*pi/2,0,0) * JyKet(j1,mi), (mi,-j1,j1)), JxKet(j2,m2))
    assert qapply(TensorProduct(1,Jy)*TensorProduct(JxKet(j1,m1),JxKet(j2,m2))) == \
                    TensorProduct(JxKet(j1,m1), Sum(hbar*mi*WignerD(j2,mi,m2,3*pi/2,0,0) * JyKet(j2,mi), (mi,-j2,j2)))

def test_jz():
    assert Commutator(Jz, Jminus).doit() == -hbar*Jminus
    # Normal Operators
    # Numerical
    assert qapply(Jz*JxKet(1,1)) == -sqrt(2)*hbar*JxKet(1,0)/2
    assert qapply(Jz*JyKet(1,1)) == -sqrt(2)*hbar*I*JyKet(1,0)/2
    assert qapply(Jz*JzKet(2,1)) == hbar*JzKet(2,1)
    assert qapply(Jz*TensorProduct(JxKet(1,1), JxKet(1,1))) == \
        -sqrt(2)*hbar*TensorProduct(JxKet(1,1),JxKet(1,0))/2 - sqrt(2)*hbar*TensorProduct(JxKet(1,0),JxKet(1,1))/2
    assert qapply(Jz*TensorProduct(JyKet(1,1), JyKet(1,1))) == \
        -sqrt(2)*hbar*I*TensorProduct(JyKet(1,1), JyKet(1,0))/2 - sqrt(2)*hbar*I*TensorProduct(JyKet(1,0), JyKet(1,1))/2
    assert qapply(Jz*TensorProduct(JzKet(1,1), JzKet(1,1))) == 2*hbar*TensorProduct(JzKet(1,1), JzKet(1,1))
    assert qapply(Jz*TensorProduct(JzKet(1,1), JzKet(1,-1))) == 0
    # Symbolic
    assert qapply(Jz*JzKet(j,m)) == hbar*m*JzKet(j,m)
    assert qapply(Jz*TensorProduct(JzKet(j1,m1), JzKet(j2,m2))) == \
        hbar*m1*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))+hbar*m2*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))
    # Uncoupled Operators
    # Numerical
    assert qapply(TensorProduct(Jz,1)*TensorProduct(JxKet(1,1),JxKet(1,-1))) == -sqrt(2)*hbar*TensorProduct(JxKet(1,0),JxKet(1,-1))/2
    assert qapply(TensorProduct(1,Jz)*TensorProduct(JxKet(1,1),JxKet(1,-1))) == -sqrt(2)*hbar*TensorProduct(JxKet(1,1),JxKet(1,0))/2
    assert qapply(TensorProduct(Jz,1)*TensorProduct(JyKet(1,1),JyKet(1,-1))) == -sqrt(2)*I*hbar*TensorProduct(JyKet(1,0),JyKet(1,-1))/2
    assert qapply(TensorProduct(1,Jz)*TensorProduct(JyKet(1,1),JyKet(1,-1))) == sqrt(2)*I*hbar*TensorProduct(JyKet(1,1),JyKet(1,0))/2
    assert qapply(TensorProduct(Jz,1)*TensorProduct(JzKet(1,1),JzKet(1,-1))) == hbar*TensorProduct(JzKet(1,1),JzKet(1,-1))
    assert qapply(TensorProduct(1,Jz)*TensorProduct(JzKet(1,1),JzKet(1,-1))) == -hbar*TensorProduct(JzKet(1,1),JzKet(1,-1))
    # Symbolic
    assert qapply(TensorProduct(1,Jz)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))) == \
        hbar*m2*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))

@XFAIL
def test_jz_failing():
    assert qapply(Jz*JxKet(j,m)) == Sum(hbar*mi*WignerD(j,mi,m,0,pi/2,0)*JzKet(j,mi), (mi,-j,j))
    assert qapply(Jz*JyKet(j,m)) == Sum(hbar*mi*WignerD(j,mi,m,3*pi/2,-pi/2,pi/2)*JzKet(j,mi), (mi,-j,j))
    assert qapply(Jz*TensorProduct(JxKet(j1,m1), JxKet(j2,m2))) == \
        TensorProduct(JxKet(j1,m1),Sum(hbar*mi*WignerD(j2,mi,m2,0,pi/2,0)*JzKet(j2,mi),(mi,-j2,j2))) + \
        TensorProduct(Sum(hbar*mi*WignerD(j1,mi,m1,0,pi/2,0)*JzKet(j1,mi),(mi,-j1,j1)),JxKet(j2,m2))
    assert qapply(Jz*TensorProduct(JyKet(j1,m1), JyKet(j2,m2))) == \
        TensorProduct(JyKet(j1,m1),Sum(hbar*mi*WignerD(j2,mi,m2,3*pi/2,-pi/2,pi/2)*JzKet(j2,mi),(mi,-j2,j2))) + \
        TensorProduct(Sum(hbar*mi*WignerD(j1,mi,m1,3*pi/2,-pi/2,pi/2)*JzKet(j1,mi),(mi,-j1,j1)),JyKet(j2,m2))
    assert qapply(TensorProduct(Jz,1)*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))) == \
        hbar*m1*TensorProduct(JzKet(j1,m1),JzKet(j2,m2))
