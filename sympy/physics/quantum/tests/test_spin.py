from __future__ import division
from sympy import cos, exp, I, Matrix, pi, S, sin, sqrt, symbols

from sympy.physics.quantum import hbar, represent, Commutator, InnerProduct
from sympy.physics.quantum import qapply
from sympy.physics.quantum.spin import (
    Jx, Jy, Jz, Jplus, Jminus, J2,
    JxBra, JyBra, JzBra,
    JxKet, JyKet, JzKet,
    Rotation, WignerD
)

from sympy.utilities.pytest import XFAIL


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
    assert represent(JxKet(S(1)/2,S(1)/2)) == sqrt(2)*Matrix([1,1])/2
    assert represent(JxKet(S(1)/2,-S(1)/2)) == sqrt(2)*Matrix([-1,1])/2
    assert represent(JxKet(1,1)) == Matrix([1,sqrt(2),1])/2
    assert represent(JxKet(1,0)) == sqrt(2)*Matrix([-1,0,1])/2
    assert represent(JxKet(1,-1)) == Matrix([1,-sqrt(2),1])/2
    assert represent(JyKet(S(1)/2,S(1)/2)) == sqrt(2)*Matrix([-1,-I])/2
    assert represent(JyKet(S(1)/2,-S(1)/2)) == sqrt(2)*Matrix([-I,-1])/2
    assert represent(JyKet(1,1)) == Matrix([1,sqrt(2)*I,-1])/2
    assert represent(JyKet(1,0)) == sqrt(2)*Matrix([I,0,I])/2
    assert represent(JyKet(1,-1)) == Matrix([-1,sqrt(2)*I,1])/2
    assert represent(JzKet(S(1)/2,S(1)/2)) == Matrix([1,0])
    assert represent(JzKet(S(1)/2,-S(1)/2)) == Matrix([0,1])
    assert represent(JzKet(1,1)) == Matrix([1,0,0])
    assert represent(JzKet(1,0)) == Matrix([0,1,0])
    assert represent(JzKet(1,-1)) == Matrix([0,0,1])

def test_rewrite():
    # Rewrite to same basis
    assert JxBra(1,1).rewrite('Jx') == JxBra(1,1)
    assert JxKet(1,1).rewrite('Jx') == JxKet(1,1)
    # Rewriting to different basis
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
    # Innerproducts of rewritten states
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

def test_innerproduct():
    j,m = symbols("j m")
    assert InnerProduct(JzBra(1,1), JzKet(1,1)).doit() == 1
    assert InnerProduct(JzBra(S(1)/2,S(1)/2), JzKet(S(1)/2,-S(1)/2)).doit() == 0
    assert InnerProduct(JzBra(j,m), JzKet(j,m)).doit() == 1
    assert InnerProduct(JzBra(1,0), JyKet(1,1)).doit() == I/sqrt(2)
    assert InnerProduct(JxBra(S(1)/2,S(1)/2), JzKet(S(1)/2,S(1)/2)).doit() == -sqrt(2)/2
    assert InnerProduct(JyBra(1,1), JzKet(1,1)).doit() == S(1)/2
    assert InnerProduct(JxBra(1,-1), JyKet(1,1)).doit() == 0

def test_rotation_small_d():
    # Symbolic tests
    beta = symbols('beta')
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
    alpha, beta, gamma = symbols('alpha beta gamma')
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
    j, m, mp, alpha, beta, gamma = symbols('j m mp alpha beta gamma')
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
    j, m = symbols('j m')
    assert Commutator(J2, Jz).doit() == 0
    assert qapply(J2*JzKet(1,1)) == 2*hbar**2*JzKet(1,1)
    assert qapply(J2*JzKet(j,m)) == j**2*hbar**2*JzKet(j,m)+j*hbar**2*JzKet(j,m)
    assert J2.matrix_element(1,1,1,1) == 2*hbar**2

def test_jx():
    assert Commutator(Jx, Jz).doit() == -I*hbar*Jy
    assert qapply(Jx*JzKet(1,1)) == sqrt(2)*hbar*JzKet(1,0)/2
    assert Jx.rewrite('plusminus') == (Jminus + Jplus)/2
    assert represent(Jx, basis=Jz, j=1) == (represent(Jplus, basis=Jz, j=1)+represent(Jminus, basis=Jz, j=1))/2

def test_jy():
    assert Commutator(Jy, Jz).doit() == I*hbar*Jx
    assert qapply(Jy*JzKet(1,1)) == I*sqrt(2)*hbar*JzKet(1,0)/2
    assert Jy.rewrite('plusminus') == (Jplus - Jminus)/(2*I)
    assert represent(Jy, basis=Jz) == (represent(Jplus, basis=Jz) - represent(Jminus, basis=Jz))/(2*I)

def test_jz():
    assert Commutator(Jz, Jminus).doit() == -hbar*Jminus
    assert qapply(Jz*JzKet(2,1)) == hbar*JzKet(2,1)
    assert Jz.rewrite('plusminus')
