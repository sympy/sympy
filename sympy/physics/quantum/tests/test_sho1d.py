"""Tests for sho1d.py"""

from sympy import Integer, Symbol, sqrt, I
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum import Commutator, adjoint
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.cartesian import X, Px
from sympy.physics.quantum.special.tensor_functions import KroneckerDelta

from sympy.physics.quantum.sho1d import (RaisingOp, LoweringOp,
										SHOKet, SHOBra, 
										Hamiltonian, NumberOp)

ad = RaisingOp('a')
a = LoweringOp('a')
k = SHOKet('k')
kz = SHOKet('0')
b = SHOBra('b')
H = Hamiltonian('H')
N = NumberOp('N')
omega = Symbol('omega')
m = Symbol('m')

def test_ad():
	assert adjoint(ad) == a
	assert Commutator(ad, a).doit() == Integer(-1)
	assert Commutator(ad, N).doit() == Integer(-1)*ad
	assert qapply(ad*k) == (sqrt(k.n + 1)*SHOKet(k.n + 1)).expand()
	assert ad().rewrite('xp').doit() == \
		(Integer(1)/sqrt(Integer(2)*hbar*m*omega))*(Integer(-1)*I*Px + m*omega*X)
	
def test_a():
	assert adjoint(a) == ad
	assert Commutator(a, ad).doit() == Integer(1)
	assert Commutator(a, N).doit() == a
	assert qapply(a*k) == (sqrt(k.n)*SHOKet(k.n-Integer(1))).expand()
	assert qapply(a*kz) == Integer(0)
	assert a().rewrite('xp').doit() == \
		(Integer(1)/sqrt(Integer(2)*hbar*m*omega))*(I*Px + m*omega*X)
		
def test_k():
	assert SHOKet('k').dual_class() == SHOBra
	assert SHOBra('b').dual_class() == SHOKet
	assert InnerProduct(b,k).doit() == KroneckerDelta(k.n, b.n)
	
def test_N():
	assert Commutator(N, ad).doit() == ad
	assert Commutator(N, a).doit() == Integer(-1)*a
	assert Commutator(N, H).doit() == Integer(0)
	assert qapply(N*k) == (k.n*k).expand()
	assert N().rewrite('a').doit() == ad*a
	assert N().rewrite('H').doit() == H/(hbar*omega) - Integer(1)/Integer(2)
	
def test_H():
	assert Commutator(H, N).doit() == Integer(0)
	assert qapply(H*k) == ((hbar*omega*(k.n + Integer(1)/Integer(2)))*k).expand()
	assert H().rewrite('a').doit() == hbar*omega*(ad*a + Integer(1)/Integer(2))
	assert H().rewrite('xp').doit() == \
		(Integer(1)/(Integer(2)*m))*(Px**2 + (m*omega*X)**2)
	assert H().rewrite('n').doit() == hbar*omega*(N + Integer(1)/Integer(2))
	