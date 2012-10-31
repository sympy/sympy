"""Tests for sho1d.py"""

from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.qexpr import *
from sympy.physics.quantum.cartesian import *

from sympy.physics.quantum.sho1d import *

ap = RaisingOp('a')
am = LoweringOp('a')
k = SHOKet('k')
b = SHOBra('b')
H = Hamiltonian('H')
N = NumberOp('N')
w = Symbol('omega')
m = Symbol('m')

def test_ap():
	assert adjoint(ap) == am
	assert Commutator(ap, am).doit() == Integer(-1)
	assert Commutator(ap, N).doit() == Integer(-1)*ap
	assert qapply(ap*k) == sqrt(k.n + 1)*SHOKet(k.n + 1)
	assert ap().rewrite('xp').doit() == \
		(Integer(1)/sqrt(Integer(2)*hbar*m*w))*(Integer(-1)*I*Px + m*w*X)
	
def test_am():
	assert adjoint(am) == ap
	assert Commutator(am, ap).doit() == Integer(1)
	assert Commutator(am, N).doit() == am
	assert qapply(am*k) == sqrt(k.n)*SHOKet(k.n-Integer(1))
	assert am().rewrite('xp').doit() == \
		(Integer(1)/sqrt(Integer(2)*hbar*m*w))*(I*Px + m*w*X)
		
def test_k():
	assert SHOKet('k').dual_class() == SHOBra
	assert SHOBra('b').dual_class() == SHOKet
	assert InnerProduct(b,k).doit() == KroneckerDelta(k.n, b.n)
	
def test_N():
	assert Commutator(N, ap).doit() == ap
	assert Commutator(N, am).doit() == Integer(-1)*am
	assert Commutator(N, H).doit() == Integer(0)
	assert qapply(N*k) == k.n*k
	assert N().rewrite('a').doit() == ap*am
	assert N().rewrite('H').doit() == H/(hbar*w) - Integer(1)/Integer(2)
	
def test_H():
	assert Commutator(H, N).doit() == Integer(0)
	assert qapply(H*k) == (hbar*w*(k.n + Integer(1)/Integer(2)))*k
	assert H().rewrite('a').doit() == hbar*w*(ap*am + Integer(1)/Integer(2))
	assert H().rewrite('am').doit() == hbar*w*(am*ap - Integer(1)/Integer(2))
	assert H().rewrite('xp').doit() == \
		(Integer(1)/(Integer(2)*m))*(Px**2 + (m*w*X)**2)
	assert H().rewrite('n').doit() == hbar*w*(N + Integer(1)/Integer(2))
	