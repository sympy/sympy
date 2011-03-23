from sympy import I, Matrix

from sympy.physics.quantum import hbar, represent, Commutator
from sympy.physics.quantum import qapply
from sympy.physics.quantum.spin import *


def test_represent():
    assert represent(Jz) == hbar*Matrix([[1,0],[0,-1]])/2
    assert represent(Jz, j=1) == hbar*Matrix([[1,0,0],[0,0,0],[0,0,-1]])

def test_jplus():
    assert Commutator(Jplus, Jminus).doit() == 2*hbar*Jz
    assert qapply(Jplus*JzKet(1,1)) == 0
    assert Jplus.matrix_element(1,1,1,1) == 0
    assert Jplus.rewrite('xyz') == Jx + I*Jy

def test_rotation():
    assert Rotation.d(1,1,1,0) == 1
