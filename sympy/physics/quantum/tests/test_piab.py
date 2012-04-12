"""Tests for piab.py"""

from sympy import DiracDelta, I, Interval, pi, S, sin, sqrt, Symbol, symbols
from sympy.functions.special.tensor_functions import KroneckerDelta
from sympy.physics.quantum import L2, qapply, hbar, represent
from sympy.physics.quantum.piab import PIABHamiltonian, PIABKet, PIABBra, m, L
from sympy.physics.quantum.cartesian import XOp, XKet
from sympy.physics.quantum.state import Wavefunction

x = Symbol('x', real=True)
n = Symbol('n', integer=True)
i, j = symbols('i,j')

def test_H():
    assert PIABHamiltonian('H').hilbert_space ==\
        L2(Interval(S.NegativeInfinity,S.Infinity))
    assert qapply(PIABHamiltonian('H')*PIABKet(n)) ==\
        (n**2*pi**2*hbar**2)/(2*m*L**2)*PIABKet(n)

def test_states():
    assert PIABKet(n).dual_class() == PIABBra
    assert PIABKet(n).hilbert_space ==\
        L2(Interval(S.NegativeInfinity,S.Infinity))
    assert represent(PIABKet(n)) == \
           Wavefunction(sqrt(2/L)*sin(n*pi*x/L), (x, 0, L))
    assert represent(PIABKet(n), basis=XKet) == \
           Wavefunction(sqrt(2/L)*sin(n*pi*x/L), (x, 0, L))
    assert represent(PIABKet(n), basis=XOp) == \
           Wavefunction(sqrt(2/L)*sin(n*pi*x/L), (x, 0, L))
    assert (PIABBra(i)*PIABKet(j)).doit() == KroneckerDelta(i, j)
    assert PIABBra(n).dual_class() == PIABKet
