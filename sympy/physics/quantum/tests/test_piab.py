"""Tests for piab.py"""

from sympy import S, Interval, symbols, I, DiracDelta, exp, sqrt, pi, sin

from sympy.physics.quantum import L2, qapply, hbar, represent
from sympy.physics.quantum import KroneckerDelta
from sympy.physics.quantum.piab import PIABHamiltonian, PIABKet, PIABBra, m, L


x, n = symbols('x,n')
i, j = symbols('i,j')


def test_H():
    assert PIABHamiltonian('H').hilbert_space ==\
        L2(Interval(S.NegativeInfinity,S.Infinity))
    assert qapply(PIABHamiltonian('H')*PIABKet(n)) ==\
        (n**2*pi**2*hbar**2)/(2*m*L**2)*PIABKet(n)


def test_states():
    assert PIABKet(n).dual_class == PIABBra
    assert PIABKet(n).hilbert_space ==\
        L2(Interval(S.NegativeInfinity,S.Infinity))
    assert represent(PIABKet(n)) == sqrt(2/L)*sin(n*pi*x/L)
    assert (PIABBra(i)*PIABKet(j)).doit() == KroneckerDelta(i, j)
    assert PIABBra(n).dual_class == PIABKet
