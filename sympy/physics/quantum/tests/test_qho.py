"""Tests for qho.py"""

from sympy import S, Interval, symbols, I, DiracDelta, exp, sqrt, pi, sin, Symbol

from sympy.physics.quantum import L2, qapply, hbar, represent
from sympy.physics.quantum import KroneckerDelta
from sympy.physics.quantum.qho import QHOHamiltonian, QHOBra, QHOKet, m, w
from sympy.physics.quantum.cartesian import XOp, XKet
from sympy.physics.quantum.state import Wavefunction


x = Symbol('x', real=True)
n = Symbol('n', integer=True)
i, j = symbols('i,j')


def test_H():
    assert QHOHamiltonian('H').hilbert_space ==\
        L2(Interval(S.NegativeInfinity,S.Infinity))
    assert qapply(QHOHamiltonian('H')*QHOKet(n)) ==\
        hbar*w*n*QHOKet(n) + hbar*w*(S(1)/2)*QHOKet(n)


def test_states():
    assert QHOKet(n).dual_class() == QHOBra
    assert QHOKet(n).hilbert_space ==\
        L2(Interval(S.NegativeInfinity,S.Infinity))
