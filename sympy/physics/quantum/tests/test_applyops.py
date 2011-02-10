from sympy import I, symbols, Symbol, sqrt, expand, Integer, srepr

from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.applyops import apply_operators as applyops
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.spin import (
    Jx, Jy, Jz, Jplus, Jminus, J2,
    JzKet, JzBra, JxKet, JxBra
)
from sympy.physics.quantum.gate import H
from sympy.physics.quantum.qubit import Qubit


j = Symbol('j')
m = Symbol('m')
jp = Symbol("j'")
mp = Symbol("m'")

z = JzKet(1,0)
po = JzKet(1,1)
mo = JzKet(1,-1)

A = Operator('A')


class Foo(Operator):

    def _apply_operator_JzKet(self, ket, **options):
        return ket


def test_basic():
    assert applyops(Jz*po) == hbar*po
    assert applyops(Jx*z) == hbar*po/sqrt(2) + hbar*mo/sqrt(2)
    assert applyops((Jplus + Jminus)*z/sqrt(2)) == hbar*po + hbar*mo
    assert applyops(Jz*(po + mo)) == hbar*po - hbar*mo
    assert applyops(Jz*po + Jz*mo) == hbar*po - hbar*mo
    assert applyops(Jminus*Jminus*po) == 2*hbar**2*mo
    assert applyops(Jplus**2*mo) == 2*hbar**2*po
    assert applyops(Jplus**2*Jminus**2*po) == 4*hbar**4*po


def test_extra():
    extra = z.dual*A*z
    assert applyops(Jz*po*extra) == hbar*po*extra
    assert applyops(Jx*z*extra) == (hbar*po/sqrt(2) + hbar*mo/sqrt(2))*extra
    assert applyops((Jplus + Jminus)*z/sqrt(2)*extra) == hbar*po*extra + hbar*mo*extra
    assert applyops(Jz*(po + mo)*extra) == hbar*po*extra - hbar*mo*extra
    assert applyops(Jz*po*extra + Jz*mo*extra) == hbar*po*extra - hbar*mo*extra
    assert applyops(Jminus*Jminus*po*extra) == 2*hbar**2*mo*extra
    assert applyops(Jplus**2*mo*extra) == 2*hbar**2*po*extra
    assert applyops(Jplus**2*Jminus**2*po*extra) == 4*hbar**4*po*extra


def test_innerproduct():
    assert applyops(po.dual*Jz*po, ip_doit=False) == hbar*(po.dual*po)
    assert applyops(po.dual*Jz*po) == hbar


def test_zero():
    assert applyops(0) == 0
    assert applyops(Integer(0)) == 0


def test_commutator():
    assert applyops(Commutator(Jx,Jy)*Jz*po) == I*hbar**3*po
    assert applyops(Commutator(J2, Jz)*Jz*po) == 0
    assert applyops(Commutator(Jz, Foo('F'))*po) == 0
    assert applyops(Commutator(Foo('F'), Jz)*po) == 0


def test_anticommutator():
    assert applyops(AntiCommutator(Jz, Foo('F'))*po) == 2*hbar*po
    assert applyops(AntiCommutator(Foo('F'), Jz)*po) == 2*hbar*po


def test_outerproduct():
    e = Jz*(mo*po.dual)*Jz*po
    assert applyops(e) == -hbar**2*mo
    assert applyops(e, ip_doit=False) == -hbar**2*(po.dual*po)*mo
    assert applyops(e).doit() == -hbar**2*mo


def test_dagger():
    lhs = Dagger(Qubit(0))*Dagger(H(0))
    rhs = Dagger(Qubit(1))/sqrt(2) + Dagger(Qubit(0))/sqrt(2)
    assert applyops(lhs, dagger=True) == rhs
