from sympy import I, symbols, Symbol, sqrt, expand
from sympy.utilities.pytest import XFAIL

from sympy.physics.quantum.operator import Operator
from sympy.physics.quantum.applyops import apply_operators as applyops
from sympy.physics.quantum.constants import hbar
from sympy.physics.quantum.spin import (
    Jx, Jy, Jz, Jplus, Jminus,
    JzKet, JzBra, JxKet, JxBra
)

j = Symbol('j')
m = Symbol('m')
jp = Symbol("j'")
mp = Symbol("m'")

z = JzKet(1,0)
po = JzKet(1,1)
mo = JzKet(1,-1)

A = Operator('A')

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

@XFAIL
def test_innerproduct():
    assert applyops(po.dual*Jz*po) == po*po
