"""Tests for cartesian.py"""

from sympy import S, Interval, symbols, I, DiracDelta, exp, sqrt, pi

from sympy.physics.quantum import qapply, represent, L2, Dagger
from sympy.physics.quantum import Commutator, hbar
from sympy.physics.quantum.cartesian import (
    XOp, PxOp, X, Px, XKet, XBra, PxKet, PxBra
)

x, y = symbols('x,y')
px, py = symbols('px py')


def test_x():
    assert X.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))
    assert Commutator(X, Px).doit() == I*hbar
    assert qapply(X*XKet(x)) == x*XKet(x)
    assert XKet(x).dual_class == XBra
    assert XBra(x).dual_class == XKet
    assert (Dagger(XKet(y))*XKet(x)).doit() == DiracDelta(x-y)
    assert (PxBra(px)*XKet(x)).doit() ==\
        exp(-I*x*px/hbar)/sqrt(2*pi*hbar)
    assert represent(XKet(x)) == x
    assert XBra(x).position == x


def test_p():
    assert Px.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))
    assert qapply(Px*PxKet(px)) == px*PxKet(px)
    assert PxKet(px).dual_class == PxBra
    assert PxBra(x).dual_class == PxKet
    assert (Dagger(PxKet(py))*PxKet(px)).doit() == DiracDelta(px-py)
    assert (XBra(x)*PxKet(px)).doit() ==\
        exp(I*x*px/hbar)/sqrt(2*pi*hbar)
    assert represent(PxKet(px)) == px
