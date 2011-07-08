"""Tests for cartesian.py"""

from sympy import S, Interval, symbols, I, DiracDelta, exp, sqrt, pi

from sympy.physics.quantum import qapply, represent, L2, Dagger
from sympy.physics.quantum import Commutator, hbar
from sympy.physics.quantum.cartesian import (
    XOp, PxOp, X, Px, XKet, XBra, PxKet, PxBra
)
from sympy.physics.quantum.operator import DifferentialOperator

x, y, x_1, x_2, x_3 = symbols('x,y,x_1,x_2,x_3')
px, py, px_1, px_2 = symbols('px py px_1 px_2')


def test_x():
    assert X.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))
    assert Commutator(X, Px).doit() == I*hbar
    assert qapply(X*XKet(x)) == x*XKet(x)
    assert XKet(x).dual_class() == XBra
    assert XBra(x).dual_class() == XKet
    assert (Dagger(XKet(y))*XKet(x)).doit() == DiracDelta(x-y)
    assert (PxBra(px)*XKet(x)).doit() ==\
        exp(-I*x*px/hbar)/sqrt(2*pi*hbar)
    assert represent(XKet(x)) == DiracDelta(x-x_1)
    assert represent(XBra(x)) == DiracDelta(-x + x_1)
    assert XBra(x).position == x
    assert represent(XOp()*XKet()) == x*DiracDelta(x-x_2)
    assert represent(XOp()*XKet()*XBra('y')) == x*DiracDelta(x - x_3)*DiracDelta(x_1 - y)
    assert represent(XBra("y")*XKet()) == DiracDelta(x - y)
    assert represent(XKet()*XBra()) == DiracDelta(x - x_2) * DiracDelta(x_1 - x)

    rep_p = represent(XOp(), basis = PxOp)
    assert rep_p == hbar*I*DiracDelta(px_1 - px_2)*DifferentialOperator(px_1)
    assert rep_p == represent(XOp(), basis = PxOp())
    assert rep_p == represent(XOp(), basis = PxKet)
    assert rep_p == represent(XOp(), basis = PxKet())

    assert represent(XOp()*PxKet(), basis = PxKet) == hbar*I*DiracDelta(px - px_2)*DifferentialOperator(px)

def test_p():
    assert Px.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))
    assert qapply(Px*PxKet(px)) == px*PxKet(px)
    assert PxKet(px).dual_class() == PxBra
    assert PxBra(x).dual_class() == PxKet
    assert (Dagger(PxKet(py))*PxKet(px)).doit() == DiracDelta(px-py)
    assert (XBra(x)*PxKet(px)).doit() ==\
        exp(I*x*px/hbar)/sqrt(2*pi*hbar)
    assert represent(PxKet(px)) == DiracDelta(px-px_1)

    rep_x = represent(PxOp(), basis = XOp)
    assert rep_x == -hbar*I*DiracDelta(x_1 - x_2)*DifferentialOperator(x_1)
    assert rep_x == represent(PxOp(), basis = XOp())
    assert rep_x == represent(PxOp(), basis = XKet)
    assert rep_x == represent(PxOp(), basis = XKet())

    assert represent(PxOp()*XKet(), basis=XKet) == -hbar*I*DiracDelta(x - x_2)*DifferentialOperator(x)
    assert represent(XBra("y")*PxOp()*XKet(), basis=XKet) == -hbar*I*DiracDelta(x-y)*DifferentialOperator(x)
