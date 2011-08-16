"""Tests for cartesian.py"""

from sympy import (
    S, Interval, symbols, I, DiracDelta, exp, sqrt, pi, Derivative, Function
)

from sympy.physics.quantum import qapply, represent, L2, Dagger
from sympy.physics.quantum import Commutator, hbar
from sympy.physics.quantum.cartesian import (
    XOp, YOp, ZOp, PxOp, X, Y, Z, Px, XKet, XBra, PxKet, PxBra,
    PositionKet3D, PositionBra3D
)
from sympy.physics.quantum.operator import DifferentialOperator
from sympy.physics.quantum.state import Wavefunction

x, y, z, x_1, x_2, x_3, y_1, z_1 = symbols('x,y,z,x_1,x_2,x_3,y_1,z_1', real=True)
px, py, px_1, px_2 = symbols('px py px_1 px_2', real=True)
f = Function('f')


def test_x():
    assert X.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))
    assert Commutator(X, Px).doit() == I*hbar
    assert qapply(X*XKet(x)) == x*XKet(x)
    assert XKet(x).dual_class() == XBra
    assert XBra(x).dual_class() == XKet
    assert (Dagger(XKet(y))*XKet(x)).doit() == DiracDelta(x-y)
    assert (PxBra(px)*XKet(x)).doit() ==\
        exp(-I*x*px/hbar)/sqrt(2*pi*hbar)
    assert represent(XKet(x)) == DiracDelta(x - x_1)
    assert represent(XOp()) == DiracDelta(x - x_1)*Wavefunction(x, x)
    assert represent(XBra(x)) == DiracDelta(x - x_1)
    assert XBra(x).position == x
    assert represent(XOp()*XKet()) == DiracDelta(x-x_1)*Wavefunction(x, x)
    assert represent(XOp()*XKet()*XBra('y')) == \
           DiracDelta(x - x_2)*DiracDelta(y - x_1)*Wavefunction(x, x)
    assert represent(XBra("y")*XKet()) == DiracDelta(y - x)
    assert represent(XKet()*XBra()) == \
           DiracDelta(x - x_1)*DiracDelta(x - x_2)

    rep_p = represent(XOp(), basis = PxOp)
    diff_op1 = DifferentialOperator(Derivative(f(px), px), f(px))
    assert rep_p == hbar*I*DiracDelta(px - px_1)*diff_op1
    assert rep_p == represent(XOp(), basis = PxOp())
    assert rep_p == represent(XOp(), basis = PxKet)
    assert rep_p == represent(XOp(), basis = PxKet())

    #assert represent(XOp()*PxKet(), basis = PxKet) == \
    #       Wavefunction(-hbar*I*DiracDelta(px - px_2, 1), px)

def test_p():
    assert Px.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))
    assert qapply(Px*PxKet(px)) == px*PxKet(px)
    assert PxKet(px).dual_class() == PxBra
    assert PxBra(x).dual_class() == PxKet
    assert (Dagger(PxKet(py))*PxKet(px)).doit() == DiracDelta(px-py)
    assert (XBra(x)*PxKet(px)).doit() ==\
        exp(I*x*px/hbar)/sqrt(2*pi*hbar)
    assert represent(PxKet(px)) == DiracDelta(px - px_1)

    rep_x = represent(PxOp(), basis = XOp)
    diff_op1 = DifferentialOperator(Derivative(f(x), x), f(x))
    assert rep_x == -hbar*I*DiracDelta(x - x_1)*diff_op1
    assert rep_x == represent(PxOp(), basis = XOp())
    assert rep_x == represent(PxOp(), basis = XKet)
    assert rep_x == represent(PxOp(), basis = XKet())

    diff_op = DifferentialOperator(Derivative(f(x), x), f(x))
    #assert represent(PxOp()*XKet(), basis=XKet) == \
    #       Wavefunction(hbar*I*DiracDelta(x - x_2, 1), x)
    #assert represent(XBra("y")*PxOp()*XKet(), basis=XKet) == \
    #       Wavefunction(hbar*I*DiracDelta(x - y, 1), x, y)

def test_3dpos():
    assert Y.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))
    assert Z.hilbert_space == L2(Interval(S.NegativeInfinity, S.Infinity))

    test_ket = PositionKet3D(x, y, z)
    assert qapply(X*test_ket) == x*test_ket
    assert qapply(Y*test_ket) == y*test_ket
    assert qapply(Z*test_ket) == z*test_ket
    assert qapply(X*Y*test_ket) == x*y*test_ket
    assert qapply(X*Y*Z*test_ket) == x*y*z*test_ket
    assert qapply(Y*Z*test_ket) == y*z*test_ket

    assert PositionKet3D() == test_ket
    assert YOp() == Y
    assert ZOp() == Z

    assert PositionKet3D.dual_class() == PositionBra3D
    assert PositionBra3D.dual_class() == PositionKet3D

    other_ket = PositionKet3D(x_1, y_1, z_1)
    assert (Dagger(other_ket)*test_ket).doit() == \
           DiracDelta(x - x_1)*DiracDelta(y - y_1)*DiracDelta(z - z_1)

    assert test_ket.position_x == x
    assert test_ket.position_y == y
    assert test_ket.position_z == z
    assert other_ket.position_x == x_1
    assert other_ket.position_y == y_1
    assert other_ket.position_z == z_1

    # TODO: Add tests for representations
