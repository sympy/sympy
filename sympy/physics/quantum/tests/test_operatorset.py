from sympy.physics.quantum.operatorset import operator_to_state, state_to_operator
from sympy.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet

def test_op_to_state():
    assert operator_to_state(XOp) == XKet
    assert operator_to_state(PxOp) == PxKet

    exception = False
    try:
        operator_to_state(XKet)
    except NotImplementedError:
        exception = True

    return exception

def test_state_to_op():
    assert state_to_operator(XKet) == XOp
    assert state_to_operator(PxKet) == PxOp

    exception = False
    try:
        state_to_operator(XOp)
    except NotImplementedError:
        exception = True

    return exception
