from sympy.physics.quantum.operatorset import operators_to_state, state_to_operators
from sympy.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet

def test_op_to_state():
    assert operators_to_state(XOp) == XKet
    assert operators_to_state(PxOp) == PxKet

    exception = False
    try:
        operators_to_state(XKet)
    except NotImplementedError:
        exception = True

    return exception

def test_state_to_op():
    assert state_to_operators(XKet) == XOp
    assert state_to_operators(PxKet) == PxOp

    exception = False
    try:
        state_to_operators(XOp)
    except NotImplementedError:
        exception = True

    return exception
