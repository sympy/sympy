from sympy import Symbol
from sympy.physics.classical import Point, ReferenceFrame, Dyad, RigidBody

def test_rigidbody():
    m = Symbol('m')
    A = ReferenceFrame('A')
    P = Point('P')
    I = Dyad([])
    B = RigidBody()
    assert B.mass == None
    assert B.mc == None
    assert B.inertia == (None, None)
    assert B.frame == None

    B.mass = m
    B.frame = A
    B.cm = P
    B.inertia = (I, B.cm)
    assert B.mass == m
    assert B.frame == A
    assert B.cm == P
    assert B.inertia == (I, B.cm)
