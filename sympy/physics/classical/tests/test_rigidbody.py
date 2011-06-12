from sympy import Symbol
from sympy.physics.classical import (Point, ReferenceFrame, InertiaDyadic,
        RigidBody)

def test_rigidbody():
    m = Symbol('m')
    A = ReferenceFrame('A')
    P = Point('P')
    I = InertiaDyadic()
    B = RigidBody()
    assert B.mass == None
    assert B.mc == None
    assert B.inertia == None
    assert B.frame == None

    B.mass = m
    B.frame = A
    B.cm = P
    B.inertia = I
    assert B.mass == m
    assert B.frame == A
    assert B.cm == P
    assert B.inertia == I
