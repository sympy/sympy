from sympy import symbols
from sympy.physics.mechanics import Point, ReferenceFrame, Dyadic, RigidBody

def test_rigidbody():
    m, m2 = symbols('m m2')
    A = ReferenceFrame('A')
    A2 = ReferenceFrame('A2')
    P = Point('P')
    P2 = Point('P2')
    I = Dyadic([])
    I2 = Dyadic([])
    B = RigidBody('B', P, A, m, (I, P))
    assert B.mass == m
    assert B.frame == A
    assert B.mc == P
    assert B.inertia == (I, B.mc)

    B.mass = m2
    B.frame = A2
    B.mc = P2
    B.inertia = (I2, B.mc)
    assert B.mass == m2
    assert B.frame == A2
    assert B.mc == P2
    assert B.inertia == (I2, B.mc)
