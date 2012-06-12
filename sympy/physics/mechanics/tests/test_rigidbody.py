from sympy import symbols
from sympy.physics.mechanics import Point, ReferenceFrame, Dyadic, RigidBody

def test_rigidbody():
    m, m2, v1,v2, v3 = symbols('m m2 v1 v2 v3')
    A = ReferenceFrame('A')
    A2 = ReferenceFrame('A2')
    P = Point('P')
    P2 = Point('P2')
    I = Dyadic([])
    I2 = Dyadic([])
    B = RigidBody('B', P, A, m, (I, P))
    assert B.mass == m
    assert B.frame == A
    assert B.masscenter == P
    assert B.inertia == (I, B.masscenter)

    B.mass = m2
    B.frame = A2
    B.masscenter = P2
    B.inertia = (I2, B.masscenter)
    assert B.mass == m2
    assert B.frame == A2
    assert B.masscenter == P2
    assert B.inertia == (I2, B.masscenter)
    assert B.masscenter == P2
    assert B.inertia == (I2, B.masscenter)

    # Testing linear momentum function assuming A2 is the inertial frame
    P2.set_vel(A2, v1 * A2.x + v2 * A2.y + v3 * A2.z)
    assert B.linmom(A2) == m2 * (v1 * A2.x + v2 * A2.y + v3 * A2.z)
