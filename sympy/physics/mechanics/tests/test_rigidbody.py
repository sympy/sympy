from sympy import symbols
from sympy.physics.mechanics import Point, ReferenceFrame, Dyadic, RigidBody
from sympy.physics.mechanics import dynamicsymbols, outer

def test_rigidbody():
    m, m2, v1,v2, v3, omega = symbols('m m2 v1 v2 v3 omega')
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
    N =  ReferenceFrame('N')
    P2.set_vel(N, v1 * N.x + v2 * N.y + v3 * N.z)
    assert B.linear_momentum(N) == m2 * (v1 * N.x + v2 * N.y + v3 * N.z)

def test_rigidbody2():
    M, v, r, omega, g, h = dynamicsymbols('M v r omega g h')
    N = ReferenceFrame('N')
    b = ReferenceFrame('b')
    b.set_ang_vel(N, omega * b.x)
    P = Point('P')
    I = outer (b.x, b.x)
    Inertia_tuple = (I, P)
    B = RigidBody('B', P, b, M, Inertia_tuple)
    P.set_vel(N, v * b.x)
    assert B.angular_momentum(P, N) == omega * b.x
    O = Point('O')
    O.set_vel(N, v * b.x)
    P.set_pos(O, r * b.y)
    assert B.angular_momentum(O, N) == omega * b.x - M*v*r*b.z
    B.set_potential_energy(M * g * h)
    assert B.potential_energy == M * g * h
    assert B.kinetic_energy(N) == (omega**2 + M * v**2) / 2
