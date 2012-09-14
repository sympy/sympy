from sympy import symbols, sin, cos
from sympy.physics.mechanics import (cross, dot, dynamicsymbols, express,
                                     ReferenceFrame, inertia, Point,
                                     kinematic_equations, Vector,
                                     inertia_of_point_mass, partial_velocity,
                                     outer, Particle,
                                     RigidBody, angular_momentum,
                                     linear_momentum, kinetic_energy,
                                     potential_energy)

Vector.simp = True
q1, q2, q3, q4, q5 = symbols('q1 q2 q3 q4 q5')
N = ReferenceFrame('N')
A = N.orientnew('A', 'Axis', [q1, N.z])
B = A.orientnew('B', 'Axis', [q2, A.x])
C = B.orientnew('C', 'Axis', [q3, B.y])

def test_dot():
    assert dot(A.x, A.x) == 1
    assert dot(A.x, A.y) == 0
    assert dot(A.x, A.z) == 0

    assert dot(A.y, A.x) == 0
    assert dot(A.y, A.y) == 1
    assert dot(A.y, A.z) == 0

    assert dot(A.z, A.x) == 0
    assert dot(A.z, A.y) == 0
    assert dot(A.z, A.z) == 1

# TODO: Add dot product tests from different frames

def test_cross():
    assert cross(A.x, A.x) == 0
    assert cross(A.x, A.y) == A.z
    assert cross(A.x, A.z) == -A.y

    assert cross(A.y, A.x) == -A.z
    assert cross(A.y, A.y) == 0
    assert cross(A.y, A.z) == A.x

    assert cross(A.z, A.x) == A.y
    assert cross(A.z, A.y) == -A.x
    assert cross(A.z, A.z) == 0

def test_cross_different_frames():
    assert cross(N.x, A.x) == sin(q1)*A.z
    assert cross(N.x, A.y) == cos(q1)*A.z
    assert cross(N.x, A.z) == -sin(q1)*A.x-cos(q1)*A.y
    assert cross(N.y, A.x) == -cos(q1)*A.z
    assert cross(N.y, A.y) == sin(q1)*A.z
    assert cross(N.y, A.z) == cos(q1)*A.x-sin(q1)*A.y
    assert cross(N.z, A.x) == A.y
    assert cross(N.z, A.y) == -A.x
    assert cross(N.z, A.z) == 0

    assert cross(N.x, A.x) == sin(q1)*A.z
    assert cross(N.x, A.y) == cos(q1)*A.z
    assert cross(N.x, A.x + A.y) == sin(q1)*A.z + cos(q1)*A.z
    assert cross(A.x + A.y, N.x) == -sin(q1)*A.z - cos(q1)*A.z

    assert cross(A.x, C.x) == sin(q3)*C.y
    assert cross(A.x, C.y) == -sin(q3)*C.x + cos(q3)*C.z
    assert cross(A.x, C.z) == -cos(q3)*C.y
    assert cross(C.x, A.x) == -sin(q3)*C.y
    assert cross(C.y, A.x) == sin(q3)*C.x - cos(q3)*C.z
    assert cross(C.z, A.x) == cos(q3)*C.y

def test_express():
    assert express(A.x, C) == cos(q3)*C.x + sin(q3)*C.z
    assert express(A.y, C) == sin(q2)*sin(q3)*C.x + cos(q2)*C.y - \
            sin(q2)*cos(q3)*C.z
    assert express(A.z, C) == -sin(q3)*cos(q2)*C.x + sin(q2)*C.y + \
            cos(q2)*cos(q3)*C.z
    assert express(A.x, N) == cos(q1)*N.x + sin(q1)*N.y
    assert express(A.y, N) == -sin(q1)*N.x + cos(q1)*N.y
    assert express(A.z, N) == N.z
    assert express(A.x, A) == A.x
    assert express(A.y, A) == A.y
    assert express(A.z, A) == A.z
    assert express(A.x, B) == B.x
    assert express(A.y, B) == cos(q2)*B.y - sin(q2)*B.z
    assert express(A.z, B) == sin(q2)*B.y + cos(q2)*B.z
    assert express(A.x, C) == cos(q3)*C.x + sin(q3)*C.z
    assert express(A.y, C) == sin(q2)*sin(q3)*C.x + cos(q2)*C.y - \
            sin(q2)*cos(q3)*C.z
    assert express(A.z, C) == -sin(q3)*cos(q2)*C.x + sin(q2)*C.y + \
            cos(q2)*cos(q3)*C.z
    # Check to make sure UnitVectors get converted properly
    assert express(N.x, N) == N.x
    assert express(N.y, N) == N.y
    assert express(N.z, N) == N.z
    assert express(N.x, A) == (cos(q1)*A.x - sin(q1)*A.y)
    assert express(N.y, A) == (sin(q1)*A.x + cos(q1)*A.y)
    assert express(N.z, A) == A.z
    assert express(N.x, B) == (cos(q1)*B.x - sin(q1)*cos(q2)*B.y + \
            sin(q1)*sin(q2)*B.z)
    assert express(N.y, B) == (sin(q1)*B.x + cos(q1)*cos(q2)*B.y - \
            sin(q2)*cos(q1)*B.z)
    assert express(N.z, B) == (sin(q2)*B.y + cos(q2)*B.z)
    assert express(N.x, C) == (
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*C.x -
            sin(q1)*cos(q2)*C.y +
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*C.z)
    assert express(N.y, C) == (
            (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*C.x +
            cos(q1)*cos(q2)*C.y +
            (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*C.z)
    assert express(N.z, C) == (-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z)

    assert express(A.x, N) == (cos(q1)*N.x + sin(q1)*N.y)
    assert express(A.y, N) == (-sin(q1)*N.x + cos(q1)*N.y)
    assert express(A.z, N) == N.z
    assert express(A.x, A) == A.x
    assert express(A.y, A) == A.y
    assert express(A.z, A) == A.z
    assert express(A.x, B) == B.x
    assert express(A.y, B) == (cos(q2)*B.y - sin(q2)*B.z)
    assert express(A.z, B) == (sin(q2)*B.y + cos(q2)*B.z)
    assert express(A.x, C) == (cos(q3)*C.x + sin(q3)*C.z)
    assert express(A.y, C) == (sin(q2)*sin(q3)*C.x + cos(q2)*C.y -
            sin(q2)*cos(q3)*C.z)
    assert express(A.z, C) == (-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z)

    assert express(B.x, N) == (cos(q1)*N.x + sin(q1)*N.y)
    assert express(B.y, N) == (-sin(q1)*cos(q2)*N.x +
            cos(q1)*cos(q2)*N.y + sin(q2)*N.z)
    assert express(B.z, N) == (sin(q1)*sin(q2)*N.x -
            sin(q2)*cos(q1)*N.y + cos(q2)*N.z)
    assert express(B.x, A) == A.x
    assert express(B.y, A) == (cos(q2)*A.y + sin(q2)*A.z)
    assert express(B.z, A) == (-sin(q2)*A.y + cos(q2)*A.z)
    assert express(B.x, B) == B.x
    assert express(B.y, B) == B.y
    assert express(B.z, B) == B.z
    assert express(B.x, C) == (cos(q3)*C.x + sin(q3)*C.z)
    assert express(B.y, C) == C.y
    assert express(B.z, C) == (-sin(q3)*C.x + cos(q3)*C.z)

    assert express(C.x, N) == (
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*N.x +
            (sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*N.y -
                sin(q3)*cos(q2)*N.z)
    assert express(C.y, N) == (
            -sin(q1)*cos(q2)*N.x + cos(q1)*cos(q2)*N.y + sin(q2)*N.z)
    assert express(C.z, N) == (
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*N.x +
            (sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))*N.y +
            cos(q2)*cos(q3)*N.z)
    assert express(C.x, A) == (cos(q3)*A.x + sin(q2)*sin(q3)*A.y -
            sin(q3)*cos(q2)*A.z)
    assert express(C.y, A) == (cos(q2)*A.y + sin(q2)*A.z)
    assert express(C.z, A) == (sin(q3)*A.x - sin(q2)*cos(q3)*A.y +
            cos(q2)*cos(q3)*A.z)
    assert express(C.x, B) == (cos(q3)*B.x - sin(q3)*B.z)
    assert express(C.y, B) == B.y
    assert express(C.z, B) == (sin(q3)*B.x + cos(q3)*B.z)
    assert express(C.x, C) == C.x
    assert express(C.y, C) == C.y
    assert express(C.z, C) == C.z == (C.z)

    #  Check to make sure Vectors get converted back to UnitVectors
    assert N.x == express((cos(q1)*A.x - sin(q1)*A.y), N)
    assert N.y == express((sin(q1)*A.x + cos(q1)*A.y), N)
    assert N.x == express((cos(q1)*B.x - sin(q1)*cos(q2)*B.y +
            sin(q1)*sin(q2)*B.z), N)
    assert N.y == express((sin(q1)*B.x + cos(q1)*cos(q2)*B.y -
        sin(q2)*cos(q1)*B.z), N)
    assert N.z == express((sin(q2)*B.y + cos(q2)*B.z), N)

    """
    These don't really test our code, they instead test the auto simplification
    (or lack thereof) of SymPy.
    assert N.x == express((
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*C.x -
            sin(q1)*cos(q2)*C.y +
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*C.z), N)
    assert N.y == express((
            (sin(q1)*cos(q3) + sin(q2)*sin(q3)*cos(q1))*C.x +
            cos(q1)*cos(q2)*C.y +
            (sin(q1)*sin(q3) - sin(q2)*cos(q1)*cos(q3))*C.z), N)
    assert N.z == express((-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z), N)
    """

    assert A.x == express((cos(q1)*N.x + sin(q1)*N.y), A)
    assert A.y == express((-sin(q1)*N.x + cos(q1)*N.y), A)

    assert A.y == express((cos(q2)*B.y - sin(q2)*B.z), A)
    assert A.z == express((sin(q2)*B.y + cos(q2)*B.z), A)

    assert A.x == express((cos(q3)*C.x + sin(q3)*C.z), A)

    # Tripsimp messes up here too.
    #print express((sin(q2)*sin(q3)*C.x + cos(q2)*C.y -
    #        sin(q2)*cos(q3)*C.z), A)
    assert A.y == express((sin(q2)*sin(q3)*C.x + cos(q2)*C.y -
            sin(q2)*cos(q3)*C.z), A)

    assert A.z == express((-sin(q3)*cos(q2)*C.x + sin(q2)*C.y +
            cos(q2)*cos(q3)*C.z), A)
    assert B.x == express((cos(q1)*N.x + sin(q1)*N.y), B)
    assert B.y == express((-sin(q1)*cos(q2)*N.x +
            cos(q1)*cos(q2)*N.y + sin(q2)*N.z), B)

    assert B.z == express((sin(q1)*sin(q2)*N.x -
            sin(q2)*cos(q1)*N.y + cos(q2)*N.z), B)

    assert B.y == express((cos(q2)*A.y + sin(q2)*A.z), B)
    assert B.z == express((-sin(q2)*A.y + cos(q2)*A.z), B)
    assert B.x == express((cos(q3)*C.x + sin(q3)*C.z), B)
    assert B.z == express((-sin(q3)*C.x + cos(q3)*C.z), B)

    """
    assert C.x == express((
            (cos(q1)*cos(q3)-sin(q1)*sin(q2)*sin(q3))*N.x +
            (sin(q1)*cos(q3)+sin(q2)*sin(q3)*cos(q1))*N.y -
                sin(q3)*cos(q2)*N.z), C)
    assert C.y == express((
            -sin(q1)*cos(q2)*N.x + cos(q1)*cos(q2)*N.y + sin(q2)*N.z), C)
    assert C.z == express((
            (sin(q3)*cos(q1)+sin(q1)*sin(q2)*cos(q3))*N.x +
            (sin(q1)*sin(q3)-sin(q2)*cos(q1)*cos(q3))*N.y +
            cos(q2)*cos(q3)*N.z), C)
    """
    assert C.x == express((cos(q3)*A.x + sin(q2)*sin(q3)*A.y -
            sin(q3)*cos(q2)*A.z), C)
    assert C.y == express((cos(q2)*A.y + sin(q2)*A.z), C)
    assert C.z == express((sin(q3)*A.x - sin(q2)*cos(q3)*A.y +
            cos(q2)*cos(q3)*A.z), C)
    assert C.x == express((cos(q3)*B.x - sin(q3)*B.z), C)
    assert C.z == express((sin(q3)*B.x + cos(q3)*B.z), C)

def test_inertia():
    N = ReferenceFrame('N')
    ixx, iyy, izz = symbols('ixx iyy izz')
    ixy, iyz, izx = symbols('ixy iyz izx')
    assert inertia(N, ixx, iyy, izz) == (ixx * (N.x | N.x) + iyy *
            (N.y | N.y) + izz * (N.z | N.z))
    assert inertia(N, 0, 0, 0) == 0 * (N.x | N.x)
    assert inertia(N, ixx, iyy, izz, ixy, iyz, izx) == (ixx * (N.x | N.x) +
            ixy * (N.x | N.y) + izx * (N.x | N.z) + ixy * (N.y | N.x) + iyy *
            (N.y | N.y) + iyz * (N.y | N.z) + izx * (N.z | N.x) + iyz * (N.z |
            N.y) + izz * (N.z | N.z))

def test_kin_eqs():
    q0, q1, q2, q3 = dynamicsymbols('q0 q1 q2 q3')
    q0d, q1d, q2d, q3d = dynamicsymbols('q0 q1 q2 q3', 1)
    u1, u2, u3 = dynamicsymbols('u1 u2 u3')
    kds = kinematic_equations([u1, u2, u3], [q0, q1, q2, q3], 'quaternion')
    assert kds == [-0.5 * q0 * u1 - 0.5 * q2 * u3 + 0.5 * q3 * u2 + q1d,
            -0.5 * q0 * u2 + 0.5 * q1 * u3 - 0.5 * q3 * u1 + q2d,
            -0.5 * q0 * u3 - 0.5 * q1 * u2 + 0.5 * q2 * u1 + q3d,
            0.5 * q1 * u1 + 0.5 * q2 * u2 + 0.5 * q3 * u3 + q0d]

def test_inertia_of_point_mass():
    r, s, t, m = symbols('r s t m')
    N = ReferenceFrame('N')

    px = r * N.x
    I = inertia_of_point_mass(m, px, N)
    assert I == m * r**2 * (N.y | N.y) + m * r**2 * (N.z | N.z)

    py = s * N.y
    I = inertia_of_point_mass(m, py, N)
    assert I == m * s**2 * (N.x | N.x) + m * s**2 * (N.z | N.z)

    pz = t * N.z
    I = inertia_of_point_mass(m, pz, N)
    assert I == m * t**2 * (N.x | N.x) + m * t**2 * (N.y | N.y)

    p = px + py + pz
    I = inertia_of_point_mass(m, p, N)
    assert I == (m * (s**2 + t**2) * (N.x | N.x) -
                 m * r * s * (N.x | N.y) -
                 m * r * t * (N.x | N.z) -
                 m * r * s * (N.y | N.x) +
                 m * (r**2 + t**2) * (N.y | N.y) -
                 m * s * t * (N.y | N.z) -
                 m * r * t * (N.z | N.x) -
                 m * s * t * (N.z | N.y) +
                 m * (r**2 + s**2) * (N.z | N.z))

def test_partial_velocity():
    q1, q2, q3, u1, u2, u3  = dynamicsymbols('q1 q2 q3 u1 u2 u3')
    u4, u5 = dynamicsymbols('u4, u5')
    r = symbols('r')

    N = ReferenceFrame('N')
    Y = N.orientnew('Y', 'Axis', [q1, N.z])
    L = Y.orientnew('L', 'Axis', [q2, Y.x])
    R = L.orientnew('R', 'Axis', [q3, L.y])
    R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)

    C = Point('C')
    C.set_vel(N, u4 * L.x + u5 * (Y.z ^ L.x))
    Dmc = C.locatenew('Dmc', r * L.z)
    Dmc.v2pt_theory(C, N, R)

    vel_list = [Dmc.vel(N), C.vel(N), R.ang_vel_in(N)]
    u_list = [u1, u2, u3, u4, u5]
    assert (partial_velocity(vel_list, u_list, N) ==
            [[- r*L.y, r*L.x, 0, L.x, cos(q2)*L.y - sin(q2)*L.z],
            [0, 0, 0, L.x, cos(q2)*L.y - sin(q2)*L.z],
            [L.x, L.y, L.z, 0, 0]])

def test_linear_momentum():
    N = ReferenceFrame('N')
    Ac = Point('Ac')
    Ac.set_vel(N, 25 * N.y)
    I = outer(N.x, N.x)
    A = RigidBody('A', Ac, N, 20, (I, Ac))
    P = Point('P')
    Pa = Particle('Pa', P, 1)
    Pa.point.set_vel(N, 10 * N.x)
    assert linear_momentum(N, A, Pa) == 10 * N.x + 500 * N.y

def test_angular_momentum_and_linear_momentum():
    m, M, l1 = symbols('m M l1')
    q1d = dynamicsymbols('q1d')
    N = ReferenceFrame('N')
    O = Point('O')
    O.set_vel(N, 0 * N.x)
    Ac = O.locatenew('Ac', l1 * N.x)
    P = Ac.locatenew('P', l1 * N.x)
    a = ReferenceFrame('a')
    a.set_ang_vel(N, q1d * N.z)
    Ac.v2pt_theory(O, N, a)
    P.v2pt_theory(O, N, a)
    Pa = Particle('Pa', P, m)
    I = outer(N.z, N.z)
    A = RigidBody('A', Ac, a, M, (I, Ac))
    assert linear_momentum(N, A, Pa) == 2 * m * q1d* l1 * N.y + M * l1 * q1d * N.y
    assert angular_momentum(O, N, A, Pa) == 4 * m * q1d * l1**2 * N.z + q1d * N.z

def test_kinetic_energy():
    m, M, l1 = symbols('m M l1')
    omega = dynamicsymbols('omega')
    N = ReferenceFrame('N')
    O = Point('O')
    O.set_vel(N, 0 * N.x)
    Ac = O.locatenew('Ac', l1 * N.x)
    P = Ac.locatenew('P', l1 * N.x)
    a = ReferenceFrame('a')
    a.set_ang_vel(N, omega * N.z)
    Ac.v2pt_theory(O, N, a)
    P.v2pt_theory(O, N, a)
    Pa = Particle('Pa', P, m)
    I = outer(N.z, N.z)
    A = RigidBody('A', Ac, a, M, (I, Ac))
    assert 0 == kinetic_energy(N, Pa, A) - (M*l1**2*omega**2/2
            + 2*l1**2*m*omega**2 + omega**2/2)

def test_potential_energy():
    m, M, l1, g, h, H = symbols('m M l1 g h H')
    omega = dynamicsymbols('omega')
    N = ReferenceFrame('N')
    O = Point('O')
    O.set_vel(N, 0 * N.x)
    Ac = O.locatenew('Ac', l1 * N.x)
    P = Ac.locatenew('P', l1 * N.x)
    a = ReferenceFrame('a')
    a.set_ang_vel(N, omega * N.z)
    Ac.v2pt_theory(O, N, a)
    P.v2pt_theory(O, N, a)
    Pa = Particle('Pa', P, m)
    I = outer(N.z, N.z)
    A = RigidBody('A', Ac, a, M, (I, Ac))
    Pa.set_potential_energy(m * g * h)
    A.set_potential_energy(M * g * H)
    assert potential_energy(A, Pa) == m * g * h + M * g * H
