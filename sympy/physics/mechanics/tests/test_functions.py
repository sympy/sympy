from sympy import sin, cos, tan, pi, symbols, Matrix
from sympy.physics.mechanics import (Particle, Point, ReferenceFrame,
                                     RigidBody, Vector)
from sympy.physics.mechanics import (angular_momentum, dynamicsymbols,
                                     inertia, inertia_of_point_mass,
                                     kinetic_energy, linear_momentum,
                                     outer, potential_energy, msubs,
                                     find_dynamicsymbols)

Vector.simp = True
q1, q2, q3, q4, q5 = symbols('q1 q2 q3 q4 q5')
N = ReferenceFrame('N')
A = N.orientnew('A', 'Axis', [q1, N.z])
B = A.orientnew('B', 'Axis', [q2, A.x])
C = B.orientnew('C', 'Axis', [q3, B.y])


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
    """A rod with length 2l, centroidal inertia I, and mass M along with a
    particle of mass m fixed to the end of the rod rotate with an angular rate
    of omega about point O which is fixed to the non-particle end of the rod.
    The rod's reference frame is A and the inertial frame is N."""
    m, M, l, I = symbols('m, M, l, I')
    omega = dynamicsymbols('omega')
    N = ReferenceFrame('N')
    a = ReferenceFrame('a')
    O = Point('O')
    Ac = O.locatenew('Ac', l * N.x)
    P = Ac.locatenew('P', l * N.x)
    O.set_vel(N, 0 * N.x)
    a.set_ang_vel(N, omega * N.z)
    Ac.v2pt_theory(O, N, a)
    P.v2pt_theory(O, N, a)
    Pa = Particle('Pa', P, m)
    A = RigidBody('A', Ac, a, M, (I * outer(N.z, N.z), Ac))
    expected = 2 * m * omega * l * N.y + M * l * omega * N.y
    assert (linear_momentum(N, A, Pa) - expected) == Vector(0)
    expected = (I + M * l**2 + 4 * m * l**2) * omega * N.z
    assert (angular_momentum(O, N, A, Pa) - expected).simplify() == Vector(0)


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
    Pa.potential_energy = m * g * h
    A.potential_energy = M * g * H
    assert potential_energy(A, Pa) == m * g * h + M * g * H


def test_msubs():
    a, b = symbols('a, b')
    x, y, z = dynamicsymbols('x, y, z')
    # Test simple substitution
    expr = Matrix([[a*x + b, x*y.diff() + y],
                   [x.diff().diff(), z + sin(z.diff())]])
    sol = Matrix([[a + b, y],
                  [x.diff().diff(), 1]])
    sd = {x: 1, z: 1, z.diff(): 0, y.diff(): 0}
    assert msubs(expr, sd) == sol
    # Test smart substitution
    expr = cos(x + y)*tan(x + y) + b*x.diff()
    sd = {x: 0, y: pi/2, x.diff(): 1}
    assert msubs(expr, sd, smart=True) == b + 1
    N = ReferenceFrame('N')
    v = x*N.x + y*N.y
    d = x*(N.x|N.x) + y*(N.y|N.y)
    v_sol = 1*N.y
    d_sol = 1*(N.y|N.y)
    sd = {x: 0, y: 1}
    assert msubs(v, sd) == v_sol
    assert msubs(d, sd) == d_sol


def test_find_dynamicsymbols():
    a, b = symbols('a, b')
    x, y, z = dynamicsymbols('x, y, z')
    expr = Matrix([[a*x + b, x*y.diff() + y],
                   [x.diff().diff(), z + sin(z.diff())]])
    # Test finding all dynamicsymbols
    sol = {x, y.diff(), y, x.diff().diff(), z, z.diff()}
    assert find_dynamicsymbols(expr) == sol
    # Test finding all but those in sym_list
    exclude = [x, y, z]
    sol = {y.diff(), x.diff().diff(), z.diff()}
    assert find_dynamicsymbols(expr, exclude) == sol
