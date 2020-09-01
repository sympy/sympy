from sympy import symbols
from sympy.physics.mechanics import Point, Particle, ReferenceFrame, inertia
from sympy.physics.units.definitions.unit_definitions import G

from sympy.testing.pytest import raises


def test_particle():
    m, m1, m2, v1, v2, v3, r, g, h = symbols('m m1 m2 v1 v2 v3 r g h')
    P = Point('P')
    P2 = Point('P2')
    p = Particle('pa', P, m)
    assert p.__str__() == 'pa'
    assert p.mass == m
    assert p.point == P
    # Test the mass setter
    p.mass = m2
    assert p.mass == m2
    # Test the point setter
    p.point = P2
    assert p.point == P2
    # Test the linear momentum function
    N = ReferenceFrame('N')
    O = Point('O')
    P2.set_pos(O, r * N.y)
    P2.set_vel(N, v1 * N.x)
    raises(TypeError, lambda: Particle(P, P, m))
    raises(TypeError, lambda: Particle('pa', m, m))
    assert p.linear_momentum(N) == m2 * v1 * N.x
    assert p.angular_momentum(O, N) == -m2 * r *v1 * N.z
    P2.set_vel(N, v2 * N.y)
    assert p.linear_momentum(N) == m2 * v2 * N.y
    assert p.angular_momentum(O, N) == 0
    P2.set_vel(N, v3 * N.z)
    assert p.linear_momentum(N) == m2 * v3 * N.z
    assert p.angular_momentum(O, N) == m2 * r * v3 * N.x
    P2.set_vel(N, v1 * N.x + v2 * N.y + v3 * N.z)
    assert p.linear_momentum(N) == m2 * (v1 * N.x + v2 * N.y + v3 * N.z)
    assert p.angular_momentum(O, N) == m2 * r * (v3 * N.x - v1 * N.z)
    p.potential_energy = m * g * h
    assert p.potential_energy == m * g * h
    # TODO make the result not be system-dependent
    assert p.kinetic_energy(
        N) in [m2*(v1**2 + v2**2 + v3**2)/2,
        m2 * v1**2 / 2 + m2 * v2**2 / 2 + m2 * v3**2 / 2]
    # Test the Gravitational Force function
    O1 = Point('O1')
    P1 = Particle('P1', O1, m1)
    O2 = Point('O2')
    O2.set_pos(O1, 10 * N.x)
    P2 = Particle('P2', O2, m2)
    assert P1.G_Force(P2) == G * m1 * m2 / 100 * N.x
    O3 = Point('O3')
    O3.set_pos(O1, -10 * N.x)
    P3 = Particle('P3', O3, m2)
    assert P1.G_Force(P2, P3) == 0
    raises(TypeError, lambda: P1.G_Force(r))

def test_parallel_axis():
    N = ReferenceFrame('N')
    m, a, b = symbols('m, a, b')
    o = Point('o')
    p = o.locatenew('p', a * N.x + b * N.y)
    P = Particle('P', o, m)
    Ip = P.parallel_axis(p, N)
    Ip_expected = inertia(N, m * b**2, m * a**2, m * (a**2 + b**2),
                          ixy=-m * a * b)
    assert Ip == Ip_expected
