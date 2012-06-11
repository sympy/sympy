from sympy import symbols
from sympy.physics.mechanics import Point, Particle, ReferenceFrame

def test_particle():
    m, m2, v1, v2, v3 = symbols('m m2 v1 v2 v3')
    P = Point('P')
    P2 = Point('P2')
    p = Particle('pa', P, m)
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
    P2.set_vel(N, v1 * N.x)
    assert p.linmom(N) == m2 * v1 * N.x
    P2.set_vel(N, v2 * N.y)
    assert p.linmom(N) == m2 * v2 * N.y
    P2.set_vel(N, v3 * N.z)
    assert p.linmom(N) == m2 * v3 * N.z
    P2.set_vel(N, v1 * N.x + v2 * N.y + v3 * N.z)
    assert p.linmom(N) == m2 * (v1 * N.x + v2 * N.y + v3 * N.z)
