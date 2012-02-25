from sympy import symbols
from sympy.physics.mechanics import Point, Particle

def test_particle():
    m, m2 = symbols('m m2')
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
