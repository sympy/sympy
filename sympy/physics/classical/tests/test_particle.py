from sympy import Symbol
from sympy.physics.classical import Point, Particle

def test_particle():
    m = Symbol('m')
    P = Point('P')
    p = Particle()
    assert p.mass == None
    assert p.point == None
    # Test the mass setter
    p.mass = m
    assert p.mass == m
    # Test the point setter
    p.point = P
    assert p.point == P
