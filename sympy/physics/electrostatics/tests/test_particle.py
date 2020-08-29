from sympy import symbols
from sympy.physics.electrostatics import Point, Point_charge
from sympy import pi
from sympy.physics.units.definitions.unit_definitions import e0

from sympy.testing.pytest import raises

def test_particle():
    q1, q2, r, E = symbols('q1 q2 r E')
    P = Point('P')
    p = Point_charge('pa', P, q1)
    S = Point('S')
    s = Point_charge('pa2', S, q2)
    assert p.__str__() == 'pa'
    assert p.charge == q1
    assert p.point == P
    assert p.electric_potential(r,Er=E) == q1/(4*pi*e0*E*r)
    assert p.electric_potential(r) == q1/(4*pi*e0*r)
    assert p.electrostatic_force(s, r) == q1*q2/(4*pi*e0*r**2)
    assert p.electrostatic_force(s, r, val = True) == 9000000000.0*q1*q2/r**2
    assert p.electrostatic_force(s, r, Er=E) == q1*q2/(4*pi*e0*E*r**2)
    assert p.electrostatic_force(s, r, Er=E, val = True) == 9000000000*q1*q2/(E*r**2)
    assert p.electric_field_intensity(r) == q1/(4*pi*e0*r**2)
    assert p.electric_field_intensity(r,val=True) == 9000000000.0*q1/r**2
    assert p.electric_field_intensity(r,E) == q1/(4*pi*e0*E*r**2)
    assert p.electric_field_intensity(r,E, val =True) == 9000000000*q1/(E*r**2)
