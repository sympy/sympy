from sympy import Matrix
from sympy.physics.units.units import PREFIXES, Unit, UnitSystem
from sympy.physics.units.mks import length, time, velocity

def test_prefix_operations():
    m = PREFIXES['m']
    k = PREFIXES['k']
    M = PREFIXES['M']
    assert m*k == 1
    assert k*k == M
    assert 1/m == k
    assert k/m == M

def test_unitsystem():
    matrix = Matrix(((1,0),(1,-1)))
    m = Unit(abbrev='m', dimension=length)
    v = Unit(abbrev='v', dimension=velocity)
    us = UnitSystem(base=(m, v))
    assert us._dim_matrix == matrix

