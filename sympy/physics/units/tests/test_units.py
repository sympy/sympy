from sympy import Matrix, eye
from sympy.physics.units.units import PREFIXES, Unit, UnitSystem
from sympy.physics.units.mks import mks, m, v, length, time, velocity
from sympy.utilities.pytest import raises

def test_prefix_operations():
    m = PREFIXES['m']
    k = PREFIXES['k']
    M = PREFIXES['M']

    assert m*k == 1
    assert k*k == M
    assert 1/m == k
    assert k/m == M

def test_unit_operations():
    m2 = m*m

    assert m2.dimension == m.dimension**2
    assert m2.factor == m.factor**2

def test_unitsystem():
    raises(AttributeError, 'UnitSystem(base=(m, m))')
    assert mks._transf_matrix == eye(3)

    matrix = Matrix(((1,1),(0,-1)))
    us = UnitSystem(base=(m, v))
    assert us._transf_matrix == matrix
