from sympy import Matrix, eye, Symbol, solve, pi
from sympy.physics.units.units import (PREFIXES, Unit, UnitSystem, set_system,
                                       Quantity as Q)
from sympy.physics.units.mks import mks, m, s, kg, J, v, length, time, velocity
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

    assert m**2 == m2

    assert m/m == 1
    assert m2/m == m
    assert m**2/m == m

    kg2 = kg**2
    assert kg2.factor == 10**6

    #assert m+m == m
    #assert m-m == m
    #raises(TypeError, 'm+1')

    assert J == kg * m**2 * s**-2
    assert J == kg * m**2 / s**2

    set_system(mks)
    assert m**2/m is m
    assert kg * m**2 * s**-2 is J
    assert kg * m**2 / s**2 is J
    set_system(None)

    x = Symbol('x')
    u = solve(x * kg - m, x)[0]
    assert u == m/kg

def test_unit_prop():
    set_system(mks)
    assert str(m/s*kg**2) == '1000 * kg**2 m s**-1'
    set_system(None)

def test_unitsystem():
    raises(AttributeError, 'UnitSystem(base=(m, m))')
    assert mks._transf_matrix == eye(3)

    matrix = Matrix(((1,1),(0,-1)))
    us = UnitSystem(base=(m, v))
    assert us._transf_matrix == matrix

    assert mks.base_factor == 10**3

def test_def_unitsystem():
    raises(TypeError, 'set_system(1)')

def test_quantity_operations():

    q1 = Q(6, m)
    q2 = Q(2, m)
    q3 = Q(12, m**2)

    assert q1 == 3*q2
    assert q1 == q2*3

    assert -q1 == Q(-6, m)
    assert q1 + q2 == Q(8, m)
    assert q1 + m == Q(7, m)
    assert q1 - q2 == Q(4, m)
    assert q1 - m == Q(5, m)
    assert q2 - q1 == Q(-4, m)
    assert q1 * q2 == q3
    assert q1 * m == Q(6, m**2)
    assert q1 / q2 == 3
    assert q1 / m == 6
    assert q3 / q1 == q2
    assert q2**2 == Q(4, m**2)
    assert q1**-1 == 1/q1

    raises(TypeError, 'q1+1')
    raises(TypeError, 'q1-1')
    raises(TypeError, '1+q1')
    raises(TypeError, '1-q1')
    raises(ValueError, 'q1+q3')

    # compute the radius using Kepler's law
    a = Symbol('a')
    T = 3.155815e7 * s
    M = 1.988435e30 * kg
    G = 6.67428 * 10**-11 * m**3 / kg / s**2

    aa = solve(T**2/a**3 - 4*pi**2 / G / M, a)[0]
    assert str(aa.simplify()) == ('149598206033.591*kg**(1/3)*kg**-1**(1/3)'
                             '*m**3**-1**-1**(1/3)'
                             '*s**2**(1/3)*s**2**-1**(1/3)')
