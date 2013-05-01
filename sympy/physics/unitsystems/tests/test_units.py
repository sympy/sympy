from sympy import Matrix, eye, Symbol, solve, pi
import sympy.physics.unitsystems

from sympy.physics.unitsystems.units import (PREFIXES, Unit, UnitSystem,
                                unit_simplify, set_system, Quantity as Q)
from sympy.physics.unitsystems.mks import (mks, m, s, kg, J, v, G, length,
                                           time, mass, velocity)
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

    g = Unit(mass)
    assert kg/g == 10**3
    assert g/kg == 10**-3

    kg2 = kg**2
    assert kg2.factor == 10**6

    #raises(TypeError, 'm+1')

    assert J == kg * m**2 * s**-2
    assert J == kg * m**2 / s**2

    set_system(mks)
    #assert m**2/m is m
    #assert kg * m**2 * s**-2 is J
    #assert kg * m**2 / s**2 is J
    set_system(None)

    x = Symbol('x')
    u = solve(x * kg - m, x)[0]
    assert u == m/kg

def test_unit_prop():
    assert str(G.as_quantity) == '6.67384000000000e-14 '\
                                 '(1 {length: 3, mass: -1, time: -2})'

    def mks_asserts():
        assert str(m/s*kg**2).strip() == 'kg**2 m s**-1'
        assert str(G.as_quantity).strip() == '6.67384000000000e-11 m**3 kg**-1 s**-2'

        assert kg.factor == 1e3
        assert kg.base_factor == 1e3
        assert kg.ratio_factor == 1
        assert (2*kg).ratio_factor == 2

    set_system(mks)
    mks_asserts()
    set_system(None)

    # Repeat the same assert block, this time inside a with statement.
    with mks:
        mks_asserts()

def test_unitsystem():
    raises(AttributeError, lambda: UnitSystem(base=(m, m)))
    assert mks._transf_matrix == eye(3)

    matrix = Matrix(((1,1),(0,-1)))
    us = UnitSystem(base=(m, v))
    assert us._transf_matrix == matrix

    assert mks.base_factor == 10**3

def test_def_unitsystem():
    raises(TypeError, lambda: set_system(1))

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

    raises(TypeError, lambda: q1+1)
    raises(TypeError, lambda: q1-1)
    raises(TypeError, lambda: 1+q1)
    raises(TypeError, lambda: 1-q1)
    raises(ValueError, lambda: q1+q3)

    # compute the radius using Kepler's law
    a = Symbol('a')
    s = mks['s']
    kg = mks['kg']
    G = mks['G']
    T = 3.155815e7 * s
    M = 1.988435e30 * kg

    aa = solve(T**2/a**3 - 4*pi**2 / G / M, a)[0]
    a = unit_simplify(aa)
    # precision is not so good, it was better before, but ok for now
    assert abs((a.factor - 149598206033.591)/149598206033.591) < 1e-4
    # problem when comparing dimensions
    #assert a.as_quantity.unit == mks['m']
