# -*- coding: utf-8 -*-
from sympy import Rational
from sympy.physics.units.definitions import (m, s, c, kg)
from sympy.physics.units.dimensions import Dimension, DimensionSystem, length, time, mass, velocity, current, \
    action
from sympy.physics.units.unitsystem import UnitSystem
from sympy.physics.units.quantities import Quantity
from sympy.utilities.pytest import raises


def test_definition():
    # want to test if the system can have several units of the same dimension
    dm = Quantity("dm", length, Rational(1, 10))

    base = (m, s)
    base_dim = (m.dimension, s.dimension)
    ms = UnitSystem(base, (c, dm), "MS", "MS system")

    assert set(ms._base_units) == set(base)
    assert set(ms._units) == set((m, s, c, dm))
    #assert ms._units == DimensionSystem._sort_dims(base + (velocity,))
    assert ms.name == "MS"
    assert ms.descr == "MS system"

    assert ms._system._base_dims == DimensionSystem.sort_dims(base_dim)
    assert set(ms._system._dims) == set(base_dim + (velocity,))


def test_error_definition():
    raises(ValueError, lambda: UnitSystem((m, s, c)))


def test_str_repr():
    assert str(UnitSystem((m, s), name="MS")) == "MS"
    assert str(UnitSystem((m, s))) == "UnitSystem((meter, second))"

    assert repr(UnitSystem((m, s))) == "<UnitSystem: (%s, %s)>" % (m, s)


def test_print_unit_base():
    A = Quantity("A", current, 1)
    Js = Quantity("Js", action, 1)
    mksa = UnitSystem((m, kg, s, A), (Js,))

    assert mksa.print_unit_base(Js) == m**2*kg*s**-1/1000


def test_extend():
    ms = UnitSystem((m, s), (c,))
    Js = Quantity("Js", action, 1)
    mks = ms.extend((kg,), (Js,))

    res = UnitSystem((m, s, kg), (c, Js))
    assert set(mks._base_units) == set(res._base_units)
    assert set(mks._units) == set(res._units)


def test_dim():
    dimsys = UnitSystem((m, kg, s), (c,))
    assert dimsys.dim == 3


def test_is_consistent():
    assert UnitSystem((m, s)).is_consistent is True
