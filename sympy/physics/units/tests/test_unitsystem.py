from sympy.physics.units import DimensionSystem, joule, second, ampere, electronvolt, convert_to, coulomb
from sympy.physics.units.systems import cgs
from sympy.utilities.pytest import warns_deprecated_sympy

from sympy import Rational, S
from sympy.physics.units.definitions import c, kg, m, s
from sympy.physics.units.definitions.dimension_definitions import (
    action, current, length, mass, time, velocity)
from sympy.physics.units.quantities import Quantity
from sympy.physics.units.unitsystem import UnitSystem
from sympy.utilities.pytest import raises


def test_definition():
    # want to test if the system can have several units of the same dimension
    dm = Quantity("dm")
    base = (m, s)
    base_dim = (m.dimension, s.dimension)
    ms = UnitSystem(base, (c, dm), "MS", "MS system")
    ms.set_quantity_dimension(dm, length)
    ms.set_quantity_scale_factor(dm, Rational(1, 10))

    assert set(ms._base_units) == set(base)
    assert set(ms._units) == {m, s, c, dm}
    # assert ms._units == DimensionSystem._sort_dims(base + (velocity,))
    assert ms.name == "MS"
    assert ms.descr == "MS system"


def test_str_repr():
    assert str(UnitSystem((m, s), name="MS")) == "MS"
    assert str(UnitSystem((m, s))) == "UnitSystem((meter, second))"

    assert repr(UnitSystem((m, s))) == "<UnitSystem: (%s, %s)>" % (m, s)


def test_print_unit_base():
    A = Quantity("A")
    A.set_global_relative_scale_factor(S.One, ampere)

    Js = Quantity("Js")
    Js.set_global_relative_scale_factor(S.One, joule*second)

    mksa = UnitSystem((m, kg, s, A), (Js,))
    with warns_deprecated_sympy():
        assert mksa.print_unit_base(Js) == m**2*kg*s**-1/1000


def test_extend():
    ms = UnitSystem((m, s), (c,))
    Js = Quantity("Js")
    Js.set_global_relative_scale_factor(1, joule*second)
    mks = ms.extend((kg,), (Js,))

    res = UnitSystem((m, s, kg), (c, Js))
    assert set(mks._base_units) == set(res._base_units)
    assert set(mks._units) == set(res._units)


def test_dim():
    dimsys = UnitSystem((m, kg, s), (c,))
    assert dimsys.dim == 3


def test_is_consistent():
    dimension_system = DimensionSystem([length, time])
    us = UnitSystem([m, s], dimension_system=dimension_system)
    assert us.is_consistent == True
