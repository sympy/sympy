# -*- coding: utf-8 -*-

from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem
from sympy.physics.unitsystems.units import Unit, UnitSystem
from sympy.physics.unitsystems.quantities import Quantity
from sympy.utilities.pytest import raises

length = Dimension(name="length", symbol="L", length=1)
mass = Dimension(name="mass", symbol="M", mass=1)
time = Dimension(name="time", symbol="T", time=1)
current = Dimension(name="current", symbol="I", current=1)
velocity = Dimension(name="velocity", symbol="V", length=1, time=-1)
action = Dimension(name="action", symbol="A", length=2, mass=1, time=-2)

m = Unit(length, abbrev="m")
s = Unit(time, abbrev="s")
kg = Unit(mass, factor=10**3, abbrev="kg")
c = Unit(velocity, abbrev="c")


def test_definition():
    # want to test if the system can have several units of the same dimension
    dm = Unit(m, factor=0.1)

    base = (m, s)
    base_dim = (m.dim, s.dim)
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
    assert str(UnitSystem((m, s))) == "(m, s)"

    assert (repr(UnitSystem((m, s))) == "<UnitSystem: (%s, %s)>"
                                        % (m.abbrev_dim, s.abbrev_dim))


def test_call():
    A = Unit(current)
    Js = Unit(action)
    mksa = UnitSystem((m, kg, s, A), (Js,))

    assert mksa(Js) == mksa.print_unit_base(Js)
    assert mksa(Js.dim) == mksa._system(Js.dim)

    q = Quantity(10, Js)

    assert mksa(q) == "%g %s" % (q.factor, mksa(Js))


def test_get_unit():
    ms = UnitSystem((m, s), (c,))

    assert ms.get_unit("s") == s
    assert ms.get_unit(s) == s
    assert ms.get_unit(Unit(time)) == s

    assert ms["s"] == ms.get_unit("s")
    raises(KeyError, lambda: ms["g"])


def test_print_unit_base():
    A = Unit(current)
    Js = Unit(action)
    mksa = UnitSystem((m, kg, s, A), (Js,))

    assert mksa.print_unit_base(Js) == "0.001 m^2 kg s^-2"


def test_extend():
    ms = UnitSystem((m, s), (c,))
    Js = Unit(action)
    mks = ms.extend((kg,), (Js,))

    res = UnitSystem((m, s, kg), (c, Js))
    assert set(mks._base_units) == set(res._base_units)
    assert set(mks._units) == set(res._units)


def test_dim():
    dimsys = UnitSystem((m, kg, s), (c,))
    assert dimsys.dim == 3


def test_is_consistent():
    assert UnitSystem((m, s)).is_consistent is True
