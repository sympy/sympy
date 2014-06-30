# -*- coding: utf-8 -*-

from __future__ import division

from sympy.physics.unitsystems.units import Unit
from sympy.physics.unitsystems.systems.mks import length, time
from sympy.physics.unitsystems.prefixes import PREFIXES
from sympy.utilities.pytest import raises

k = PREFIXES['k']


def test_definition():
    u = Unit(length, factor=10, abbrev="dm")

    assert u.dim == length
    assert u._factor == 10
    assert u._abbrev == "dm"
    assert u.prefix is None

    km = Unit(length, prefix=k)
    assert km.prefix == k

    v = Unit(u, factor=5)
    assert v.dim == length
    assert v._factor == 5 * 10


def test_error_definition():
    raises(TypeError, lambda: Unit("m"))


def test_factor():
    u = Unit(length, factor=10, abbrev="dm")
    assert u.factor == 10

    u = Unit(length, factor=5, prefix=k)
    assert u.factor == 5000


def test_abbrev():
    u = Unit(length)
    assert u.abbrev == ""

    u = Unit(length, abbrev="m")
    assert u.abbrev == "m"

    u = Unit(length, abbrev="m", prefix=k)
    assert u.abbrev == "km"


def test_abbrev_dim():
    u = Unit(length, factor=10)
    assert u.abbrev_dim == "(10 L)"


def test_str():
    u = Unit(length, factor=10)
    assert str(u) == u.abbrev_dim

    u = Unit(length, factor=10, abbrev="m")
    assert str(u) == "m"


def test_repr():
    u = Unit(length, factor=10, abbrev="m")
    assert repr(u) == u.abbrev_dim


def test_eq():
    u = Unit(length, factor=10, abbrev="dm")

    v = Unit(length, factor=10)
    assert (u == v) is True

    v = Unit(time, factor=10, abbrev="ds")
    assert (u == v) is False

    v = Unit(length, factor=1, abbrev="dm")
    assert (u == v) is False


def test_add_sub():
    u = Unit(length, factor=10)
    v = Unit(length, factor=5)
    w = Unit(time, factor=2)

    assert u.add(v) == Unit(length, factor=15)
    assert u.sub(v) == Unit(length, factor=5)

    raises(ValueError, lambda: u.add(w))
    raises(ValueError, lambda: u.sub(w))
    raises(TypeError, lambda: u.add(1))
    raises(TypeError, lambda: u.sub(1))


def test_pow():
    u = Unit(length, factor=10)

    assert u.pow(0) == 1
    assert u.pow(1) == u
    assert u.pow(2) == Unit(length.pow(2), factor=100)
    assert u.pow(-1) == Unit(length.pow(-1), factor=0.1)


def test_mul():
    u = Unit(length, factor=10)

    assert u.mul(1) == u

    assert u.mul(Unit(time, factor=2)) == Unit(length.mul(time), factor=20)
    assert u.mul(Unit(length.pow(-1), factor=2)) == 20


def test_div():
    u = Unit(length, factor=10)

    assert u.rdiv(1) == u.pow(-1)
    assert u.div(1) == u

    assert u.div(Unit(time, factor=2)) == Unit(length.div(time), factor=5)
    assert u.div(Unit(length, factor=2)) == 5


def test_is_compatible():
    u = Unit(length, factor=10)

    assert u.is_compatible(Unit(length)) is True
    assert u.is_compatible(Unit(time)) is False
    assert u.is_compatible(2) is False


def test_as_quantity():
    from sympy.physics.unitsystems.quantities import Quantity

    u = Unit(length, factor=10)
    q = Quantity(10, Unit(length))

    assert u.as_quantity == q
