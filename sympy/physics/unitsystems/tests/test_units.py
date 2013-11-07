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


def test_factor():
    u = Unit(length, factor=10, abbrev="dm")
    assert u.factor == 10

    u = Unit(length, factor=5, prefix=k)
    assert u.factor == 5000
