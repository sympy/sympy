# -*- coding: utf-8 -*-

from __future__ import division

from sympy.physics.unitsystems.units import Unit
from sympy.physics.unitsystems.systems.mks import length
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
