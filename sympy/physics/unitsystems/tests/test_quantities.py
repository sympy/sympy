# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Symbol

from sympy.physics.unitsystems.definitions import s, m, kg, c, SpeedOfLight, speed_of_light
from sympy.physics.unitsystems.dimensions import length, time
from sympy.physics.unitsystems.quantities import Quantity
from sympy.physics.unitsystems.prefixes import PREFIXES
from sympy.utilities.pytest import raises

k = PREFIXES["k"]


def test_definition():
    q = Quantity("s10", time, 10, "sabbr")

    assert q.factor == 10
    assert q.dimension == time
    assert q.abbrev == Symbol("sabbr")


def test_str_repr():
    assert str(kg) == "kiloGram"

def test_eq():
    # simple test
    assert 10*m == 10*m
    assert 10*m != 10*s


def test_convert_to():
    q = Quantity("q1", length, 5000)
    assert q.convert_to(m) == 5000*m

    assert speed_of_light.convert_to(m / s) == 299792458 * m / s
    assert (2*speed_of_light).convert_to(m / s) == 2 * 299792458 * m / s
    assert day.convert_to(s) == 86400*s
    raises(ValueError, lambda: q.convert_to(s))
    raises(ValueError, lambda: SpeedOfLight.convert_to(m))