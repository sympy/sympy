# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Symbol, Add, Number, S

from sympy.physics.unitsystems.definitions import s, m, kg, speed_of_light, day
from sympy.physics.unitsystems.dimensions import length, time
from sympy.physics.unitsystems.quantities import Quantity
from sympy.physics.unitsystems.prefixes import PREFIXES, kilo

k = PREFIXES["k"]


def test_str_repr():
    assert str(kg) == "kilogram"

def test_eq():
    # simple test
    assert 10*m == 10*m
    assert 10*m != 10*s


def test_convert_to():
    q = Quantity("q1", length, 5000)
    assert q.convert_to(m) == 5000*m

    assert speed_of_light.convert_to(m / s) == 299792458 * m / s
    # TODO: eventually support this kind of conversion:
    # assert (2*speed_of_light).convert_to(m / s) == 2 * 299792458 * m / s
    assert day.convert_to(s) == 86400*s

    # Wrong dimension to convert:
    assert q.convert_to(s) == q
    assert speed_of_light.convert_to(m) == speed_of_light


def test_Quantity_definition():
    q = Quantity("s10", time, 10, "sabbr")

    assert q.scale_factor == 10
    assert q.dimension == time
    assert q.abbrev == Symbol("sabbr")

    u = Quantity("u", length, 10, abbrev="dam")

    assert u.dimension == length
    assert u.scale_factor == 10
    assert u.abbrev == Symbol("dam")

    km = Quantity("km", length, kilo)
    assert km.scale_factor == 1000
    assert km.func(*km.args) == km
    assert km.func(*km.args).args == km.args

    v = Quantity("u", length, 5*kilo)
    assert v.dimension == length
    assert v.scale_factor == 5 * 1000


def test_abbrev():
    u = Quantity("u", length, 1)
    assert u.name == Symbol("u")
    assert u.abbrev == Symbol("u")

    u = Quantity("u", length, 2, "om")
    assert u.name == Symbol("u")
    assert u.abbrev == Symbol("om")
    assert u.scale_factor == 2
    assert isinstance(u.scale_factor, Number)

    u = Quantity("u", length, 3*kilo, "ikm")
    assert u.abbrev == Symbol("ikm")
    assert u.scale_factor == 3000


def test_print():
    u = Quantity("unitname", length, 10, "dam")
    assert repr(u) == "unitname"
    assert str(u) == "unitname"


def test_Quantity_eq():
    u = Quantity("u", length, 10, "dam")

    v = Quantity("v1", length, 10)
    assert u != v

    v = Quantity("v2", time, 10, "ds")
    assert u != v

    v = Quantity("v3", length, 1, "dm")
    assert u != v


def test_add_sub():
    u = Quantity("u", length, 10)
    v = Quantity("v", length, 5)
    w = Quantity("w", time, 2)

    assert isinstance(u + v, Add)
    assert (u + v.convert_to(u)) == (1 + S.Half)*u
    # TODO: eventually add this:
    # assert (u + v).convert_to(u) == (1 + S.Half)*u
    assert isinstance(u - v, Add)
    assert (u - v.convert_to(u)) == S.Half*u
    # TODO: eventually add this:
    # assert (u - v).convert_to(u) == S.Half*u


def test_check_unit_consistency():
    return  # TODO remove
    u = Quantity("u", length, 10)
    v = Quantity("v", length, 5)
    w = Quantity("w", time, 2)

    # TODO: no way of checking unit consistency:
    #raises(ValueError, lambda: check_unit_consistency(u + w))
    #raises(ValueError, lambda: check_unit_consistency(u - w))
    #raises(TypeError, lambda: check_unit_consistency(u + 1))
    #raises(TypeError, lambda: check_unit_consistency(u - 1))


def test_mul_div():
    u = Quantity("u", length, 10)

    assert 1 / u == u**(-1)
    assert u / 1 == u

    v1 = u / Quantity("t", time, 2)
    v2 = Quantity("v", length / time, 5)

    # Pow only supports structural equality:
    assert v1 != v2
    assert v1 == v2.convert_to(v1)

    # TODO: decide whether to allow such expression in the future
    # (requires somehow manipulating the core).
    #assert u / Quantity(length, 2) == 5

    assert u * 1 == u

    ut1 = u * Quantity("t", time, 2)
    ut2 = Quantity("ut", length*time, 20)

    # Mul only supports structural equality:
    assert ut1 != ut2
    assert ut1 == ut2.convert_to(ut1)

    # Mul only supports structural equality:
    assert u * Quantity("lp1", length**-1, 2) != 20

    assert u**0 == 1
    assert u**1 == u
    # TODO: Pow only support structural equality:
    assert u ** 2 != Quantity("u2", length ** 2, 100)
    assert u ** -1 != Quantity("u3", length ** -1, 0.1)

    # TODO: conversion to dimensional power:
    # assert u ** 2 == Quantity("u2", length ** 2, 100).convert_to(u)
    # assert u ** -1 == Quantity("u3", length ** -1, 0.1).convert_to(u)
