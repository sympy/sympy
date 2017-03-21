# -*- coding: utf-8 -*-

from __future__ import division

from sympy import Symbol, Number, Add, S

from sympy.physics.unitsystems.units import Unit
from sympy.physics.unitsystems.dimensions import length, time
from sympy.physics.unitsystems.prefixes import PREFIXES
from sympy.utilities.pytest import raises

k = PREFIXES['k']


def test_definition():
    u = Unit("u", length, factor=10, abbrev="dam")

    assert u.dimension == length
    assert u.factor == 10
    assert u.abbrev == Symbol("dam")
    assert u.prefix is None

    km = Unit("km", length, 1, prefix=k)
    assert km.prefix == k
    assert km.prefix.factor == 1000
    assert km.factor == 1000
    assert km.func(*km.args) == km
    assert km.func(*km.args).args == km.args

    v = Unit("u", length, factor=5, prefix=k)
    assert v.dimension == length
    assert v.factor == 5 * 1000


def test_abbrev():
    u = Unit("u", length, 1)
    assert u.name == Symbol("u")
    assert u.abbrev == Symbol("")

    u = Unit("u", length, 2, abbrev="m")
    assert u.name == Symbol("u")
    assert u.abbrev == Symbol("m")
    assert u.factor == 2
    assert isinstance(u.factor, Number)

    u = Unit("u", length, 3, abbrev="m", prefix=k)
    assert u.abbrev == Symbol("km")
    assert u.factor == 3000


def test_print():
    u = Unit("unitname", length, factor=10, abbrev="m")
    assert repr(u) == "unitname"
    assert str(u) == "unitname"


def test_eq():
    u = Unit("u", length, factor=10, abbrev="dam")

    v = Unit("v1", length, factor=10)
    assert u == v

    v = Unit("v2", time, factor=10, abbrev="ds")
    assert u != v

    v = Unit("v3", length, factor=1, abbrev="dm")
    assert u != v


def test_add_sub():
    u = Unit("u", length, factor=10)
    v = Unit("v", length, factor=5)
    w = Unit("w", time, factor=2)

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
    u = Unit("u", length, factor=10)
    v = Unit("v", length, factor=5)
    w = Unit("w", time, factor=2)

    # TODO: no way of checking unit consistency:
    #raises(ValueError, lambda: check_unit_consistency(u + w))
    #raises(ValueError, lambda: check_unit_consistency(u - w))
    #raises(TypeError, lambda: check_unit_consistency(u + 1))
    #raises(TypeError, lambda: check_unit_consistency(u - 1))


def test_mul_div():
    u = Unit("u", length, factor=10)

    assert 1 / u == u**(-1)
    assert u / 1 == u

    v1 = u / Unit("t", time, factor=2)
    v2 = Unit("v", length / time, factor=5)

    # Pow only supports structural equality:
    assert v1 != v2
    assert v1 == v2.convert_to(v1)

    # TODO: decide whether to allow such expression in the future
    # (requires somehow manipulating the core).
    #assert u / Unit(length, factor=2) == 5

    assert u * 1 == u

    ut1 = u * Unit("t", time, factor=2)
    ut2 = Unit("ut", length*time, factor=20)

    # Mul only supports structural equality:
    assert ut1 != ut2
    assert ut1 == ut2.convert_to(ut1)

    # Mul only supports structural equality:
    assert u * Unit("lp1", length**-1, factor=2) != 20

    assert u**0 == 1
    assert u**1 == u
    # TODO: Pow only support structural equality:
    assert u ** 2 != Unit("u2", length ** 2, factor=100)
    assert u ** -1 != Unit("u3", length ** -1, factor=0.1)

    # TODO: conversion to dimensional power:
    # assert u ** 2 == Unit("u2", length ** 2, factor=100).convert_to(u)
    # assert u ** -1 == Unit("u3", length ** -1, factor=0.1).convert_to(u)
