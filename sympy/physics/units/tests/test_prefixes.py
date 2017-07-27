# -*- coding: utf-8 -*-
from sympy import symbols, log
from sympy.physics.units import Quantity, Dimension
from sympy.physics.units.prefixes import PREFIXES, Prefix, prefix_unit, kilo, \
    kibi


def test_prefix_operations():
    m = PREFIXES['m']
    k = PREFIXES['k']
    M = PREFIXES['M']

    dodeca = Prefix('dodeca', 'dd', 1, base=12)

    assert m * k == 1
    assert k * k == M
    assert 1 / m == k
    assert k / m == M

    assert dodeca * dodeca == 144
    assert 1 / dodeca == 1 / 12
    assert k / dodeca == 1000 / 12
    assert dodeca / dodeca == 1

    m = Quantity("meter", 1, 6)
    assert dodeca / m == 12 / m


def test_prefix_unit():
    length = Dimension("length")
    m = Quantity("meter", length, 1, abbrev="m")

    pref = {"m": PREFIXES["m"], "c": PREFIXES["c"], "d": PREFIXES["d"]}

    res = [Quantity("millimeter", length, PREFIXES["m"], "mm"),
           Quantity("centimeter", length, PREFIXES["c"], "cm"),
           Quantity("decimeter", length, PREFIXES["d"], "dm")]

    prefs = prefix_unit(m, pref)
    assert set(prefs) == set(res)
    assert set(map(lambda x: x.abbrev, prefs)) == set(symbols("mm,cm,dm"))


def test_bases():
    assert kilo.base == 10
    assert kibi.base == 2


def test_repr():
    assert eval(repr(kilo)) == kilo
    assert eval(repr(kibi)) == kibi
