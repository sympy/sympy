# -*- coding: utf-8 -*-
from sympy import symbols
from sympy.physics.units.prefixes import PREFIXES, prefix_unit


def test_prefix_operations():
    m = PREFIXES['m']
    k = PREFIXES['k']
    M = PREFIXES['M']

    assert m * k == 1
    assert k * k == M
    assert 1 / m == k
    assert k / m == M


def test_prefix_unit():
    from sympy.physics.units import Quantity, Dimension

    length = Dimension("length")
    m = Quantity("meter", length, 1, abbrev="m")

    pref = {"m": PREFIXES["m"], "c": PREFIXES["c"], "d": PREFIXES["d"]}

    res = [Quantity("millimeter", length, PREFIXES["m"], "mm"),
           Quantity("centimeter", length, PREFIXES["c"], "cm"),
           Quantity("decimeter", length, PREFIXES["d"], "dm")]

    prefs = prefix_unit(m, pref)
    assert set(prefs) == set(res)
    assert set(map(lambda x: x.abbrev, prefs)) == set(symbols("mm,cm,dm"))
