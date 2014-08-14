# -*- coding: utf-8 -*-

from sympy.physics.unitsystems.prefixes import PREFIXES, prefix_unit


def test_prefix_operations():
    m = PREFIXES['m']
    k = PREFIXES['k']
    M = PREFIXES['M']

    assert m * k == 1
    assert k * k == M
    assert 1 / m == k
    assert k / m == M


def test_prefix_unit():
    from sympy.physics.unitsystems import Unit, Dimension

    length = Dimension(length=1)
    m = Unit(length, abbrev="m")

    pref = {"m": PREFIXES["m"], "c": PREFIXES["c"], "d": PREFIXES["d"]}

    res = [Unit(length, abbrev="m", prefix=PREFIXES["m"]),
           Unit(length, abbrev="m", prefix=PREFIXES["c"]),
           Unit(length, abbrev="m", prefix=PREFIXES["d"])]

    assert set(prefix_unit(m, pref)) == set(res)
