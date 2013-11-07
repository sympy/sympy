# -*- coding: utf-8 -*-

from sympy.physics.unitsystems.prefixes import PREFIXES


def test_prefix_operations():
    m = PREFIXES['m']
    k = PREFIXES['k']
    M = PREFIXES['M']

    assert m * k == 1
    assert k * k == M
    assert 1 / m == k
    assert k / m == M
