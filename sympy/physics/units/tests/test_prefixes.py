# -*- coding: utf-8 -*-
from sympy import symbols, log, Mul, Symbol
from sympy.physics.units import Quantity, Dimension, meter
from sympy.physics.units.prefixes import PREFIXES, Prefix, prefix_unit, kilo, \
    kibi

x = Symbol('x')

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
    assert dodeca * m == 12 * m
    assert dodeca / m == 12 / m

    expr1 = kilo * 3
    assert isinstance(expr1, Mul)
    assert (expr1).args == (3, kilo)

    expr2 = kilo * x
    assert isinstance(expr2, Mul)
    assert (expr2).args == (x, kilo)

    expr3 = kilo / 3
    assert isinstance(expr3, Mul)
    assert (expr3).args == (1/3, kilo)

    expr4 = kilo / x
    assert isinstance(expr4, Mul)
    assert (expr4).args == (1/x, kilo)


def test_prefix_unit():
    length = Dimension("length")
    m = Quantity("meter", length, 1, abbrev="m")

    pref = {"m": PREFIXES["m"], "c": PREFIXES["c"], "d": PREFIXES["d"]}

    res = [Quantity("millimeter", length, PREFIXES["m"], abbrev="mm"),
           Quantity("centimeter", length, PREFIXES["c"], abbrev="cm"),
           Quantity("decimeter", length, PREFIXES["d"], abbrev="dm")]

    prefs = prefix_unit(m, pref)
    assert set(prefs) == set(res)
    assert set(map(lambda x: x.abbrev, prefs)) == set(symbols("mm,cm,dm"))


def test_bases():
    assert kilo.base == 10
    assert kibi.base == 2


def test_repr():
    assert eval(repr(kilo)) == kilo
    assert eval(repr(kibi)) == kibi
