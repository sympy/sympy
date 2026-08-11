from __future__ import annotations
from sympy import srepr, sstr
from sympy.core.mul import Mul
from sympy.core.numbers import Rational
from sympy.core.singleton import S
from sympy.core.symbol import (Symbol, symbols)
from sympy.physics.units import Quantity, length, meter, W
from sympy.physics.units.prefixes import PREFIXES, Prefix, prefix_unit, kilo, \
    kibi
from sympy.physics.units.systems import SI

x = Symbol('x')


def test_prefix_operations():
    m = PREFIXES['m']
    k = PREFIXES['k']
    M = PREFIXES['M']

    dodeca = Prefix('dodeca', 'dd', 1, base=12)

    assert m * k is S.One
    assert m * W == W / 1000
    assert k * k == M
    assert 1 / m == k
    assert k / m == M

    assert dodeca * dodeca == 144
    assert 1 / dodeca == S.One / 12
    assert k / dodeca == S(1000) / 12
    assert dodeca / dodeca is S.One

    m = Quantity("fake_meter")
    SI.set_quantity_dimension(m, S.One)
    SI.set_quantity_scale_factor(m, S.One)

    assert dodeca * m == 12 * m
    assert dodeca / m == 12 / m

    expr1 = kilo * 3
    assert isinstance(expr1, Mul)
    assert expr1.args == (3, kilo)

    expr2 = kilo * x
    assert isinstance(expr2, Mul)
    assert expr2.args == (x, kilo)

    expr3 = kilo / 3
    assert isinstance(expr3, Mul)
    assert expr3.args == (Rational(1, 3), kilo)
    assert expr3.args == (S.One/3, kilo)

    expr4 = kilo / x
    assert isinstance(expr4, Mul)
    assert expr4.args == (1/x, kilo)


def test_prefix_unit():
    m = Quantity("fake_meter", abbrev="m")
    m.set_global_relative_scale_factor(1, meter)

    pref = {"m": PREFIXES["m"], "c": PREFIXES["c"], "d": PREFIXES["d"]}

    q1 = Quantity("millifake_meter", abbrev="mm")
    q2 = Quantity("centifake_meter", abbrev="cm")
    q3 = Quantity("decifake_meter", abbrev="dm")

    SI.set_quantity_dimension(q1, length)

    SI.set_quantity_scale_factor(q1, PREFIXES["m"])
    SI.set_quantity_scale_factor(q1, PREFIXES["c"])
    SI.set_quantity_scale_factor(q1, PREFIXES["d"])

    res = [q1, q2, q3]

    prefs = prefix_unit(m, pref)
    assert set(prefs) == set(res)
    assert {v.abbrev for v in prefs} == set(symbols("mm,cm,dm"))


def test_bases():
    assert kilo.base == 10
    assert kibi.base == 2


def test_prefix_printing():
    assert str(kilo) == "k"
    assert repr(kilo) == "kilo"
    assert sstr(kilo) == "kilo"
    assert sstr(kilo, abbrev=True) == "k"
    assert srepr(kilo) == "kilo"


def test_binary_prefix_abbreviations():
    # The binary prefixes previously all had the abbreviation "Y" (a copy-paste
    # from yotta); they must carry their own IEC abbreviations.  The abbreviation
    # is a name, so it must be a Symbol, not sympified (which would turn "Ei"
    # into the exponential-integral function and "E" into Euler's number).
    from sympy.physics.units.prefixes import (
        BIN_PREFIXES, kibi, mebi, gibi, tebi, pebi, exbi, exa)
    expected = {
        kibi: ("Ki", 2**10), mebi: ("Mi", 2**20), gibi: ("Gi", 2**30),
        tebi: ("Ti", 2**40), pebi: ("Pi", 2**50), exbi: ("Ei", 2**60),
    }
    for prefix, (abbrev, scale) in expected.items():
        assert str(prefix.abbrev) == abbrev
        assert isinstance(prefix.abbrev, Symbol)
        assert prefix.scale_factor == scale
    # BIN_PREFIXES is keyed by the correct abbreviations:
    assert {str(p.abbrev) for p in BIN_PREFIXES.values()} == set(BIN_PREFIXES)
    # a decimal prefix whose abbreviation collides with a constant stays a name:
    assert isinstance(exa.abbrev, Symbol)
    assert str(exa.abbrev) == "E"


def test_binary_prefix_operations():
    # Products/quotients of binary prefixes must simplify to binary prefixes,
    # just as decimal prefixes simplify to decimal prefixes.
    from sympy.physics.units.prefixes import (
        kibi, mebi, gibi, pebi, exbi, mega)
    assert kibi*kibi == mebi
    assert kibi*mebi == gibi
    assert mebi*kibi == gibi
    assert gibi/kibi == mebi
    assert exbi/pebi == kibi
    # decimal prefixes still simplify among themselves:
    assert kilo*kilo == mega
    # a mixed decimal/binary product has no prefix, so a number is returned:
    assert kilo*kibi == 1024000
    # 1/kibi has no (negative) binary prefix, so a plain fraction is returned:
    assert 1/kibi == Rational(1, 1024)
