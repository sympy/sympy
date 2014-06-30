# -*- coding: utf-8 -*-

from __future__ import division

from sympy.physics.unitsystems.quantities import Quantity
from sympy.physics.unitsystems.units import Unit
from sympy.physics.unitsystems.prefixes import PREFIXES
from sympy.physics.unitsystems.systems import mks
from sympy.utilities.pytest import raises

m, s, c = mks["m"], mks["s"], mks["c"]
k = PREFIXES["k"]


def test_definition():
    q = Quantity(10, s)

    assert q.factor == 10
    assert q.unit == s

    assert Quantity(10) == 10


def test_error_definition():
    raises(NotImplementedError, lambda: Quantity("10 m"))
    raises(TypeError, lambda: Quantity(10, m.dim))


def test_str_repr():
    q = Quantity(10, m)

    assert str(q) == "%g %s" % (10, str(m))
    assert repr(q) == "%g %s" % (10, repr(m))


def test_eq():
    # stupid test
    assert (Quantity(10, m) == Quantity(10, m)) is True
    assert (Quantity(10, m) == Quantity(10, s)) is False


def test_operations():
    q = Quantity(10, m)
    p = Quantity(5, s)

    assert -q == Quantity(-10, m)

    assert q.add(Quantity(20, m)) == Quantity(30, m)
    assert q.add(Quantity(20, Unit(m.dim, 10))) == Quantity(30, Unit(m.dim, 11))

    assert q.sub(Quantity(20, m)) == Quantity(-10, m)
    assert q.sub(Quantity(20, Unit(m.dim, 10))) == Quantity(-10, Unit(m.dim, -9))

    assert q.pow(2) == Quantity(10**2, m.pow(2))

    assert q.mul(p) == Quantity(5*10, m.mul(s))

    assert q.div(p) == Quantity(10/5, m.div(s))


def test_error_operations():
    q = Quantity(10, m)
    p = Quantity(5, s)

    raises(ValueError, lambda: q.add(p))
    raises(TypeError, lambda: q.add(1))

    raises(ValueError, lambda: q.sub(p))
    raises(TypeError, lambda: q.sub(1))


def test_convert_to():
    q = Quantity(5, Unit(m.dim, prefix=k))
    assert q.convert_to(m) == Quantity(5000, m)
