# -*- coding: utf-8 -*-

from __future__ import division

from sympy.physics.unitsystems.quantities import Quantity
from sympy.physics.unitsystems.systems import mks
from sympy.utilities.pytest import raises

m, s, c = mks["m"], mks["s"], mks["c"]


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
