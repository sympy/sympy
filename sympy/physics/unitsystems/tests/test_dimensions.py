# -*- coding: utf-8 -*-

from sympy.physics.unitsystems.dimensions import Dimension
from sympy.utilities.pytest import raises


def test_definition():

    length = Dimension(name="length", symbol="L", length=1)

    assert length.get('length') == 1
    assert length.get('time') is None
    assert length.name == "length"
    assert length.symbol == "L"


def test_error_definition():
    # tuple with more or less than two entries
    raises(ValueError, lambda: Dimension(("length", 1, 2)))
    raises(ValueError, lambda: Dimension(["length"]))

    # non-number power
    raises(TypeError, lambda: Dimension(length="a"))

    # non-dict/list/tuple as positional arg
    raises(TypeError, lambda: Dimension("length"))

    # non-number with named argument
    raises(TypeError, lambda: Dimension(length=(1, 2)))


def test_str():
    assert str(Dimension(length=1)) == "{length: 1}"
    assert str(Dimension(length=1, symbol="L")) == "L"
    assert str(Dimension(length=1, name="length")) == "length"
    assert str(Dimension(length=1, symbol="L", name="length")) == 'L'


def test_properties():
    assert Dimension(length=1).is_dimensionless is False
    assert Dimension().is_dimensionless is True
    assert Dimension(length=0).is_dimensionless is True
