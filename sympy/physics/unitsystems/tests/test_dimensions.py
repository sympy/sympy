# -*- coding: utf-8 -*-

from sympy import sympify
from sympy.physics.unitsystems.dimensions import Dimension
from sympy.utilities.pytest import raises


def test_definition():

    length = Dimension(name="length", symbol="L", length=1)

    assert length.get('length') == 1
    assert length.get('time') is None
    assert length.name == "length"
    assert length.symbol == "L"

    halflength = Dimension(name="length", length=0.5)
    assert halflength.get('length') == sympify("1/2")


def test_dict_properties():
    dic = {"length": sympify(1), "time": sympify(2)}
    d = Dimension(dic)

    assert d["length"] == 1
    assert set(d.items()) == set(dic.items())

    assert len(d) == 2

    assert d.get("length") == 1
    assert d.get("mass") is None

    assert ("length" in d) is True
    assert ("mass" in d) is False


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
    assert str(Dimension(length=1)) == "{'length': 1}"
    assert str(Dimension(length=1, symbol="L")) == "L"
    assert str(Dimension(length=1, name="length")) == "length"
    assert str(Dimension(length=1, symbol="L", name="length")) == 'L'


def test_properties():
    assert Dimension(length=1).is_dimensionless is False
    assert Dimension().is_dimensionless is True
    assert Dimension(length=0).is_dimensionless is True

    assert Dimension(length=1).has_integer_powers is True
    assert Dimension(length=-1).has_integer_powers is True
    assert Dimension(length=1.5).has_integer_powers is False


def test_add_sub():
    length = Dimension(length=1)

    assert length.add(length) == length
    assert length.sub(length) == length
    assert -length == length

    raises(TypeError, lambda: length.add(1))
    raises(TypeError, lambda: length.sub(1))
    raises(ValueError, lambda: length.add(Dimension(time=1)))
    raises(ValueError, lambda: length.sub(Dimension(time=1)))


def test_mul_div_exp():
    length = Dimension(length=1)
    time = Dimension(time=1)
    velocity = length.div(time)

    assert length.pow(2) == Dimension(length=2)
    assert length.mul(length) == length.pow(2)
    assert length.mul(time) == Dimension(length=1, time=1)
    assert velocity == Dimension(length=1, time=-1)
    assert velocity.pow(2) == Dimension(length=2, time=-2)

    raises(TypeError, lambda: length.pow("a"))
