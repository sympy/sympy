# -*- coding: utf-8 -*-

from sympy import sympify, Symbol, S, sqrt
from sympy.physics.units.dimensions import Dimension
from sympy.physics.units.dimensions import length, time
from sympy.utilities.pytest import raises


def test_definition():
    assert length.get_dimensional_dependencies() == {"length": 1}
    assert length.name == Symbol("length")
    assert length.symbol == Symbol("L")

    halflength = sqrt(length)
    assert halflength.get_dimensional_dependencies() == {'length': S.Half}


def test_error_definition():
    # tuple with more or less than two entries
    raises(TypeError, lambda: Dimension(("length", 1, 2)))
    raises(TypeError, lambda: Dimension(["length"]))

    # non-number power
    raises(TypeError, lambda: Dimension({"length": "a"}))

    # non-number with named argument
    raises(TypeError, lambda: Dimension({"length": (1, 2)}))


def test_str():
    assert str(Dimension("length")) == "Dimension(length)"
    assert str(Dimension("length", "L")) == "Dimension(length, L)"


def test_properties():
    assert length.is_dimensionless is False
    assert (length/length).is_dimensionless is True
    assert Dimension("undefined").is_dimensionless is True

    assert length.has_integer_powers is True
    assert (length**(-1)).has_integer_powers is True
    assert (length**1.5).has_integer_powers is False


def test_add_sub():
    assert length + length == length
    assert length - length == length
    assert -length == length

    raises(TypeError, lambda: length + 1)
    raises(TypeError, lambda: length - 1)
    raises(ValueError, lambda: length + time)
    raises(ValueError, lambda: length - time)


def test_mul_div_exp():
    velo = length / time

    assert (length * length) == length ** 2

    assert (length * length).get_dimensional_dependencies() == {"length": 2}
    assert (length ** 2).get_dimensional_dependencies() == {"length": 2}
    assert (length * time).get_dimensional_dependencies() == { "length": 1, "time": 1}
    assert velo.get_dimensional_dependencies() == { "length": 1, "time": -1}
    assert (velo ** 2).get_dimensional_dependencies() == {"length": 2, "time": -2}

    assert (length / length).get_dimensional_dependencies() == {}
    assert (velo / length * time).get_dimensional_dependencies() == {}
    assert (length ** -1).get_dimensional_dependencies() == {"length": -1}
    assert (velo ** -1.5).get_dimensional_dependencies() == {"length": -1.5, "time": 1.5}

    length_a = length**"a"
    assert length_a.get_dimensional_dependencies() == {"length": Symbol("a")}

    assert length != 1
    assert length / length != 1
