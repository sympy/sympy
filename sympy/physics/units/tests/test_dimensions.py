# -*- coding: utf-8 -*-
import warnings

from sympy.utilities.exceptions import SymPyDeprecationWarning

from sympy import S, Symbol, sqrt
from sympy.physics.units.dimensions import Dimension, length, time, dimsys_default
from sympy.utilities.pytest import raises


def test_Dimension_definition():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)
        assert length.get_dimensional_dependencies() == {"length": 1}
    assert dimsys_default.get_dimensional_dependencies(length) == {"length": 1}
    assert length.name == Symbol("length")
    assert length.symbol == Symbol("L")

    halflength = sqrt(length)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)
        assert halflength.get_dimensional_dependencies() == {"length": S.Half}
    assert dimsys_default.get_dimensional_dependencies(halflength) == {"length": S.Half}


def test_Dimension_error_definition():
    # tuple with more or less than two entries
    raises(TypeError, lambda: Dimension(("length", 1, 2)))
    raises(TypeError, lambda: Dimension(["length"]))

    # non-number power
    raises(TypeError, lambda: Dimension({"length": "a"}))

    # non-number with named argument
    raises(TypeError, lambda: Dimension({"length": (1, 2)}))

    # symbol should by Symbol or str
    raises(AssertionError, lambda: Dimension("length", symbol=1))


def test_Dimension_error_regisration():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)

        # tuple with more or less than two entries
        raises(IndexError, lambda: length._register_as_base_dim())

        one = Dimension(1)
        raises(TypeError, lambda: one._register_as_base_dim())


def test_str():
    assert str(Dimension("length")) == "Dimension(length)"
    assert str(Dimension("length", "L")) == "Dimension(length, L)"


def test_Dimension_properties():
    assert dimsys_default.is_dimensionless(length) is False
    assert dimsys_default.is_dimensionless(length/length) is True
    assert dimsys_default.is_dimensionless(Dimension("undefined")) is True

    assert length.has_integer_powers(dimsys_default) is True
    assert (length**(-1)).has_integer_powers(dimsys_default) is True
    assert (length**1.5).has_integer_powers(dimsys_default) is False


def test_Dimension_add_sub():
    assert length + length == length
    assert length - length == length
    assert -length == length

    raises(TypeError, lambda: length + 1)
    raises(TypeError, lambda: length - 1)
    raises(ValueError, lambda: length + time)
    raises(ValueError, lambda: length - time)


def test_Dimension_mul_div_exp():
    velo = length / time

    assert (length * length) == length ** 2

    assert dimsys_default.get_dimensional_dependencies(length * length) == {"length": 2}
    assert dimsys_default.get_dimensional_dependencies(length ** 2) == {"length": 2}
    assert dimsys_default.get_dimensional_dependencies(length * time) == { "length": 1, "time": 1}
    assert dimsys_default.get_dimensional_dependencies(velo) == { "length": 1, "time": -1}
    assert dimsys_default.get_dimensional_dependencies(velo ** 2) == {"length": 2, "time": -2}

    assert dimsys_default.get_dimensional_dependencies(length / length) == {}
    assert dimsys_default.get_dimensional_dependencies(velo / length * time) == {}
    assert dimsys_default.get_dimensional_dependencies(length ** -1) == {"length": -1}
    assert dimsys_default.get_dimensional_dependencies(velo ** -1.5) == {"length": -1.5, "time": 1.5}

    length_a = length**"a"
    assert dimsys_default.get_dimensional_dependencies(length_a) == {"length": Symbol("a")}

    assert length != 1
    assert length / length != 1

    length_0 = length ** 0
    assert dimsys_default.get_dimensional_dependencies(length_0) == {}
