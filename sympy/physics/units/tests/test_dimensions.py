from sympy.physics.units.systems.si import dimsys_SI
from sympy.utilities.pytest import warns_deprecated_sympy

from sympy import S, Symbol, sqrt
from sympy.physics.units.dimensions import Dimension
from sympy.physics.units.definitions.dimension_definitions import (
    length, time
)
from sympy.physics.units import foot
from sympy.utilities.pytest import raises


def test_Dimension_definition():
    assert dimsys_SI.get_dimensional_dependencies(length) == {"length": 1}
    assert length.name == Symbol("length")
    assert length.symbol == Symbol("L")

    halflength = sqrt(length)
    assert dimsys_SI.get_dimensional_dependencies(halflength) == {"length": S.Half}


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


def test_str():
    assert str(Dimension("length")) == "Dimension(length)"
    assert str(Dimension("length", "L")) == "Dimension(length, L)"


def test_Dimension_properties():
    assert dimsys_SI.is_dimensionless(length) is False
    assert dimsys_SI.is_dimensionless(length/length) is True
    assert dimsys_SI.is_dimensionless(Dimension("undefined")) is False

    assert length.has_integer_powers(dimsys_SI) is True
    assert (length**(-1)).has_integer_powers(dimsys_SI) is True
    assert (length**1.5).has_integer_powers(dimsys_SI) is False


def test_Dimension_add_sub():
    assert length + length == length
    assert length - length == length
    assert -length == length

    raises(TypeError, lambda: length + foot)
    raises(TypeError, lambda: foot + length)
    raises(TypeError, lambda: length - foot)
    raises(TypeError, lambda: foot - length)

    # issue 14547 - only raise error for dimensional args; allow
    # others to pass
    x = Symbol('x')
    e = length + x
    assert e == x + length and e.is_Add and set(e.args) == {length, x}
    e = length + 1
    assert e == 1 + length == 1 - length and e.is_Add and set(e.args) == {length, 1}


def test_Dimension_mul_div_exp():
    assert 2*length == length*2 == length/2 == length
    assert 2/length == 1/length
    x = Symbol('x')
    m = x*length
    assert m == length*x and m.is_Mul and set(m.args) == {x, length}
    d = x/length
    assert d == x*length**-1 and d.is_Mul and set(d.args) == {x, 1/length}
    d = length/x
    assert d == length*x**-1 and d.is_Mul and set(d.args) == {1/x, length}

    velo = length / time

    assert (length * length) == length ** 2

    assert dimsys_SI.get_dimensional_dependencies(length * length) == {"length": 2}
    assert dimsys_SI.get_dimensional_dependencies(length ** 2) == {"length": 2}
    assert dimsys_SI.get_dimensional_dependencies(length * time) == { "length": 1, "time": 1}
    assert dimsys_SI.get_dimensional_dependencies(velo) == { "length": 1, "time": -1}
    assert dimsys_SI.get_dimensional_dependencies(velo ** 2) == {"length": 2, "time": -2}

    assert dimsys_SI.get_dimensional_dependencies(length / length) == {}
    assert dimsys_SI.get_dimensional_dependencies(velo / length * time) == {}
    assert dimsys_SI.get_dimensional_dependencies(length ** -1) == {"length": -1}
    assert dimsys_SI.get_dimensional_dependencies(velo ** -1.5) == {"length": -1.5, "time": 1.5}

    length_a = length**"a"
    assert dimsys_SI.get_dimensional_dependencies(length_a) == {"length": Symbol("a")}

    assert length != 1
    assert length / length != 1

    length_0 = length ** 0
    assert dimsys_SI.get_dimensional_dependencies(length_0) == {}
