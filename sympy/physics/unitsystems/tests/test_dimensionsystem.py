# -*- coding: utf-8 -*-

from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem

length = Dimension(name="length", symbol="L", length=1)
time = Dimension(name="time", symbol="T", time=1)
velocity = Dimension(name="velocity", symbol="V", length=1, time=-1)


def test_definition():

    base = (length, time)
    ms = DimensionSystem(base, (velocity,), "MS", "MS system")

    assert ms._base_dims == base
    assert ms._dims == DimensionSystem._sort_dims(base + (velocity,))
    assert ms.name == "MS"
    assert ms.descr == "MS system"


def test_sort_dims():

    assert (DimensionSystem._sort_dims((length, velocity, time))
                                       == (length, time, velocity))
