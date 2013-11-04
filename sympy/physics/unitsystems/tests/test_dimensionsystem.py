# -*- coding: utf-8 -*-

from sympy import Matrix, eye
from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem
from sympy.utilities.pytest import raises

length = Dimension(name="length", symbol="L", length=1)
mass = Dimension(name="mass", symbol="M", mass=1)
time = Dimension(name="time", symbol="T", time=1)
velocity = Dimension(name="velocity", symbol="V", length=1, time=-1)
acceleration = Dimension(name="acceleration", symbol="A", length=2, mass=1,
                         time=-2)


def test_definition():

    base = (length, time)
    ms = DimensionSystem(base, (velocity,), "MS", "MS system")

    assert ms._base_dims == base
    assert ms._dims == DimensionSystem._sort_dims(base + (velocity,))
    assert ms.name == "MS"
    assert ms.descr == "MS system"

    assert ms._can_transf_matrix is None
    assert ms._list_can_dims is None


def test_error_definition():
    raises(ValueError,
           lambda: DimensionSystem((Dimension(time=1, name="time", symbol="T"),
                                    Dimension(length=1, symbol="L"),
                                    Dimension(mass=1, name="mass"),
                                    Dimension(current=1))))


def test_str_repr():
    assert str(DimensionSystem((length, time), name="MS")) == "MS"
    dimsys = DimensionSystem((Dimension(length=1, symbol="L"),
                              Dimension(time=1, symbol="T")))
    assert str(dimsys) == "(L, T)"
    dimsys = DimensionSystem((Dimension(length=1, name="length", symbol="L"),
                              Dimension(time=1, name="time")))
    assert str(dimsys) == "(L, time)"

    assert (repr(DimensionSystem((length, time), name="MS"))
            == "<DimensionSystem: ({length: 1}, {time: 1})>")


def test_get_dim():
    ms = DimensionSystem((length, time), (velocity,))

    assert ms.get_dim("L") == length
    assert ms.get_dim("length") == length
    assert ms.get_dim(length) == length
    assert ms.get_dim(Dimension(length=1)) == length

    assert ms["L"] == ms.get_dim("L")
    raises(KeyError, lambda: ms["M"])


def test_extend():
    ms = DimensionSystem((length, time), (velocity,))

    mks = ms.extend((mass,), (acceleration,))

    assert mks._base_dims == DimensionSystem((length, time, mass),
                                        (velocity, acceleration))._base_dims
    assert mks._dims == DimensionSystem((length, time, mass),
                                        (velocity, acceleration))._dims


def test_sort_dims():

    assert (DimensionSystem._sort_dims((length, velocity, time))
                                       == (length, time, velocity))


def test_list_dims():
    dimsys = DimensionSystem((length, time, mass))

    assert dimsys.list_can_dims == ("length", "mass", "time")


def test_dim_can_vector():
    dimsys = DimensionSystem((length, mass, time), (velocity, acceleration))

    assert dimsys.dim_can_vector(length) == Matrix([1, 0, 0])
    assert dimsys.dim_can_vector(velocity) == Matrix([1, 0, -1])

    dimsys = DimensionSystem((length, velocity, acceleration), (mass, time))

    assert dimsys.dim_can_vector(length) == Matrix([1, 0, 0])
    assert dimsys.dim_can_vector(velocity) == Matrix([1, 0, -1])


def test_dim_vector():
    dimsys = DimensionSystem((length, mass, time), (velocity, acceleration))

    assert dimsys.dim_vector(length) == Matrix([1, 0, 0])
    assert dimsys.dim_vector(velocity) == Matrix([1, 0, -1])

    dimsys = DimensionSystem((length, velocity, acceleration), (mass, time))

    assert dimsys.dim_vector(length) == Matrix([0, 1, 0])
    assert dimsys.dim_vector(velocity) == Matrix([0, 0, 1])
    assert dimsys.dim_vector(time) == Matrix([0, 1, -1])


def test_can_transf_matrix():
    dimsys = DimensionSystem((length, mass, time))

    assert dimsys._can_transf_matrix is None
    assert dimsys.can_transf_matrix == eye(3)
    assert dimsys._can_transf_matrix == eye(3)

    dimsys = DimensionSystem((length, velocity, acceleration))
    assert dimsys.can_transf_matrix == Matrix(((0, 1, 0), (1, 0, 1),
                                               (0, -2, -1)))
