# -*- coding: utf-8 -*-

from sympy import Matrix, eye
from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem, length, time, velocity, mass, current, \
    action
from sympy.utilities.pytest import raises


def test_definition():

    base = (length, time)
    ms = DimensionSystem(base, (velocity,), "MS", "MS system")

    assert ms._base_dims == DimensionSystem.sort_dims(base)
    assert set(ms._dims) == set(base + (velocity,))
    assert ms.name == "MS"
    assert ms.descr == "MS system"


def test_error_definition():
    raises(ValueError,
           lambda: DimensionSystem((Dimension("time", {"time": 1}, symbol="T"),
                                    Dimension({"length": 1}, symbol="L"),
                                    Dimension("mass", {"mass": 1}),
                                    Dimension(current=1))))

    raises(ValueError, lambda: DimensionSystem((length, time, velocity)))


def test_str_repr():
    assert str(DimensionSystem((length, time), name="MS")) == "MS"
    dimsys = DimensionSystem((Dimension({"length": 1}, symbol="L"),
                              Dimension({"time": 1}, symbol="T")))
    assert str(dimsys) == '(Dimension({length: 1}, L), Dimension({time: 1}, T))'
    dimsys = DimensionSystem((Dimension("length", {"length": 1}, symbol="L"),
                              Dimension("time", {"time": 1})))
    assert str(dimsys) == '(Dimension(length, {length: 1}, L), Dimension(time, {time: 1}))'

    assert (repr(DimensionSystem((length, time), name="MS"))
            == '<DimensionSystem: (Dimension(length, {length: 1}, L), Dimension(time, {time: 1}, T))>')


def test_call():
    mksa = DimensionSystem((length, time, mass, current), (action,))
    assert mksa(action) == mksa.print_dim_base(action)


def test_get_dim():
    ms = DimensionSystem((length, time), (velocity,))

    assert ms.get_dim("L") == length
    assert ms.get_dim("length") == length
    assert ms.get_dim(length) == length
    assert ms.get_dim(Dimension({"length": 1})) == length

    assert ms["L"] == ms.get_dim("L")
    raises(KeyError, lambda: ms["M"])


def test_extend():
    ms = DimensionSystem((length, time), (velocity,))

    mks = ms.extend((mass,), (action,))

    res = DimensionSystem((length, time, mass), (velocity, action))
    assert mks._base_dims == res._base_dims
    assert set(mks._dims) == set(res._dims)


def test_sort_dims():

    assert (DimensionSystem.sort_dims((length, velocity, time))
                                      == (length, time, velocity))


def test_list_dims():
    dimsys = DimensionSystem((length, time, mass))

    assert dimsys.list_can_dims == ("length", "mass", "time")


def test_dim_can_vector():
    dimsys = DimensionSystem((length, mass, time), (velocity, action))

    assert dimsys.dim_can_vector(length) == [1, 0, 0]
    assert dimsys.dim_can_vector(velocity) == [1, 0, -1]

    dimsys = DimensionSystem((length, velocity, action), (mass, time))

    assert dimsys.dim_can_vector(length) == [1, 0, 0]
    assert dimsys.dim_can_vector(velocity) == [1, 0, -1]


def test_dim_vector():
    dimsys = DimensionSystem((length, mass, time), (velocity, action))

    assert dimsys.dim_vector(length) == Matrix([1, 0, 0])
    assert dimsys.dim_vector(velocity) == Matrix([1, 0, -1])

    dimsys = DimensionSystem((length, velocity, action), (mass, time))

    assert dimsys.dim_vector(length) == Matrix([0, 1, 0])
    assert dimsys.dim_vector(velocity) == Matrix([0, 0, 1])
    assert dimsys.dim_vector(time) == Matrix([0, 1, -1])


def test_inv_can_transf_matrix():
    dimsys = DimensionSystem((length, mass, time))

    assert dimsys.inv_can_transf_matrix == eye(3)

    dimsys = DimensionSystem((length, velocity, action))
    assert dimsys.inv_can_transf_matrix == Matrix(((2, 1, 1), (1, 0, 0),
                                                   (-2, 0, -1)))


def test_can_transf_matrix():
    dimsys = DimensionSystem((length, mass, time))

    assert dimsys.can_transf_matrix == eye(3)

    dimsys = DimensionSystem((length, velocity, action))
    assert dimsys.can_transf_matrix == Matrix(((0, 1, 0), (1, 0, 1),
                                               (0, -2, -1)))


def test_is_consistent():
    assert DimensionSystem((length, time)).is_consistent is True
    #assert DimensionSystem((length, time, velocity)).is_consistent is False


def test_print_dim_base():
    mksa = DimensionSystem((length, time, mass, current), (action,))
    assert mksa.print_dim_base(action) == "L^2 M T^-2"


def test_dim():
    dimsys = DimensionSystem((length, mass, time), (velocity, action))
    assert dimsys.dim == 3
