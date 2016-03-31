# -*- coding: utf-8 -*-

from sympy import Matrix, eye
from sympy.physics.unitsystems.dimensions import Dimension, DimensionSystem
from sympy.utilities.pytest import raises

length = Dimension(name="length", symbol="L", length=1)
mass = Dimension(name="mass", symbol="M", mass=1)
time = Dimension(name="time", symbol="T", time=1)
velocity = Dimension(name="velocity", symbol="V", length=1, time=-1)
action = Dimension(name="action", symbol="A", length=2, mass=1, time=-2)


def test_definition():

    base = (length, time)
    ms = DimensionSystem(base, (velocity,), "MS", "MS system")

    assert ms._base_dims == DimensionSystem.sort_dims(base)
    assert set(ms._dims) == set(base + (velocity,))
    assert ms.name == "MS"
    assert ms.descr == "MS system"

    assert ms._can_transf_matrix is None

    # be sure that there is no duplicates in ms._dims
    dims = (Dimension(length=1), Dimension(length=1, symbol="l"))
    ms = DimensionSystem(base, dims)

    assert set(ms._dims) == set(base)


def test_error_definition():
    raises(ValueError,
           lambda: DimensionSystem((Dimension(time=1, name="time", symbol="T"),
                                    Dimension(length=1, symbol="L"),
                                    Dimension(mass=1, name="mass"),
                                    Dimension(current=1))))

    raises(ValueError, lambda: DimensionSystem((length, time, velocity)))


def test_str_repr():
    assert str(DimensionSystem((length, time), name="MS")) == "MS"
    dimsys = DimensionSystem((Dimension(length=1, symbol="L"),
                              Dimension(time=1, symbol="T")))
    assert str(dimsys) == "(L, T)"
    dimsys = DimensionSystem((Dimension(length=1, name="length", symbol="L"),
                              Dimension(time=1, name="time")))
    assert str(dimsys) == "(L, time)"

    assert (repr(DimensionSystem((length, time), name="MS"))
            == "<DimensionSystem: ({'length': 1}, {'time': 1})>")


def test_call():
    current = Dimension(name="current", symbol="I", current=1)
    mksa = DimensionSystem((length, time, mass, current), (action,))
    assert mksa(action) == mksa.print_dim_base(action)


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

    assert dimsys.dim_can_vector(length) == Matrix([1, 0, 0])
    assert dimsys.dim_can_vector(velocity) == Matrix([1, 0, -1])

    dimsys = DimensionSystem((length, velocity, action), (mass, time))

    assert dimsys.dim_can_vector(length) == Matrix([1, 0, 0])
    assert dimsys.dim_can_vector(velocity) == Matrix([1, 0, -1])


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

    assert dimsys._can_transf_matrix is None
    assert dimsys.can_transf_matrix == eye(3)
    assert dimsys._can_transf_matrix == eye(3)

    dimsys = DimensionSystem((length, velocity, action))
    assert dimsys.can_transf_matrix == Matrix(((0, 1, 0), (1, 0, 1),
                                               (0, -2, -1)))


def test_is_consistent():
    assert DimensionSystem((length, time)).is_consistent is True
    #assert DimensionSystem((length, time, velocity)).is_consistent is False


def test_print_dim_base():
    current = Dimension(name="current", symbol="I", current=1)
    mksa = DimensionSystem((length, time, mass, current), (action,))
    assert mksa.print_dim_base(action) == "L^2 M T^-2"


def test_dim():
    dimsys = DimensionSystem((length, mass, time), (velocity, action))
    assert dimsys.dim == 3
