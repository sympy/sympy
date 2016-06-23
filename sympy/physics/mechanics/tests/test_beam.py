from sympy import Symbol
from sympy.physics.mechanics.beam import Beam, PointLoad, DistributedLoad
from sympy.physics.mechanics import Point
from sympy.printing import sstr


def test_Beam():
    E = Symbol('E')
    I = Symbol('I')
    b = Beam(1, E, I)
    assert b.length == 1
    assert b.E == E
    assert b.I == I


def test_PointLoad():
    P1 = Point('P1')
    P2 = Point('P2')
    Load_1 = PointLoad(location = P1, value = -4)
    assert Load_1.location == P1
    assert Load_1.value == -4
    assert Load_1.moment is False

    Load_2 = PointLoad(location = P1, value = 5, moment=True)
    assert Load_2.location == P1
    assert Load_2.value == 5
    assert Load_2.moment is True


def test_DistributedLoad():
    P1 = Point('P1')
    P2 = Point('P2')
    x = Symbol('x')
    Load_1 = DistributedLoad(start = P1, end = P2, value = -4*x**2)
    assert Load_1.start == P1
    assert Load_1.end == P2
    assert Load_1.value == -4*x**2
