from sympy import Symbol
from sympy.physics.mechanics.beam import Beam, PointLoad, DistributedLoad
from sympy.physics.mechanics import Point
from sympy.printing import sstr


def test_Beam():
    E = Symbol('E')
    E_1 = Symbol('E_1')
    I = Symbol('I')
    I_1 = Symbol('I_1')
    b = Beam(1, E, I)
    assert b._length == 1
    assert b._E == E
    assert b._I == I

    # Test the length setter
    b.length = 4
    assert b.length == 4

    # Test the E setter
    b.E = E_1
    assert b.E == E_1

    # Test the I setter
    b.I = I_1
    assert b.I is I_1


def test_PointLoad():
    P1 = Point('P1')
    P2 = Point('P2')
    Load_1 = PointLoad(location = P1, value = -4)
    assert Load_1.location == P1
    assert Load_1.value == -4
    assert Load_1.moment is False

    # Test the location setter
    Load_1.location = P2
    assert Load_1.location == P2

    # Test the value setter
    Load_1.value = 4
    assert Load_1.value == 4

    # Test the moment setter
    Load_1.moment = True
    assert Load_1.moment is True

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

    # Test the start setter
    Load_1.start = P2
    assert Load_1.start == P2

    # Test the end setter
    Load_1.end = P1
    assert Load_1.end == P1

    # Test the value setter
    Load_1.value = 4*x**2 + x
    assert Load_1.value == 4*x**2 + x
