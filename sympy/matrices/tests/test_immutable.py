from sympy import *
from sympy.utilities.pytest import raises, XFAIL

IM = ImmutableMatrix([[1,2,3], [4,5,6], [7,8,9]])
ieye = ImmutableMatrix(eye(3))

def test_immutable_creation():

    assert IM.shape == (3,3)
    assert IM[1,2] == 6
    assert IM[2,2] == 9

def test_immutability():
    raises(TypeError, "IM[2,2] = 5")

def test_slicing():
    IM[1,:] == ImmutableMatrix([[4,5,6]])
    IM[:2, :2] == ImmutableMatrix([[1,2],[4,5]])
