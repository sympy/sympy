from sympy import I
from sympy.physics.paulialgebra import Pauli

def test_Pauli():
    sigma1=Pauli(1)
    sigma2=Pauli(2)
    sigma3=Pauli(3)

    assert sigma1 == sigma1
    assert sigma1 != sigma2

    assert sigma1*sigma2 == I*sigma3
    assert sigma3*sigma1 == I*sigma2
    assert sigma2*sigma3 == I*sigma1

    assert sigma1*sigma1 == 1
    assert sigma2*sigma2 == 1
    assert sigma3*sigma3 == 1

    assert sigma1**0 == 1
    assert sigma1**1 == sigma1
    assert sigma1**2 == 1
    assert sigma1**3 == sigma1
    assert sigma1**4 == 1

    assert sigma3**2 == 1


    assert sigma1*2*sigma1 == 2
    #assert sigma1*sigma3*sigma1 == -sigma3     XXX should work
