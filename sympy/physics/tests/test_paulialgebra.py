from sympy import I, symbols
from sympy.physics.paulialgebra import Pauli
from sympy.utilities.pytest import XFAIL

sigma1 = Pauli(1)
sigma2 = Pauli(2)
sigma3 = Pauli(3)

tau1 = symbols("tau1", commutative = False)


def test_Pauli():

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


def test_evaluate_pauli_product():
    from sympy.physics.paulialgebra import evaluate_pauli_product

    assert evaluate_pauli_product(I*sigma2*sigma3) == -sigma1

    # Check issue 6471
    assert evaluate_pauli_product(-I*4*sigma1*sigma2) == 4*sigma3

    # Acting on non-Mul objects should return the input
    assert evaluate_pauli_product(1) == 1
    # After consecutive multiplication, one or no Pauli should remain
    assert evaluate_pauli_product(I*sigma1*sigma2*sigma1*sigma2) == -I
    # Must respect non-commuting properties of other symbols
    assert evaluate_pauli_product(I*sigma1*sigma2*tau1*sigma1*sigma3) == I*sigma3*tau1*sigma2


@XFAIL
def test_Pauli_should_work():
    assert sigma1*sigma3*sigma1 == -sigma3
