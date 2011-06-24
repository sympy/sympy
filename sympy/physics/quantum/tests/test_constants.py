from sympy import Float

from sympy.physics.quantum.constants import hbar

def test_hbar():
    assert hbar.is_commutative == True
    assert hbar.is_real == True
    assert hbar.is_positive == True
    assert hbar.is_negative == False
    assert hbar.is_irrational == True

    assert hbar.evalf() == Float(1.05457162e-34)
