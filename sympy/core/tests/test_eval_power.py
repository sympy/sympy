from sympy.core import *

def test_rational():
    a = Rational(1, 5)

    assert a**Rational(1, 2) == a**Rational(1, 2)
    assert 2 * a**Rational(1, 2) == 2 * a**Rational(1, 2)

    assert a**Rational(3, 2) == a * a**Rational(1, 2)
    assert 2 * a**Rational(3, 2) == 2*a * a**Rational(1, 2)

    assert a**Rational(17, 3) == a**5 * a**Rational(2, 3)
    assert 2 * a**Rational(17, 3) == 2*a**5 * a**Rational(2, 3)