from sympy import (IV, oo, Rational)

from sympy.utilities.pytest import raises
from sympy.utilities.pytest import raises, XFAIL


def test_IV():
    assert IV(1, 2).delta == Rational(1, 1)
    assert IV(1, 2).mid == Rational(3, 2)

    assert IV(1, 2) + 1 == IV(2, 3)
    assert 1 + IV(1, 2) == IV(2, 3)
    assert IV(1, 2) + IV(2, 3) == IV(3, 5)

    assert -IV(1, 2) == IV(-2, -1)

    assert IV(1, 2) - 1 == IV(0, 1)
    assert 1 - IV(1, 2) == IV(-1, 0)
    assert IV(2, 3) - IV(1, 2) == IV(0, 2)

    assert IV(1, 2)*2 == IV(2, 4)
    assert 2*IV(1, 2) == IV(2, 4)
    assert IV(1, 2)*IV(2, 3) == IV(2, 6)

    assert IV(1, 2)/2 == IV(Rational(1, 2), 1)
    assert 2/IV(2, 3) == IV(Rational(2, 3), 1)
    assert 1/IV(-1, 1) == IV(-oo, oo)
    assert IV(1, 2)/IV(-1, 1) == IV(-oo, oo)
    assert IV(1, 2)/IV(2, 3) == IV(Rational(1, 3), 1)

    assert IV(-1, 1)**2 == IV(0, 1)
    assert IV(1, 2)**2 == IV(1, 4)
