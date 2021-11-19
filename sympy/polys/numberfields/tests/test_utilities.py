"""Tests on utilities. """

from sympy.core.numbers import (I, Rational)
from sympy.core.singleton import S
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.testing.pytest import raises
from sympy.polys.numberfields.utilities import isolate
from sympy.printing.lambdarepr import IntervalPrinter
from sympy.abc import x


def test_IntervalPrinter():
    ip = IntervalPrinter()
    assert ip.doprint(x**Rational(1, 3)) == "x**(mpi('1/3'))"
    assert ip.doprint(sqrt(x)) == "x**(mpi('1/2'))"


def test_isolate():
    assert isolate(1) == (1, 1)
    assert isolate(S.Half) == (S.Half, S.Half)

    assert isolate(sqrt(2)) == (1, 2)
    assert isolate(-sqrt(2)) == (-2, -1)

    assert isolate(sqrt(2), eps=Rational(1, 100)) == (Rational(24, 17), Rational(17, 12))
    assert isolate(-sqrt(2), eps=Rational(1, 100)) == (Rational(-17, 12), Rational(-24, 17))

    raises(NotImplementedError, lambda: isolate(I))
