"""Tests for quotient rings."""

from sympy import QQ, ZZ
from sympy.abc import x, y

from sympy.polys.polyerrors import NotReversible

from sympy.utilities.pytest import raises

def test_QuotientRingElement():
    R = QQ[x]/[x**10]
    X = R.convert(x)

    assert X*(X+1) == R.convert(x**2 + x)
    assert X*x == R.convert(x**2)
    assert x*X == R.convert(x**2)
    assert X+x == R.convert(2*x)
    assert x+X == 2*X
    assert X**2 == R.convert(x**2)
    assert 1/(1 - X) == R.convert(sum(x**i for i in range(10)))
    assert X**10 == R.zero
    assert X != x

    raises(NotReversible, lambda: 1/X)

def test_QuotientRing():
    I = QQ[x].ideal(x**2 + 1)
    R = QQ[x]/I

    assert R == QQ[x]/[x**2 + 1]
    assert R == QQ[x]/QQ[x].ideal(x**2 + 1)
    assert R != QQ[x]

    assert R.convert(1)/x == -x + I
    assert -1 + I == x**2 + I
    assert R.convert(1, ZZ) == 1 + I
    assert R.convert(R.convert(x), R) == R.convert(x)

    X = R.convert(x)
    Y = QQ[x].convert(x)
    assert -1 + I == X**2 + I
    assert -1 + I == Y**2 + I
    assert R.to_sympy(X) == x

    raises(ValueError, lambda: QQ[x]/QQ[x, y].ideal(x))

    R = QQ.poly_ring(x, order="ilex")
    I = R.ideal(x)
    assert R.convert(1) + I == (R/I).convert(1)
