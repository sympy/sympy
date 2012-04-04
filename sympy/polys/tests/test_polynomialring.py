"""Tests for the PolynomialRing classes. """

from sympy.polys.domains import QQ, ZZ, PolynomialRing
from sympy.polys.polyerrors import ExactQuotientFailed, CoercionFailed

from sympy.abc import x, y

from sympy.utilities.pytest import raises

def test_build_order():
    R = QQ.poly_ring(x, y, order=(("lex", x), ("ilex", y)))
    assert R.order((1, 5)) == ((1,), (-5,))

def test_globalring():
    Qxy = QQ.frac_field(x, y)
    R   = QQ[x, y]
    X = R.convert(x)
    Y = R.convert(y)

    assert x in R
    assert 1/x not in R
    assert 1/(1 + x) not in R
    assert Y in R
    assert X.ring == R
    assert X * (Y**2 + 1) == R.convert(x * (y**2 + 1))
    assert X * y == X * Y == R.convert(x * y) == x * Y
    assert X + y == X + Y == R.convert(x + y) == x + Y
    assert X + 1 == R.convert(x + 1)
    raises(ExactQuotientFailed, 'X/Y')
    raises(ExactQuotientFailed, 'x/Y')
    raises(ExactQuotientFailed, 'X/y')
    assert X**2 / X == X

    assert R.from_GlobalPolynomialRing(ZZ[x, y].convert(x), ZZ[x, y]) == X
    assert R.from_FractionField(Qxy.convert(x), Qxy) == X
    assert R.from_FractionField(Qxy.convert(x)/y, Qxy) is None

def test_localring():
    Qxy = QQ.frac_field(x, y)
    R   = QQ.poly_ring(x, y, order="ilex")
    X = R.convert(x)
    Y = R.convert(y)

    assert x in R
    assert 1/x not in R
    assert 1/(1 + x) in R
    assert Y in R
    assert X.ring == R
    assert X*(Y**2+1)/(1 + X) == R.convert(x*(y**2 + 1)/(1 + x))
    assert X*y == X*Y
    raises(ExactQuotientFailed, 'X/Y')
    raises(ExactQuotientFailed, 'x/Y')
    raises(ExactQuotientFailed, 'X/y')
    assert X + y == X + Y == R.convert(x + y) == x + Y
    assert X + 1 == R.convert(x + 1)
    assert X**2 / X == X

    assert R.from_GlobalPolynomialRing(ZZ[x, y].convert(x), ZZ[x, y]) == X
    assert R.from_FractionField(Qxy.convert(x), Qxy) == X
    raises(CoercionFailed, 'R.from_FractionField(Qxy.convert(x)/y, Qxy)')

def test_conversion():
    L = QQ.poly_ring(x, y, order="ilex")
    G = QQ[x, y]

    assert L.convert(x) == L.convert(G.convert(x), G)
    assert G.convert(x) == G.convert(L.convert(x), L)
    raises(CoercionFailed, 'G.convert(L.convert(1/(1+x)), L)')
