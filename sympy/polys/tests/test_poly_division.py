import pytest
from sympy import ZZ, QQ, symbols
from sympy.polys.polyerrors import ExactQuotientFailed

x, y = symbols("x y")

D = ZZ
R = D[x]
xp = R.gens[0]

def test_poly_division_examples():
    assert xp / xp == 1
    with pytest.raises(ExactQuotientFailed):
        xp / (2*xp)
    with pytest.raises(ExactQuotientFailed):
        xp / 2
    with pytest.raises(ExactQuotientFailed):
        xp / (2*D.one)
    with pytest.raises(ExactQuotientFailed):
        xp / (2*R.one)
    with pytest.raises(ExactQuotientFailed):
        2 / xp

def test_poly_division_extended():
    assert 0 / xp == 0
    with pytest.raises(ZeroDivisionError):
        xp / 0
    assert (-xp) / (-xp) == 1
    with pytest.raises(ExactQuotientFailed):
        (-xp) / 2

    QQx = QQ[x]
    qxp = QQx.gens[0]
    assert (qxp / 2) == qxp / 2

    ZZx = ZZ[x]
    zxp = ZZx.gens[0]
    with pytest.raises(ExactQuotientFailed):
        zxp / 3

    ZZxy = ZZ[x, y]
    f = ZZxy(x + y)
    g = ZZxy(x + y)
    assert f / g == 1
    with pytest.raises(ExactQuotientFailed):
        f / (2*g)
    with pytest.raises(ExactQuotientFailed):
        f / 2

    f2 = ZZxy(x - y)
    g2 = ZZxy(x - y)
    assert f2 / g2 == 1
    with pytest.raises(ExactQuotientFailed):
        f2 / (3*g2)
    assert 0 / f2 == 0