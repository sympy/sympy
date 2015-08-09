from sympy import symbols, differenceDelta, oo
from sympy.utilities.pytest import raises

n, m, k = symbols('n m k', integer=True)


def test_differenceDelta():
    e = n*(n + 1)
    e2 = e * k

    assert e.differenceDelta() == 2*n + 2
    assert e2.differenceDelta(n, 2) == k*(4*n + 6)
    assert e.differenceDelta(n, 2) == differenceDelta(e, n, 2)

    raises(ValueError, lambda: e2.differenceDelta())
    raises(ValueError, lambda: e2.differenceDelta(n, oo))
