from sympy import symbols, differenceDelta, oo, Sum, harmonic
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


def test_differenceDelta__Sum():
    e = Sum(1/k, (k, 1, n))
    assert e.differenceDelta(n) == 1/(n + 1)
    assert e.differenceDelta(n, 5) == Sum(1/k, (k, n + 1, n + 5))

    e = Sum(1/k, (k, 1, 3*n))
    assert e.differenceDelta(n) == Sum(1/k, (k, 3*n + 1, 3*n + 3))

    e = n * Sum(1/k, (k, 1, n))
    assert e.differenceDelta(n) == 1 + Sum(1/k, (k, 1, n))

    e = Sum(1/k, (k, 1, n), (m, 1, n))
    assert e.differenceDelta(n) == harmonic(n)
