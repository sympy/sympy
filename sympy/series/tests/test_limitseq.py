from sympy import symbols, oo, Sum, harmonic
from sympy import difference_delta as dd
from sympy.utilities.pytest import raises

n, m, k = symbols('n m k', integer=True)


def test_difference_delta():
    e = n*(n + 1)
    e2 = e * k

    assert dd(e) == 2*n + 2
    assert dd(e2, n, 2) == k*(4*n + 6)

    raises(ValueError, lambda: dd(e2))
    raises(ValueError, lambda: dd(e2, n, oo))


def test_difference_delta__Sum():
    e = Sum(1/k, (k, 1, n))
    assert dd(e, n) == 1/(n + 1)
    assert dd(e, n, 5) == Sum(1/k, (k, n + 1, n + 5))

    e = Sum(1/k, (k, 1, 3*n))
    assert dd(e, n) == Sum(1/k, (k, 3*n + 1, 3*n + 3))

    e = n * Sum(1/k, (k, 1, n))
    assert dd(e, n) == 1 + Sum(1/k, (k, 1, n))

    e = Sum(1/k, (k, 1, n), (m, 1, n))
    assert dd(e, n) == harmonic(n)


def test_difference_delta__Add():
    e = n + n*(n + 1)
    assert dd(e, n) == 2*n + 3
    assert dd(e, n, 2) == 4*n + 8

    e = n + Sum(1/k, (k, 1, n))
    assert dd(e, n) == 1 + 1/(n + 1)
    assert dd(e, n, 5) == 5 + Sum(1/k, (k, n + 1, n + 5))
