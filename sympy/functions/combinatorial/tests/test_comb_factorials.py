from sympy import (Symbol, symbols, factorial, factorial2, binomial,
                   rf, ff, gamma, polygamma, EulerGamma, O, pi, nan,
                   oo, zoo, simplify, expand_func)
from sympy.functions.combinatorial.factorials import subfactorial
from sympy.utilities.pytest import XFAIL, raises


def test_rf_eval_apply():
    x, y = symbols('x,y')

    assert rf(nan, y) == nan

    assert rf(x, y) == rf(x, y)

    assert rf(oo, 0) == 1
    assert rf(-oo, 0) == 1

    assert rf(oo, 6) == oo
    assert rf(-oo, 7) == -oo

    assert rf(oo, -6) == oo
    assert rf(-oo, -7) == oo

    assert rf(x, 0) == 1
    assert rf(x, 1) == x
    assert rf(x, 2) == x*(x + 1)
    assert rf(x, 3) == x*(x + 1)*(x + 2)
    assert rf(x, 5) == x*(x + 1)*(x + 2)*(x + 3)*(x + 4)

    assert rf(x, -1) == 1/(x - 1)
    assert rf(x, -2) == 1/((x - 1)*(x - 2))
    assert rf(x, -3) == 1/((x - 1)*(x - 2)*(x - 3))

    assert rf(1, 100) == factorial(100)


def test_ff_eval_apply():
    x, y = symbols('x,y')

    assert ff(nan, y) == nan

    assert ff(x, y) == ff(x, y)

    assert ff(oo, 0) == 1
    assert ff(-oo, 0) == 1

    assert ff(oo, 6) == oo
    assert ff(-oo, 7) == -oo

    assert ff(oo, -6) == oo
    assert ff(-oo, -7) == oo

    assert ff(x, 0) == 1
    assert ff(x, 1) == x
    assert ff(x, 2) == x*(x - 1)
    assert ff(x, 3) == x*(x - 1)*(x - 2)
    assert ff(x, 5) == x*(x - 1)*(x - 2)*(x - 3)*(x - 4)

    assert ff(x, -1) == 1/(x + 1)
    assert ff(x, -2) == 1/((x + 1)*(x + 2))
    assert ff(x, -3) == 1/((x + 1)*(x + 2)*(x + 3))

    assert ff(100, 100) == factorial(100)


def test_factorial():
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True, positive=True)

    assert factorial(-2) == zoo
    assert factorial(0) == 1
    assert factorial(7) == 5040
    assert factorial(n).func == factorial
    assert factorial(2*n).func == factorial

    assert factorial(n).is_integer
    assert factorial(n).is_positive is None
    assert factorial(k).is_positive

    assert factorial(oo) == oo


def test_factorial_diff():
    n = Symbol('n', integer=True)

    assert factorial(n).diff(n) == \
        gamma(1 + n)*polygamma(0, 1 + n)
    assert factorial(n**2).diff(n) == \
        2*n*gamma(1 + n**2)*polygamma(0, 1 + n**2)


def test_factorial_series():
    n = Symbol('n', integer=True)

    assert factorial(n).series(n, 0, 3) == \
        1 - n*EulerGamma + n**2*(EulerGamma**2/2 + pi**2/12) + O(n**3)


def test_factorial_rewrite():
    n = Symbol('n', integer=True)

    assert factorial(n).rewrite(gamma) == gamma(n + 1)


def test_factorial2():
    n = Symbol('n', integer=True)

    assert factorial2(-1) == 1
    assert factorial2(0) == 1
    assert factorial2(7) == 105
    assert factorial2(8) == 384
    assert factorial2(n).func == factorial2


def test_binomial():
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True)
    u = Symbol('v', negative=True)
    v = Symbol('m', positive=True)

    assert binomial(0, 0) == 1
    assert binomial(1, 1) == 1
    assert binomial(10, 10) == 1
    assert binomial(1, 2) == 0
    assert binomial(1, -1) == 0
    assert binomial(-1, 1) == -1
    assert binomial(-10, 1) == -10
    assert binomial(-10, 7) == -11440
    assert binomial(n, -1) == 0
    assert binomial(n, 0) == 1
    assert expand_func(binomial(n, 1)) == n
    assert expand_func(binomial(n, 2)) == n*(n - 1)/2
    assert expand_func(binomial(n, n - 2)) == n*(n - 1)/2
    assert expand_func(binomial(n, n - 1)) == n
    assert binomial(n, 3).func == binomial
    assert binomial(n, 3).expand(func=True) ==  n**3/6 - n**2/2 + n/3
    assert expand_func(binomial(n, 3)) ==  n*(n - 2)*(n - 1)/6
    assert binomial(n, n) == 1
    assert binomial(n, n + 1) == 0
    assert binomial(n, u) == 0
    assert binomial(n, v).func == binomial
    assert binomial(n, k).func == binomial
    assert binomial(n, n + v) == 0

    assert expand_func(binomial(n, n-3)) == n*(n - 2)*(n - 1)/6


def test_binomial_diff():
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True)

    assert binomial(n, k).diff(n) == \
        (-polygamma(0, 1 + n - k) + polygamma(0, 1 + n))*binomial(n, k)
    assert binomial(n**2, k**3).diff(n) == \
        2*n*(-polygamma(
            0, 1 + n**2 - k**3) + polygamma(0, 1 + n**2))*binomial(n**2, k**3)

    assert binomial(n, k).diff(k) == \
        (-polygamma(0, 1 + k) + polygamma(0, 1 + n - k))*binomial(n, k)
    assert binomial(n**2, k**3).diff(k) == \
        3*k**2*(-polygamma(
            0, 1 + k**3) + polygamma(0, 1 + n**2 - k**3))*binomial(n**2, k**3)


def test_binomial_rewrite():
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True)

    assert binomial(n, k).rewrite(
        factorial) == factorial(n)/(factorial(k)*factorial(n - k))
    assert binomial(
        n, k).rewrite(gamma) == gamma(n + 1)/(gamma(k + 1)*gamma(n - k + 1))


@XFAIL
def test_factorial_simplify_fail():
    # simplify(factorial(x + 1).diff(x) - ((x + 1)*factorial(x)).diff(x))) == 0
    from sympy.abc import x
    assert simplify(x*polygamma(0, x + 1) - x*polygamma(0, x + 2) +
    polygamma(0, x + 1) - polygamma(0, x + 2) + 1) == 0


def test_subfactorial():
    assert all(subfactorial(i) == ans for i, ans in enumerate(
        [1, 0, 1, 2, 9, 44, 265, 1854, 14833, 133496]))
    raises(ValueError, lambda: subfactorial(0.1))
    raises(ValueError, lambda: subfactorial(-2))
