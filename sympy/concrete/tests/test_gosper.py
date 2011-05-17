"""Tests for Gosper's algorithm for hypergeometric summation. """

from sympy.concrete.gosper import gosper_normal, gosper_term, gosper_sum
from sympy import S, Poly, symbols, factorial, binomial, gamma

n, m, k, i, j = symbols('n,m,k,i,j', integer=True)

def test_gosper_normal():
    assert gosper_normal(4*n + 5, 2*(4*n + 1)*(2*n + 3), n) == \
        (Poly(S(1)/4, n), Poly(n + S(3)/2), Poly(n + S(1)/4))

def test_gosper_term():
    assert gosper_term((4*k + 1)*factorial(k)/factorial(2*k + 1), k) == (-k - S(1)/2)/(k + S(1)/4)

def test_gosper_sum():
    assert gosper_sum(1, (k, 0, n)) == 1 + n
    assert gosper_sum(k, (k, 0, n)) == n*(1 + n)/2
    assert gosper_sum(k**2, (k, 0, n)) == n*(1 + n)*(1 + 2*n)/6
    assert gosper_sum(k**3, (k, 0, n)) == n**2*(1 + n)**2/4

    assert gosper_sum(2**k, (k, 0, n)) == 2*2**n - 1

    assert gosper_sum(factorial(k), (k, 0, n)) is None
    assert gosper_sum(binomial(n, k), (k, 0, n)) is None

    assert gosper_sum(factorial(k)/k**2, (k, 0, n)) is None
    assert gosper_sum((k - 3)*factorial(k), (k, 0, n)) is None

    assert gosper_sum(k*factorial(k), k) == factorial(k)
    assert gosper_sum(k*factorial(k), (k, 0, n)) == n*factorial(n) + factorial(n) - 1

    assert gosper_sum((-1)**k*binomial(n, k), (k, 0, n)) == 0
    assert gosper_sum((-1)**k*binomial(n, k), (k, 0, m)) == -(-1)**m*(m - n)*binomial(n, m)/n

    assert gosper_sum((4*k + 1)*factorial(k)/factorial(2*k + 1), (k, 0, n)) == \
        (2*factorial(2*n + 1) - factorial(n))/factorial(2*n + 1)

def test_gosper_sum_indefinite():
    assert gosper_sum(k, k) == k*(k - 1)/2
    assert gosper_sum(k**2, k) == k*(k - 1)*(2*k - 1)/6

    assert gosper_sum(1/(k*(k + 1)), k) == -1/k
    assert gosper_sum(-(27*k**4 + 158*k**3 + 430*k**2 + 678*k + 445)*gamma(2*k + 4)/(3*(3*k + 7)*gamma(3*k + 6)), k) == \
        (3*k + 5)*(k**2 + 2*k + 5)*gamma(2*k + 4)/gamma(3*k + 6)

def test_gosper_sum_parametric():
    assert gosper_sum(binomial(S(1)/2, m - j + 1)*binomial(S(1)/2, m + j), (j, 1, n)) == \
        n*(1 + m - n)*(-1 + 2*m + 2*n)*binomial(S(1)/2, 1 + m - n)*binomial(S(1)/2, m + n)/(m*(1 + 2*m))

def test_gosper_sum_iterated():
    f1 = binomial(2*k, k)/4**k
    f2 = (1 + 2*n)*binomial(2*n, n)/4**n
    f3 = (1 + 2*n)*(3 + 2*n)*binomial(2*n, n)/(3*4**n)
    f4 = (1 + 2*n)*(3 + 2*n)*(5 + 2*n)*binomial(2*n, n)/(15*4**n)
    f5 = (1 + 2*n)*(3 + 2*n)*(5 + 2*n)*(7 + 2*n)*binomial(2*n, n)/(105*4**n)

    assert gosper_sum(f1, (k, 0, n)) == f2
    assert gosper_sum(f2, (n, 0, n)) == f3
    assert gosper_sum(f3, (n, 0, n)) == f4
    assert gosper_sum(f4, (n, 0, n)) == f5

