from sympy.concrete.zeilberger import zb_recur, zb_sum
from sympy import factorial, binomial, simplify, symbols, sqrt, pi, gamma
from sympy import combsimp, summation, Rational

k, x, i, j, r, y = symbols('k, x, i, j, r, y')
n = symbols('n', integer = True, positive = True)

def test_zb_recur():

    F_1 = (-1)**k * binomial(x - k + 1, k) * binomial(x - 2*k, n - k)
    F_2 = binomial(n, k) / 2**n
    F_3 = (-1)**k * binomial(n, k) * x * binomial(x + n, n) / (k + x)
    F_4 = binomial(n - k - 1, k)
    F_5 = (factorial(n - i) * factorial(n - j) * factorial(i - 1) *
        factorial(j - 1)) / (factorial(n - 1) * factorial(k - 1) *
        factorial(n - i - j + k) * factorial(i - k) * factorial(j - k))
    F_6 = binomial(n, k) * binomial(x, k + r) / binomial(n + x, n + r)

    R_1 = k*(-2*k + x + 1)*(-k + x + 2)/(2*(k - n - 2)*(k - n - 1))
    R_2 = k/(2*(k - n - 1))
    R_3 = k*(k + x)/(k - n - 1)
    R_4 = k*(k - n)/((2*k - n)*(2*k - n - 1))
    R_5 = -(k - 1)/n
    R_6 = k*(k + r)/((k - n - 1)*(n + x + 1))

    assert simplify(zb_recur(F_1, n, k, J = 2)[1]) == simplify(R_1)
    assert simplify(zb_recur(F_2, n, k)[1]) == simplify(R_2)
    assert simplify(zb_recur(F_3, n, k)[1]) == simplify(R_3)
    assert simplify(zb_recur(F_4, n, k, J = 2)[1]) == simplify(R_4)
    assert simplify(zb_recur(F_5, n, k)[1]) == simplify(R_5)
    assert simplify(zb_recur(F_6, n, k)[1]) == simplify(R_6)

def test_zb_sum():

    F_1 = (-1)**k * binomial(n + k, 2 * k) * binomial(2 * k, k) / (k + 1)
    F_2 = binomial(n, k)
    F_3 = (-1)**k * binomial(n, k) * x / (k + x)
    F_4 = binomial(n, k)**2
    F_5 = k * binomial(2 * n + 1, 2 * k + 1)
    F_6 = (-1)**k * binomial(n, k) / binomial(x + k, k)

    R_1 = 0
    R_2 = 2**n
    R_3 = gamma(n + 1)*gamma(x + 1)/gamma(n + x + 1)
    R_4 = 4**n*gamma(n + Rational(1/2))/(sqrt(pi)*gamma(n + 1))
    R_5 = 4*2**(2*n - 3)*gamma(n + Rational(1/2))/gamma(n - Rational(1/2))
    R_6 = x/(n + x)

    assert combsimp(zb_sum(F_1, (k, 0, n), J = 2)[0]) == R_1
    assert combsimp(zb_sum(F_2, (k, 0, n))[0]) == R_2
    assert combsimp(zb_sum(F_3, (k, 0, n))[0]) == R_3
    assert combsimp(zb_sum(F_4, (k, 0, n))[0]) == R_4
    assert combsimp(zb_sum(F_5, (k, 0, n))[0]) == R_5
    assert combsimp(zb_sum(F_6, (k, 0, n))[0]) == R_6

def test_new_summation():

    F_1 = binomial(n + k, k) / 2**k
    F_2 = binomial(x, k) * binomial(y, n - k)

    R_1 = 2**n
    R_2 = (-1)**n*gamma(n - x - y)/(gamma(n + 1)*gamma(-x - y))

    assert summation(F_1, (k, 0, n)) == R_1
    assert summation(F_2, (k, 0, n)) == R_2
