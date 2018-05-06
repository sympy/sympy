from sympy.concrete.zeilberger import zb_recur
from sympy import factorial, binomial, simplify, symbols

n, k, x, i, j, r = symbols('n, k, x, i, j, r')

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

    assert simplify(zb_recur(F_1, 2, n, k)[1]) == simplify(R_1)
    assert simplify(zb_recur(F_2, 1, n, k)[1]) == simplify(R_2)
    assert simplify(zb_recur(F_3, 1, n, k)[1]) == simplify(R_3)
    assert zb_recur(F_4, 1, n, k) == "try higher order"
    assert simplify(zb_recur(F_4, 2, n, k)[1]) == simplify(R_4)
    assert simplify(zb_recur(F_5, 1, n, k)[1]) == simplify(R_5)
    assert simplify(zb_recur(F_6, 1, n, k)[1]) == simplify(R_6)
