from sympy.functions.combinatorial.factorials import binomial, factorial
from sympy.functions.special.gamma_functions import gamma
from sympy.concrete.zeil import zeil
from sympy.abc import k, n, N, x, i
from sympy import degree


def test_zeil():
    fs = {binomial(n, k),
          binomial(n, k) * (-1)**(n - k),
          (-1)**(k - 1 + i) * binomial(n, k) * binomial(k - 1, i),
          binomial(n, k) / 2**n,
          factorial(n + k),
          factorial(k),
          binomial(n, k)**2,
          2**k * k / binomial(n, k),
          gamma(n + k) / gamma(2 * k)**2,
          (-1)**k * binomial(x - k + 1, k) * binomial(x - 2 * k, n - k)
          }

    for f in fs:
        res = zeil(f, k, n, N)
        assert res

        rec, cert = res
        assert not rec.is_zero

        d = degree(rec, N)
        left = sum(f.subs(n, n + i) * rec.coeff(N, i) for i in range(d+1))
        right = cert * f
        right = right.subs(k, k + 1) - right

        expr = (left - right).combsimp()

        assert expr.is_zero
