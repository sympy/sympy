from __future__ import annotations

from sympy import symbols, factorial, harmonic, S, cancel
from sympy.concrete.karr_field import KarrField
from sympy.concrete.karr_solver import karr_sum

def test_karr_field():
    x, y = symbols('x y')
    field = KarrField()
    field.add_base(x)
    field.add_pi(y, S(2)) # y represents 2**x, sigma(y) = 2*y

    elem1 = field.element(x + y)
    assert elem1.shift().expr == x + 1 + 2*y

    elem2 = field.element(x * y)
    assert elem2.shift().expr == cancel((x + 1) * (2*y))


def test_karr_sum():
    k, n = symbols('k n', integer=True)

    # 1. Arithmetic progression: sum_{k=1}^n k
    res1 = karr_sum(k, (k, 1, n))
    assert res1 == cancel(n*(n + 1)/2)

    # 2. Geometric progression: sum_{k=1}^n 2^k
    res2 = karr_sum(2**k, (k, 1, n))
    assert res2 == cancel(2**(n + 1) - 2)

    # 3. Factorial summation: sum_{k=1}^n k * k!
    res3 = karr_sum(k * factorial(k), (k, 1, n))
    assert res3 == cancel(factorial(n + 1) - 1)

    # 4. Rational function summation: sum_{k=1}^n 1/(k*(k+1))
    res4 = karr_sum(S(1)/(k*(k+1)), (k, 1, n))
    assert res4 == cancel(n/(n + 1))

    # 5. Harmonic number summation: sum_{k=1}^n H_k
    res5 = karr_sum(harmonic(k), (k, 1, n))
    assert cancel((res5 - ((n + 1)*harmonic(n) - n)).subs(harmonic(n+1), harmonic(n) + 1/(n+1))) == 0
