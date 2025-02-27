from sympy.series import approximants
from sympy.series.approximants import pade_approximants, pade_approximant
from sympy.core.symbol import symbols
from sympy.functions.combinatorial.factorials import binomial
from sympy.functions.combinatorial.numbers import (fibonacci, lucas)
from sympy import exp, sin, log


def test_approximants():
    x, t = symbols("x,t")
    g = [lucas(k) for k in range(16)]
    assert list(approximants(g)) == (
        [2, -4/(x - 2), (5*x - 2)/(3*x - 1), (x - 2)/(x**2 + x - 1)] )
    g = [lucas(k)+fibonacci(k+2) for k in range(16)]
    assert list(approximants(g)) == (
        [3, -3/(x - 1), (3*x - 3)/(2*x - 1), -3/(x**2 + x - 1)] )
    g = [lucas(k)**2 for k in range(16)]
    assert list(approximants(g)) == (
        [4, -16/(x - 4), (35*x - 4)/(9*x - 1), (37*x - 28)/(13*x**2 + 11*x - 7),
        (50*x**2 + 63*x - 52)/(37*x**2 + 19*x - 13),
        (-x**2 - 7*x + 4)/(x**3 - 2*x**2 - 2*x + 1)] )
    p = [sum(binomial(k,i)*x**i for i in range(k+1)) for k in range(16)]
    y = approximants(p, t, simplify=True)
    assert next(y) == 1
    assert next(y) == -1/(t*(x + 1) - 1)


def test_pade_approximants():
    x = symbols("x")

    # test exp
    all_order_5_exp_pade_approximations = [
    x**5/120 + x**4/24 + x**3/6 + x**2/2 + x + 1,
    (-x**4 - 8*x**3 - 36*x**2 - 96*x - 120)/(24*(x - 5)),
    (x**3 + 9*x**2 + 36*x + 60)/(3*(x**2 - 8*x + 20)),
    3*(-x**2 - 8*x - 20)/(x**3 - 9*x**2 + 36*x - 60),
    24*(x + 5)/(x**4 - 8*x**3 + 36*x**2 - 96*x + 120),
    -120/(x**5 - 5*x**4 + 20*x**3 - 60*x**2 + 120*x - 120)
    ]

    exp_pade = pade_approximants(exp(x), x, 5)

    for p_true, (num, denom) in zip(all_order_5_exp_pade_approximations, exp_pade):
        assert (p_true - num/denom).simplify() == 0

    # test sin
    all_order_6_sin_pade_approximations = [
    x**5/120 - x**3/6 + x,
    x*(60 - 7*x**2)/(3*(x**2 + 20)),
    360*x/(7*x**4 + 60*x**2 + 360)
    ]

    sin_pade = pade_approximants(sin(x), x, 6)

    for p_true, (num, denom) in zip(all_order_6_sin_pade_approximations, sin_pade):
        assert (p_true - num/denom).simplify() == 0


def test_pade_approximant():
    x = symbols("x")

    assert pade_approximant(exp(x), x, 2) == (x**2/4 + 3*x/2 + 3)/(x**2/4 - 3*x/2 + 3)
    assert pade_approximant(exp(sin(x**3)), x, 3, 4) == (x**3/2 + 1)/(1 - x**3/2)
    assert pade_approximant(sin(x), x, 0, 4) is None
    assert (pade_approximant(log(1 + x), x, 3, 3)
             - x*(11*x**2 + 60*x + 60)/(3*(x**3 + 12*x**2 + 30*x + 20))).simplify() == 0
