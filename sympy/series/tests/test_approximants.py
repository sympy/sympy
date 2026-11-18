from __future__ import annotations
from sympy.series import approximants
from sympy.series.approximants import pade_approximants_gcdex, pade_approximant, pade_approximant_gcdex
from sympy.core.symbol import symbols
from sympy.functions.combinatorial.factorials import binomial
from sympy.functions.combinatorial.numbers import (fibonacci, lucas)
from sympy import exp, sin, cos, log, Poly, cancel
from sympy.testing.pytest import raises


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


def test_pade_approximants_gcdex():
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

    exp_series = Poly(1 + x + x**2/2 + x**3/6 + x**4/24 + x**5/120, x, domain='QQ')

    exp_pade = pade_approximants_gcdex(exp_series, 5)

    for p_true, (num, denom) in zip(all_order_5_exp_pade_approximations, exp_pade):
        assert (p_true - num/denom).simplify() == 0

    # test sin
    all_order_6_sin_pade_approximations = [
    x**5/120 - x**3/6 + x,
    x*(60 - 7*x**2)/(3*(x**2 + 20)),
    360*x/(7*x**4 + 60*x**2 + 360)
    ]

    sin_series = Poly(x - x**3/6 + x**5/120, x, domain='QQ')

    sin_pade = pade_approximants_gcdex(sin_series, 6)

    for p_true, (num, denom) in zip(all_order_6_sin_pade_approximations, sin_pade):
        assert (p_true - num/denom).simplify() == 0


def test_pade_approximant_gcdex():
    x = symbols('x')
    poly_1 = Poly(1 + x + x**2/2 + x**3/6 + x**4/24, x, domain='QQ')
    poly_2 = Poly(1 + x**3 + x**6/2, x, domain='QQ')
    poly_3 = Poly(x - x**3/6 + x**5/120, x, domain='QQ')
    poly_4 = Poly(2 + 2*x/5 - 2*x**2/25 + 8*x**3/375
                   - 4*x**4/625 + 32*x**5/15625 - 32*x**6/46875, x, domain='QQ')

    numerator, denominator = pade_approximant_gcdex(poly_1, 2)
    assert numerator is not None and denominator is not None
    assert numerator/numerator.LC() == x**2 + 6*x + 12
    assert (numerator/denominator).equals((x**2/4 + 3*x/2 + 3)/(x**2/4 - 3*x/2 + 3))

    numerator, denominator = pade_approximant_gcdex(poly_2, 3, 4)
    assert numerator is not None and denominator is not None
    assert numerator/numerator.LC() == x**3 + 2
    assert (numerator/denominator).equals((x**3/2 + 1)/(1 - x**3/2))

    numerator, denominator = pade_approximant_gcdex(poly_4, 3)
    assert numerator is not None and denominator is not None
    assert numerator/numerator.LC() - (x**3 + 330*x**2/17 + 1500*x/17 + 1875/17) == 0
    assert (numerator/denominator - (2 + 8*x/5 + 44*x**2/125 + 34*x**3/1875)\
                                    /(1 + 3*x/5 + 12*x**2/125 + 2*x**3/625)).simplify() == 0

    with raises(ValueError):
        pade_approximant_gcdex(poly_3, 0, 4)


def test_pade_approximant():
    x = symbols("x")

    assert pade_approximant(exp(x), x, 0, 2)\
            == (x**2/4 + 3*x/2 + 3)/(x**2/4 - 3*x/2 + 3)

    assert pade_approximant(exp(sin(x**3)), x, 0, 3, 4)\
            == (x**3/2 + 1)/(1 - x**3/2)

    assert (
        pade_approximant(log(1 + x), x, 0, 3) - \
        x*(11*x**2 + 60*x + 60)/(3*(x**3 + 12*x**2 + 30*x + 20))
        ).simplify() == 0

    assert (
        pade_approximant(log(x), x, 1, 3) - \
        (x-1)*(11*(x-1)**2 + 60*(x-1) + 60)/(3*((x-1)**3 + 12*(x-1)**2 + 30*(x-1) + 20))
        ).simplify() == 0

    assert (
        pade_approximant(exp(x)/x, x, 0, 2, 3) - \
        (x**2/4 + 3*x/2 + 3)/(x*(x**2/4 - 3*x/2 + 3))
    ).simplify() == 0

    with raises(ValueError):
        pade_approximant(sin(x), x, 0, 0, 4)

    with raises(ValueError):
        pade_approximant(1/x, x, 0, 1, 0)


def test_pade_approximant_gcdex_degree_jump():
    # The EEA degree sequence can skip m: for 1/(1 + x) truncated at order 2
    # it is 2, 0, so no convergent has numerator degree exactly m = n = 1.
    # The [1/1] approximant is still 1/(1 + x) and must be the one returned.
    x = symbols("x")
    f = Poly(1 - x + x**2, x, domain='QQ')
    num, den = pade_approximant_gcdex(f, 1, 1)
    assert num.degree() <= 1
    assert den.degree() <= 1
    assert (num/den - 1/(1 + x)).simplify() == 0


def test_pade_approximant_degree_bounds():
    x = symbols("x")
    for m, n in [(1, 1), (2, 1), (1, 2), (3, 2), (2, 3), (5, 1)]:
        approx = pade_approximant(1/(1 + x), x, 0, m, n)
        # the expression as returned, before any cancellation
        num, den = approx.as_numer_denom()
        assert Poly(num, x).degree() <= m
        assert Poly(den, x).degree() <= n
        approx = cancel(approx)
        num, den = approx.as_numer_denom()
        assert Poly(num, x).degree() <= m
        assert Poly(den, x).degree() <= n
        assert (approx - 1/(1 + x)).simplify() == 0

    # expansion around a non-zero center
    t = symbols("t")
    approx = cancel(pade_approximant(1/(1 + x), x, 1, 1, 1).subs(x, t + 1))
    num, den = approx.as_numer_denom()
    assert Poly(num, t).degree() <= 1
    assert Poly(den, t).degree() <= 1
    assert (approx - 1/(t + 2)).simplify() == 0


def test_pade_approximant_no_solution_degree_bounds():
    # No [2/1] approximant exists for 1 + x**3, but the fallback returned by
    # the Euclidean algorithm still has to respect both degree bounds.
    x = symbols("x")
    num, den = pade_approximant_gcdex(Poly(1 + x**3, x), 2, 1)
    assert num.degree() <= 2
    assert den.degree() <= 1
    assert den.eval(0) == 0
    assert num/den == 1

    # cos(x) has no [1/1] approximant; the documented fallback value is 1.
    assert pade_approximant(cos(x), x, 0, 1) == 1


def test_pade_approximant_gcdex_zero_polynomial():
    # The zero polynomial already meets both degree bounds; its first EEA
    # candidate is 0/1.
    x = symbols("x")
    num, den = pade_approximant_gcdex(Poly(0, x), 1, 1)
    assert num.is_zero
    assert den.as_expr() == 1
