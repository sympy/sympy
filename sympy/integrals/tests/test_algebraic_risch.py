from __future__ import annotations

from sympy import symbols, sqrt, oo, S, log, cancel
from sympy.integrals.algebraic_risch import (
    Place, Divisor, valuation, coates_torsion_divisor,
    integrate_algebraic_risch
)

def test_place_valuation():
    x, y = symbols('x y')
    Q_x = x**3 - x
    # Non-branch point: place P = (2, sqrt(6))
    P = Place(S(2), sqrt(S(6)))
    assert valuation(x - 2, P, x, y, Q_x) == 1
    assert valuation(y - sqrt(S(6)), P, x, y, Q_x) == 1

    # Branch point: place P = (1, 0)
    P_br = Place(S(1), S.Zero)
    assert valuation(x - 1, P_br, x, y, Q_x) == 1
    # y = sqrt(x**3 - x) = sqrt((x-1)(x**2+x)). Since x-1 = u**2, y is proportional to u.
    # Valuation in u is 1, so valuation is 1/2.
    assert valuation(y, P_br, x, y, Q_x) == S.Half


def test_divisor_arithmetic():
    P1 = Place(S.Zero, S.Zero)
    P2 = Place(S.One, S.Zero)

    div1 = Divisor({P1: 1, P2: -1})
    div2 = Divisor({P1: 2, P2: 3})

    assert div1 + div2 == Divisor({P1: 3, P2: 2})
    assert div1 - div2 == Divisor({P1: -1, P2: -4})
    assert 3 * div1 == Divisor({P1: 3, P2: -3})


def test_coates_torsion_divisor():
    x, y = symbols('x y')
    Q_x = x**2 + 1
    P_inf_1 = Place(oo, oo)
    P_inf_2 = Place(oo, -oo)

    div = Divisor({P_inf_2: 1, P_inf_1: -1})
    res = coates_torsion_divisor(div, x, y, Q_x)
    assert res is not None
    m, g = res
    assert m == 1
    # g should be proportional to x + y
    assert cancel(g / (x + y)) == 1 or cancel(g / (x + y)) == -1


def test_integrate_algebraic_risch():
    x, y = symbols('x y')

    # 1. Purely algebraic differential: (3*x**2 + 1)/(2*y) on y**2 = x**3 + x
    # Integral is y
    Q1 = x**3 + x
    res1 = integrate_algebraic_risch((3*x**2 + 1)/(2*y), x, y, Q1)
    assert res1 == y

    # 2. Logarithmic differential: 1/y on y**2 = x**2 + 1
    # Integral is log(x + y)
    Q2 = x**2 + 1
    res2 = integrate_algebraic_risch(1/y, x, y, Q2)
    assert res2 == log(x + y) or res2 == -log(x + y) or res2 == log(1/(x + y))
