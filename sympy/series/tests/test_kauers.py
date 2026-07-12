from __future__ import annotations
from sympy.series.kauers import finite_diff
from sympy.series.kauers import finite_diff_kauers
from sympy.abc import x, y, z, m, n, w, k
from sympy.core.numbers import pi, S
from sympy.functions.elementary.trigonometric import (cos, sin)
from sympy.concrete.summations import Sum


def test_finite_diff():
    assert finite_diff(x**2 + 2*x + 1, x) == 2*x + 3
    assert finite_diff(y**3 + 2*y**2 + 3*y + 5, y) == 3*y**2 + 7*y + 6
    assert finite_diff(z**2 - 2*z + 3, z) == 2*z - 1
    assert finite_diff(w**2 + 3*w - 2, w) == 2*w + 4
    assert finite_diff(sin(x), x, pi/6) == -sin(x) + sin(x + pi/6)
    assert finite_diff(cos(y), y, pi/3) == -cos(y) + cos(y + pi/3)
    assert finite_diff(x**2 - 2*x + 3, x, 2) == 4*x
    assert finite_diff(n**2 - 2*n + 3, n, 3) == 6*n + 3


def test_finite_diff_constant():
    assert finite_diff(S(5), x) == 0
    assert finite_diff(S(0), x) == 0
    assert finite_diff(S(-3), y) == 0


def test_finite_diff_linear():
    assert finite_diff(3*x + 2, x) == 3
    assert finite_diff(-x + 7, x) == -1


def test_finite_diff_high_degree():
    assert finite_diff(x**4, x) == 4*x**3 + 6*x**2 + 4*x + 1
    assert finite_diff(x**5, x) == 5*x**4 + 10*x**3 + 10*x**2 + 5*x + 1


def test_finite_diff_negative_increment():
    assert finite_diff(x**2, x, -1) == -2*x + 1
    assert finite_diff(x**3, x, -1) == -3*x**2 + 3*x - 1


def test_finite_diff_iterated():
    second = finite_diff(finite_diff(x**2, x), x)
    assert second == 2
    third = finite_diff(finite_diff(finite_diff(x**3, x), x), x)
    assert third == 6


def test_finite_diff_kauers():
    assert finite_diff_kauers(Sum(x**2, (x, 1, n))) == (n + 1)**2
    assert finite_diff_kauers(Sum(y, (y, 1, m))) == (m + 1)
    assert finite_diff_kauers(Sum((x*y), (x, 1, m), (y, 1, n))) == (m + 1)*(n + 1)
    assert finite_diff_kauers(Sum((x*y**2), (x, 1, m), (y, 1, n))) == (n + 1)**2*(m + 1)


def test_finite_diff_kauers_constant_summand():
    assert finite_diff_kauers(Sum(S(1), (k, 1, n))) == 1


def test_finite_diff_kauers_higher_power():
    assert finite_diff_kauers(Sum(k**3, (k, 1, n))) == (n + 1)**3


def test_finite_diff_kauers_reciprocal():
    assert finite_diff_kauers(Sum(1/k, (k, 1, n))) == S(1)/(n + 1)


def test_finite_diff_kauers_triple_nested():
    expr = Sum(x*y*z, (x, 1, m), (y, 1, n), (z, 1, k))
    assert finite_diff_kauers(expr) == (m + 1)*(n + 1)*(k + 1)
