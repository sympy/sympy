from __future__ import annotations

from sympy import (And, Eq, I, Integral, Ne, Piecewise, Rational, S, cos,
                   exp, integrate, oo, pi, sign, sin, sqrt, symbols)
from sympy.integrals.trigonometry import trigintegrate_definite
from sympy.testing.pytest import raises


x = symbols('x')


def test_issue_20370():
    a = symbols('a', positive=True)
    result = integrate(1/(1 + a*cos(x)), (x, 0, 2*pi))
    assert result == Piecewise(
        (2*pi/sqrt(1 - a**2), a**2 < 1),
        (oo, Eq(a**2, 1)), (S.NaN, True))
    assert result.subs(a, Rational(1, 2)) == 4*sqrt(3)*pi/3
    assert result.subs(a, 1) is oo
    assert result.subs(a, 2) is S.NaN


def test_real_parameter_and_specialization():
    a = symbols('a', real=True)
    result = integrate(1/(1 + a*cos(x)), (x, 0, 2*pi))
    for value in [-2, -1, -S.Half, 0, S.Half, 1, 2]:
        assert result.subs(a, value) == integrate(
            1/(1 + value*cos(x)), (x, 0, 2*pi))


def test_affine_sine_cosine_denominator():
    a, b, c = symbols('a b c', real=True)
    delta = a**2 - b**2 - c**2
    result = integrate(1/(a + b*cos(x) + c*sin(x)), (x, 0, 2*pi))
    assert result == Piecewise(
        (2*pi*sign(a)/sqrt(delta), delta > 0),
        (oo*sign(a), And(Eq(delta, 0), Ne(a, 0))), (S.NaN, True))
    assert result.subs({a: 6, b: 3, c: 4}) == 2*pi/sqrt(11)
    assert result.subs({a: -6, b: 3, c: 4}) == -2*pi/sqrt(11)
    assert result.subs({a: 5, b: 3, c: 4}) is oo
    assert result.subs({a: -5, b: 3, c: 4}) is -oo
    assert result.subs({a: 4, b: 3, c: 4}) is S.NaN
    assert result.subs({a: 2, b: 0, c: 0}) == pi


def test_period_orientation_and_phase():
    phase, start = symbols('phase start', real=True)
    k = symbols('k', positive=True)
    n = symbols('n', positive=True, integer=True)
    f = 1/(6 + 3*cos(k*x + phase) + 4*sin(k*x + phase))
    end = start + 2*pi*n/k
    assert integrate(f, (x, start, end)) == 2*pi*n/(k*sqrt(11))
    assert integrate(f, (x, end, start)) == -2*pi*n/(k*sqrt(11))
    assert integrate(f.subs(k, -k), (x, start, end)) == 2*pi*n/(k*sqrt(11))
    assert integrate(1/(2 + cos(x)), (x, -pi, pi)) == 2*pi/sqrt(3)


def test_improper_integrals_and_scalar_numerator():
    d = symbols('d', real=True)
    result = integrate(d/(1 + cos(x)), (x, 0, 2*pi))
    assert result.subs(d, 0) == 0
    assert result.subs(d, 2) is oo
    assert result.subs(d, -2) is -oo
    assert integrate(1/(1 + cos(x)), (x, 2*pi, 0)) is -oo
    assert integrate(1/(-1 + sin(x)), (x, 0, 2*pi)) is -oo
    assert integrate(1/cos(x), (x, 0, 2*pi)) is S.NaN
    assert integrate(1/(1 + 2*sin(x)), (x, 0, 2*pi)) is S.NaN
    assert integrate(7/(6 + 3*cos(x) + 4*sin(x)), (x, 0, 2*pi)) == 14*pi/sqrt(11)


def test_convergence_hints():
    a = symbols('a', positive=True)
    f = 1/(1 + a*cos(x))
    assert integrate(f, (x, 0, 2*pi), conds='separate') == (
        2*pi/sqrt(1 - a**2), 1 - a**2 > 0)
    assert integrate(f, (x, 0, 2*pi), conds='none') == 2*pi/sqrt(1 - a**2)
    assert integrate(1/(2 + cos(x)), (x, 0, 2*pi), conds='separate') == (
        2*pi/sqrt(3), S.true)
    y = symbols('y')
    raises(ValueError, lambda: integrate(f, (x, 0, 2*pi), (y, 0, 1), conds='separate'))


def test_nested_integral():
    a = symbols('a', positive=True)
    assert integrate(1/(1 + a*cos(x)), (x, 0, 2*pi), (a, 0, S.Half)) == pi**2/3


def test_definite_rule_scope():
    a = symbols('a')
    r = symbols('r', real=True)
    # Unsupported inputs must reach the other integration methods.
    for f, lo, hi in [
        (exp(x), 0, 2*pi),
        (cos(x), 0, 2*pi),
        (x/(2 + cos(x)), 0, 2*pi),
        (1/(2 + cos(x)**2), 0, 2*pi),
        (1/(2 + cos(x) + sin(2*x)), 0, 2*pi),
        (1/(2 + cos(x**2)), 0, 2*pi),
        (1/(2 + cos(a*x)), 0, 2*pi),
        (1/(2 + cos(r*x)), 0, 2*pi),
        (1/(a + cos(x)), 0, 2*pi),
        (1/(2 + cos(x + I)), 0, 2*pi),
        (1/(2 + cos(x)), 0, pi),
        (1/(2 + cos(x)), 0, oo),
        (1/(2 + cos(x)), 0, 2*pi*I),
        (1/(2 + cos(x)), 0, 0),
    ]:
        assert trigintegrate_definite(f, x, S(lo), S(hi)) is None
    assert integrate(1/(2 + cos(x)), (x, 0, pi)) == pi/sqrt(3)
    f = 1/(2 + cos(x))
    assert isinstance(integrate(f, x, meijerg=True), Integral)
