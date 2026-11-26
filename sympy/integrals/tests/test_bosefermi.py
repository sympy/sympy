"""Tests for the Bose-Einstein and Fermi-Dirac definite integrals."""
from __future__ import annotations

from sympy import (
    S, exp, gamma, integrate, Integral, log, oo, pi, simplify, symbols, zeta)
from sympy.functions.special.zeta_functions import dirichlet_eta, polylog
from sympy.integrals.bosefermi import bose_fermi_definite

x = symbols('x', positive=True)
a, p = symbols('a p', positive=True)


def test_fermi_dirac_zero_offset():
    assert integrate(1/(exp(x) + 1), (x, 0, oo)) == log(2)
    assert integrate(x**2/(exp(x) + 1), (x, 0, oo)) == gamma(3)*dirichlet_eta(3)
    assert simplify(
        integrate(x**p/(exp(x) + 1), (x, 0, oo))
        - gamma(p + 1)*dirichlet_eta(p + 1)) == 0


def test_bose_einstein_zero_offset():
    assert integrate(x**3/(exp(x) - 1), (x, 0, oo)) == pi**4/15
    assert integrate(x**2/(exp(x) - 1), (x, 0, oo)) == 2*zeta(3)
    assert simplify(
        integrate(x**p/(exp(a*x) - 1), (x, 0, oo))
        - gamma(p + 1)*zeta(p + 1)/a**(p + 1)) == 0


def test_nonzero_offset():
    assert simplify(
        integrate(x**2/(exp(x + 1) + 1), (x, 0, oo))
        + gamma(3)*polylog(3, -exp(-1))) == 0
    assert simplify(
        integrate(x**2/(exp(x + 1) - 1), (x, 0, oo))
        - gamma(3)*polylog(3, exp(-1))) == 0
    assert integrate(1/(exp(x + 1) + 1), (x, 0, oo)) == log(1 + exp(-1))


def test_scaled_denominators():
    assert integrate(x/(2*(exp(x) - 1)), (x, 0, oo)) == pi**2/12
    assert integrate(x/(exp(x)/2 + 1), (x, 0, oo)) == -polylog(2, -S(2))
    assert integrate(x**2/(exp(x) + 2), (x, 0, oo)) == -polylog(3, -S(2))


def test_float_offset():
    result = integrate(x**2/(exp(2*x + 0.3) + 1), (x, 0, oo))
    assert abs(result.evalf() - 0.17095710982756174) < 1e-10


def test_unsupported_forms():
    assert bose_fermi_definite(x**2/(exp(x**2) - 1), x) is None
    assert bose_fermi_definite(x/(x**2 + 1), x) is None
    assert bose_fermi_definite(exp(-x), x) is None


def test_divergent_integrals_are_rejected():
    # p = 0 for bosons and a negative offset for bosons both diverge
    assert bose_fermi_definite(1/(exp(x) - 1), x) is None
    assert bose_fermi_definite(x/(exp(x - 1) - 1), x) is None


def test_meijerg_false_leaves_the_integral_unevaluated():
    result = integrate(x**3/(exp(x) - 1), (x, 0, oo), meijerg=False)
    assert isinstance(result, Integral)
    assert simplify(result.function - x**3/(exp(x) - 1)) == 0


def test_conditions_are_returned():
    result, cond = integrate(
        x**p/(exp(a*x) - 1), (x, 0, oo), conds='separate')
    assert simplify(result - gamma(p + 1)*zeta(p + 1)/a**(p + 1)) == 0
    assert cond == True
