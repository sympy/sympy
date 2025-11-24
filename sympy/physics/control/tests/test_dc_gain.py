"""Tests for TransferFunction.dc_gain."""

from sympy import symbols, Rational, oo
from sympy.physics.control.lti import TransferFunction


def test_dc_gain_constant_num_den():
    s = symbols('s')
    G = TransferFunction(5, s + 2, s)
    assert G.dc_gain() == Rational(5, 2)


def test_dc_gain_polynomial_num_den():
    s = symbols('s')
    G = TransferFunction(s**2 + 3*s + 2, (s + 1)*(s + 2), s)
    assert G.dc_gain() == 1


def test_dc_gain_zero_numerator():
    s = symbols('s')
    G = TransferFunction(0, s + 5, s)
    assert G.dc_gain() == 0


def test_dc_gain_den_zero_at_zero_returns_infinity():
    s = symbols('s')
    G = TransferFunction(1, s, s)
    assert G.dc_gain() == oo


def test_dc_gain_negative_power_gives_infinity():
    s = symbols('s')
    G = TransferFunction(1/s, 1, s)
    assert G.dc_gain() == oo
