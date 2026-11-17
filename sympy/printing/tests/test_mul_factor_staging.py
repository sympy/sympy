"""Regression tests for staged numerator and denominator factors in printers."""

from __future__ import annotations

from sympy import Float, Mul, Pow, S
from sympy.printing import (ccode, fcode, glsl_code, jscode, julia_code,
                            maple_code, octave_code, pycode, rcode, sstr)
from sympy.abc import x, y


def _nested_reciprocal(exponent):
    inner = Pow(y, exponent, evaluate=False)
    outer = Pow(inner, -1, evaluate=False)
    return Mul(x, outer, evaluate=False)


def test_nested_reciprocal_in_code_printers():
    expr = _nested_reciprocal(-1)
    assert sstr(expr) == 'x/(1/y)'
    assert pycode(expr) == 'x/(1/y)'
    assert ccode(expr) == 'x/(1.0/y)'
    assert julia_code(expr) == 'x ./ (1 ./ y)'
    assert octave_code(expr) == 'x./(1./y)'
    assert pycode(-expr) == '-x/(1/y)'
    assert eval(pycode(expr), {'x': 3.0, 'y': 4.0}) == 12.0
    assert eval(pycode(-expr), {'x': 3.0, 'y': 4.0}) == -12.0


def test_nested_reciprocal_in_other_code_printers():
    # Root cause 2 in gh-30844 affects all of these CodePrinter subclasses.
    expr = _nested_reciprocal(-1)
    assert jscode(expr) == 'x/(1/y)'
    assert rcode(expr) == 'x/(1.0/y)'
    assert glsl_code(expr) == 'x/(1.0/y)'
    assert fcode(expr).strip() == 'x/(1d0/y)'
    assert maple_code(expr) == 'x/(1/y)'


def test_nested_reciprocal_sqrt():
    expr = _nested_reciprocal(-S.Half)
    assert pycode(expr) == 'x/(1/math.sqrt(y))'
    assert julia_code(expr) == 'x ./ (1 ./ sqrt(y))'
    assert octave_code(expr) == 'x./(1./sqrt(y))'


def test_float_reciprocal_grouping():
    expr = _nested_reciprocal(Float(-1))
    assert sstr(expr) == 'x/y**(-1.0)'
    assert pycode(expr) == 'x/y**(-1.0)'
    assert ccode(expr) == 'x/(1.0/y)'
    assert julia_code(expr) == 'x ./ (1 ./ y)'
    assert octave_code(expr) == 'x./(1./y)'


def test_float_half_reciprocal_grouping():
    expr = _nested_reciprocal(Float(-0.5))
    assert pycode(expr) == 'x/y**(-0.5)'
    assert ccode(expr) == 'x/pow(y, -0.5)'
    assert julia_code(expr) == 'x ./ (1 ./ sqrt(y))'
    assert octave_code(expr) == 'x./(1./sqrt(y))'


def test_nonspecial_negative_powers_remain_powers():
    expr = _nested_reciprocal(-2)
    assert sstr(expr) == 'x/y**(-2)'
    assert pycode(expr) == 'x/y**(-2)'
    assert ccode(expr) == 'x/pow(y, -2)'
    assert julia_code(expr) == 'x ./ y .^ (-2)'
    assert octave_code(expr) == 'x./y.^(-2)'
