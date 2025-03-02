from sympy.core.numbers import Infinity, NegativeInfinity, ComplexInfinity
from sympy.assumptions.ask import Q
from sympy.assumptions.refine import refine
from sympy.core.numbers import (I, Rational, nan, pi)
from sympy.core import S
from sympy.core.symbol import Symbol
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.abc import w, x, y, z

def test_refine_infinity():
    expr = (sqrt(z) + S.Infinity)**2
    result = refine(expr)
    assert result == S.Infinity, f"Expected oo, but got {result}"

def test_refine_infinity_with_finite_term():
    z = Symbol('z', real=True)
    expr = S.Infinity + sqrt(z) + 5
    result = refine(expr)
    assert result == S.Infinity, f"Expected Infinity, got {result}"

def test_refine_negative_infinity():
    z = Symbol('z', real=True)
    expr = S.NegativeInfinity + sqrt(z)
    result = refine(expr)
    assert result == S.NegativeInfinity, f"Expected NegativeInfinity, got {result}"

def test_refine_add_nan():
    z = Symbol('z', real=True)
    expr = sqrt(z) + nan
    result = refine(expr)
    assert result is nan, f"Expected nan, got {result}"

def test_refine_add_complex_infinity():
    z = Symbol('z', real=True)
    expr = S.ComplexInfinity + sqrt(z)
    result = refine(expr)
    assert result == S.ComplexInfinity, f"Expected zoo, got {result}"

def test_refine_add_both_infinities():
    z = Symbol('z', real=True)
    expr = S.Infinity + S.NegativeInfinity + sqrt(z)
    result = refine(expr)
    assert result is nan, f"Expected nan, got {result}"

def test_refine_add_positive_infinity():
    z = Symbol('z', real=True)
    expr = S.Infinity + sqrt(z)
    result = refine(expr)
    assert result == S.Infinity, f"Expected Infinity, got {result}"

def test_refine_add_negative_infinity():
    z = Symbol('z', real=True)
    expr = S.NegativeInfinity + sqrt(z)
    result = refine(expr)
    assert result == S.NegativeInfinity, f"Expected NegativeInfinity, got {result}"

def test_refine_directed_infinity():
    from sympy import I, S, sqrt, refine, Symbol
    z = Symbol('z', real=True)
    expr = I * S.Infinity + sqrt(z)
    result = refine(expr)
    assert result == S.ComplexInfinity, f"Expected ComplexInfinity, got {result}"

def test_refine_infinite_symbol():
    from sympy import S, refine, Symbol
    x = Symbol('x', real=True)
    x._assumptions = x._assumptions.copy()
    x._assumptions['infinite'] = True
    expr = x + 1
    result = refine(expr)
    assert result == S.Infinity, f"Expected Infinity for an infinite symbol, got {result}"

def test_refine_accumulation_bounds_with_infinity():
    from sympy import S, refine
    from sympy.calculus.accumulationbounds import AccumulationBounds
    expr = S.Infinity + AccumulationBounds(1, 2)
    result = refine(expr)
    assert result == S.Infinity, f"Expected Infinity, got {result}"
