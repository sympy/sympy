from __future__ import annotations

import sympy as sp
from sympy import Rational
from sympy.printing.latex import latex
from sympy.printing.pretty import pretty
from sympy.printing.mathml import mathml


def test_latex_root_notation_false():
    expr = sp.Symbol('x') ** Rational(1, 3)
    assert latex(expr, root_notation=False) == r"x^{\frac{1}{3}}"


def test_latex_root_notation_true():
    expr = sp.Symbol('x') ** Rational(1, 3)
    # Default includes root notation, which renders as a radical
    # The exact output may vary; we just ensure it's not the exponent form
    assert "sqrt" in latex(expr, root_notation=True)


def test_pretty_root_notation_false():
    expr = sp.Symbol('x') ** Rational(1, 3)
    # pretty prints without radicals when root_notation=False
    out = pretty(expr, root_notation=False)
    # Should contain a fraction representation rather than a radical
    assert "1/3" in out and "sqrt" not in out


def test_mathml_root_notation_false():
    expr = sp.Symbol('x') ** Rational(1, 3)
    out = mathml(expr, root_notation=False)
    # Ensure no <msqrt> tag appears
    assert "<msqrt>" not in out
