from sympy import sqrt, oo, refine, S
from sympy.abc import z

def test_refine_infinity():
    expr = (sqrt(z) + oo)**2

    result = refine(expr)

    assert result == oo, f"Expected oo, but got {result}"

def test_refine_infinity_with_finite_term():
    expr = sqrt(z) + 5 + oo
    result = refine(expr)

    assert result == oo, f"Expected oo, but got {result}"

def test_refine_negative_infinity():
    expr = sqrt(z) - oo
    result = refine(expr)
    assert result == S.NegativeInfinity
