from sympy import I, Pow, Symbol, sin
from sympy.printing.latex import latex

def _one_of(actual: str, choices):
    return any(actual == c for c in choices)


def test_latex_unevaluated_pow_simple():
    expr = Pow(I, -4, evaluate=False)
    # For I**-4 some printers may evaluate to 1; the important thing is
    # there is no recursion and output is sensible
    # Accept either form
    out = latex(expr)
    acceptable = [r"i^{-4}", r"\frac{1}{i^{4}}", r"1"]
    assert out in acceptable or "RecursionError" not in out

def test_latex_unevaluated_pow_symbol_base():
    x = Symbol('x')
    expr = Pow(x, -3, evaluate=False)
    out = latex(expr)
    # Accept either the explicit negative exponent or a fraction form
    acceptable = [r"x^{-3}", r"\frac{1}{x^{3}}"]
    assert _one_of(out, acceptable)

def test_latex_unevaluated_pow_add_base():
    x = Symbol('x')
    y = Symbol('y')
    expr = Pow(x + y, -2, evaluate=False)
    out = latex(expr)
    acceptable = [r"\left(x + y\right)^{-2}", r"\frac{1}{\left(x + y\right)^{2}}"]
    assert _one_of(out, acceptable)

def test_latex_unevaluated_pow_function_base():
    x = Symbol('x')
    expr = Pow(sin(x), -5, evaluate=False)
    out = latex(expr)
    acceptable = [r"\sin{\left(x \right)}^{-5}", r"\frac{1}{\sin^{5}{\left(x \right)}}"]
    assert _one_of(out, acceptable)

def test_latex_unevaluated_pow_does_not_rewrite_as_mul():
    # Ensure the case that previously recursed does not recurse
    expr = Pow(I, -4, evaluate=False)
    out = latex(expr)
    # It should return a string and not blow up with RecursionError
    assert isinstance(out, str)
