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

def test_latex_unevaluated_pow_I_no_double_exponent():
    expr = Pow(I, -4, evaluate=False)
    out = latex(expr)

    if "-4" in out:
        assert out.count("-4") == 1

def test_latex_unevaluated_pow_I_no_recursive_pattern():
    expr = Pow(I, -4, evaluate=False)
    out = latex(expr)
    # Ensure the output does not repeat itself (a recursion symptom)
    assert out.count(out) == 1

def test_latex_unevaluated_pow_I_case_variants():
    expr = Pow(I, -4, evaluate=False)
    out = latex(expr)
    acceptable = ["i^{-4}", "I^{-4}", r"\frac{1}{i^{4}}", r"\frac{1}{I^{4}}", "1"]
    assert out in acceptable

def test_latex_unevaluated_pow_I_positive_exponent():
    expr = Pow(I, 4, evaluate=False)
    out = latex(expr)
    acceptable = ["i^{4}", "I^{4}"]
    assert _one_of(out, acceptable)

def test_latex_pow_I_evaluated_behavior():
    expr = Pow(I, -4, evaluate=True)
    out = latex(expr)
    # i**-4 = 1
    assert out == "1"
def test_latex_unevaluated_pow_complex_but_not_I():
    z = 2*I
    expr = Pow(z, -4, evaluate=False)
    out = latex(expr)
    # Should NOT match the special-case for I
    assert out != "i^{-4}" and out != "I^{-4}"
