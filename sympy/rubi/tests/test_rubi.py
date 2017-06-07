from sympy.core.symbol import symbols
from sympy.rubi.rubi import rubi_integrate

a, b, c, d, e, f, x = symbols('a b c d e f x')

def test_rubi_integrate():
    expr = a
    assert rubi_integrate(expr, x) == x**(a + 1)/(a + 1)
