from sympy.core.symbol import symbols
from sympy.rubi.rubi import rubi_integrate
from sympy.functions import log

a, b, c, d, e, f, x = symbols('a b c d e f x')

def test_rubi_integrate():
    expr = x**a
    assert rubi_integrate(expr, x) == x**(a + 1)/(a + 1)
    expr = x
    assert rubi_integrate(expr, x) == (1/2)*x**2
    expr = 1/x
    assert rubi_integrate(expr, x) == log(x)
