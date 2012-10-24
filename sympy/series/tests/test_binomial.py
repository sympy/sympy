from sympy.series import binomial_expand
from sympy.abc import x, y
from sympy import sin, Pow

def test_binomial():
    """
    Testing function of binomial_expand(function).
    """
    assert binomial_expand(sin(x) ** 3) == Pow(sin(x), 3)
    assert binomial_expand(1) == 1
    assert binomial_expand((1+x)** 0.5).subs(x, x ** 2) == binomial_expand((1 + x ** 2) ** 0.5)
