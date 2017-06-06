from sympy.rubi.utility_function import *
from sympy.core.symbol import symbols, S

a, b, c, d, x = symbols('a b c d x')

def test_NonzeroQ():
    assert NonzeroQ(S(1)) == True

def test_FreeQ():
    l = [a*b, x, a + b]
    assert FreeQ(l, x) == False

    l = [a*b, a + b]
    assert FreeQ(l, x) == True
