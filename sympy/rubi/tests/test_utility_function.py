from sympy.rubi.utility_function import *
from sympy.core.symbol import symbols, S
from sympy.functions import log, sin, cos

a, b, c, d, x = symbols('a b c d x')

def test_NonzeroQ():
    assert NonzeroQ(S(1)) == True

def test_FreeQ():
    l = [a*b, x, a + b]
    assert FreeQ(l, x) == False

    l = [a*b, a + b]
    assert FreeQ(l, x) == True

def test_List():
    assert List(a, b, c) == [a, b, c]

def test_Log():
    assert Log(a) == log(a)
