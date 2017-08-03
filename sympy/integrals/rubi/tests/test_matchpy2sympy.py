import sys
from sympy.external import import_module
matchpy = import_module("matchpy")

if not matchpy:
    #bin/test will not execute any tests now
    disabled = True

if sys.version_info[:2] < (3, 6):
    disabled = True

from sympy.integrals.rubi.matchpy2sympy import matchpy2sympy
from sympy.integrals.rubi.symbol import VariableSymbol, Integer
from sympy import symbols, S, I
from sympy.integrals.rubi.operation import List, LinearQ, Add, Mul, Pow, PosQ

a, b, x = symbols('a b x', real=True, imaginary=False)

def test_matchpy2sympy():
    expr = VariableSymbol('a')
    assert matchpy2sympy(expr) == a
    expr = Integer('1')
    assert matchpy2sympy(expr) == S(1)
    expr = Integer('-1')
    assert matchpy2sympy(expr) == -S(1)
    expr = VariableSymbol('I')
    assert matchpy2sympy(expr) == I
    expr = LinearQ(VariableSymbol('a'), VariableSymbol('x'))
    assert matchpy2sympy(expr) == False
    expr = Mul(Integer('2'), VariableSymbol('x'))
    assert matchpy2sympy(expr) == S(2)*x
    expr = LinearQ(Add(Mul(Integer('2'), VariableSymbol('x')), Pow(VariableSymbol('y'), Integer(2))), VariableSymbol('x'))
    assert matchpy2sympy(expr) == True
    expr = PosQ(Integer('0'))
    assert matchpy2sympy(expr) == False
