from sympy import S, Symbol, sin, cos
from sympy.integrals.rubi.SineIntegrationRules import intsin5

def test_intsin5():
    x = Symbol("x")
    assert intsin5(S(1), S(2), S(3), S(2), x) == 9*x/2 - 9*sin(2*x + 1)*cos(2*x + 1)/4
    assert intsin5(S(1), S(2), S(3), -S(2), x) == -cos(2*x + 1)/(18*sin(2*x + 1))
