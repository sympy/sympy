from sympy import abs, exp, Function, symbols

from sympy.printing import ccode
from sympy.utilities.pytest import XFAIL

x, y = symbols('xy')
g = Function('g')

def test_Pow():
    assert ccode(x**3) == "pow(x,3)"
    assert ccode(x**(y**3)) == "pow(x,(pow(y,3)))"
    assert ccode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "1/(y + pow(x,2))*pow((3.50000000000000*g(x)),(-x + pow(y,x)))"

def test_Exp1():
    assert ccode(exp(1)) == "exp(1)"
