from sympy import abs, exp, Function, Piecewise, symbols

from sympy.printing import ccode
from sympy.utilities.pytest import XFAIL

x, y = symbols('xy')
g = Function('g')

def test_ccode_Pow():
    assert ccode(x**3) == "pow(x,3)"
    assert ccode(x**(y**3)) == "pow(x,(pow(y,3)))"
    assert ccode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "1/(y + pow(x,2))*pow((3.50000000000000*g(x)),(-x + pow(y,x)))"

def test_ccode_Exp1():
    assert ccode(exp(1)) == "exp(1)"

def test_ccode_Piecewise():
    p = ccode(Piecewise((x,x<1),(x**2,True)))
    s = \
"""\
if (x < 1) {
x
}
else {
pow(x,2)
}\
"""
    assert p == s
