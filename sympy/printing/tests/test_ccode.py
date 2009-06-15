from sympy import sin, cos, abs, exp, pi, oo, symbols, ceiling
from sympy import Function, Piecewise, Rational, Integer

from sympy.printing import ccode
from sympy.utilities.pytest import XFAIL

x, y = symbols('xy')
g = Function('g')

def test_printmethod():
    class fabs(abs):
        def _ccode_(self):
            return "fabs(%s)" % ccode(self.args[0])
    assert ccode(fabs(x)) == "fabs(x)"

def test_ccode_Pow():
    assert ccode(x**3) == "pow(x,3)"
    assert ccode(x**(y**3)) == "pow(x,(pow(y,3)))"
    assert ccode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "pow((3.5*g(x)),(-x + pow(y,x)))/(y + pow(x,2))"

def test_ccode_constants():
    assert ccode(exp(1)) == "M_E"
    assert ccode(pi) == "M_PI"
    assert ccode(oo) == "HUGE_VAL"
    assert ccode(-oo) == "-HUGE_VAL"

def test_ccode_Rational():
    assert ccode(Rational(3,7)) == "3.0/7.0"
    assert ccode(Rational(18,9)) == "2"
    assert ccode(Rational(3,-7)) == "-3.0/7.0"
    assert ccode(Rational(-3,-7)) == "3.0/7.0"

def test_ccode_Integer():
    assert ccode(Integer(67)) == "67"
    assert ccode(Integer(-1)) == "-1"

def test_ccode_functions():
    assert ccode(sin(x) ** cos(x)) == "pow(sin(x),cos(x))"

def test_ccode_exceptions():
    assert ccode(ceiling(x)) == "ceil(x)"
    assert ccode(abs(x)) == "fabs(x)"

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
