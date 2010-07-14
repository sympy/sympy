from sympy import sin, cos, abs, exp, pi, oo, symbols, ceiling, raises, sqrt
from sympy import Function, Piecewise, Rational, Integer, GoldenRatio, EulerGamma, Catalan

from sympy.printing.ccode import ccode, CCodePrinter
from sympy.utilities.pytest import XFAIL

x, y, z = symbols('xyz')
g = Function('g')

def test_printmethod():
    class fabs(abs):
        def _ccode(self, printer):
            return "fabs(%s)" % printer._print(self.args[0])
    assert ccode(fabs(x)) == "fabs(x)"

def test_ccode_sqrt():
    assert ccode(sqrt(x)) == "sqrt(x)"
    assert ccode(x**0.5) == "sqrt(x)"
    assert ccode(x**Rational(1,2)) == "sqrt(x)"

def test_ccode_Pow():
    assert ccode(x**3) == "pow(x, 3)"
    assert ccode(x**(y**3)) == "pow(x, pow(y, 3))"
    assert ccode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "pow(3.5*g(x), -x + pow(y, x))/(y + pow(x, 2))"

def test_ccode_constants_mathh():
    assert ccode(exp(1)) == "M_E"
    assert ccode(pi) == "M_PI"
    assert ccode(oo) == "HUGE_VAL"
    assert ccode(-oo) == "-HUGE_VAL"

def test_ccode_constants_other():
    assert ccode(2*GoldenRatio) == "double const GoldenRatio = 1.61803398874989\n2*GoldenRatio"
    assert ccode(2*Catalan) == "double const Catalan = 0.915965594177219\n2*Catalan"
    assert ccode(2*EulerGamma) == "double const EulerGamma = 0.577215664901533\n2*EulerGamma"

def test_ccode_Rational():
    assert ccode(Rational(3,7)) == "3.0/7.0"
    assert ccode(Rational(18,9)) == "2"
    assert ccode(Rational(3,-7)) == "-3.0/7.0"
    assert ccode(Rational(-3,-7)) == "3.0/7.0"

def test_ccode_Integer():
    assert ccode(Integer(67)) == "67"
    assert ccode(Integer(-1)) == "-1"

def test_ccode_functions():
    assert ccode(sin(x) ** cos(x)) == "pow(sin(x), cos(x))"

def test_ccode_exceptions():
    assert ccode(ceiling(x)) == "ceil(x)"
    assert ccode(abs(x)) == "fabs(x)"

def test_ccode_boolean():
    assert ccode(x&y) == "x&&y"
    assert ccode(x|y) == "x||y"
    assert ccode(~x) == "!x"
    assert ccode(x&y&z) == "x&&y&&z"
    assert ccode(x|y|z) == "x||y||z"
    assert ccode((x&y)|z) == "x&&y||z"
    assert ccode((x|y)&z) == "(x||y)&&z"

def test_ccode_Piecewise():
    p = ccode(Piecewise((x,x<1),(x**2,True)))
    s = \
"""\
if (x < 1) {
   x
}
else {
   pow(x, 2)
}\
"""
    assert p == s

def test_ccode_Piecewise_deep():
    p = ccode(2*Piecewise((x,x<1),(x**2,True)))
    s = \
"""\
if (x < 1) {
   2*x
}
else {
   2*pow(x, 2)
}\
"""
    assert p == s

def test_ccode_settings():
    raises(TypeError, 'ccode(sin(x),method="garbage")')

def test_ccode_Indexed():
    from sympy.tensor import Indexed, Idx
    from sympy import symbols
    i,j,k,n,m,o = symbols('i j k n m o', integer=True)

    p = CCodePrinter()
    p._not_c = set()

    x = Indexed('x')(Idx(j, n))
    assert p._print_IndexedElement(x) == 'x[j]'
    A = Indexed('A')(Idx(i, m), Idx(j, n))
    assert p._print_IndexedElement(A) == 'A[%s]'% str(j + n*i)
    B = Indexed('B')(Idx(i, m), Idx(j, n), Idx(k, o))
    assert p._print_IndexedElement(B) == 'B[%s]'% str(k + i*n*o + j*o)

    assert p._not_c == set([Indexed('x'), Indexed('A'), Indexed('B')])
