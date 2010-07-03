from sympy.utilities.pytest import XFAIL, raises
from sympy import Symbol, symbols, Function, Integer, Matrix, nan, oo, Abs, \
    Rational, Real, S, WildFunction
from sympy.geometry import Point, Circle, Ellipse
from sympy.printing import srepr

x, y = symbols('x,y')

# eval(srepr(expr)) == expr has to succeed in the right environment. The right
# environment is the scope of "from sympy import *" for most cases.
ENV = {}
exec "from sympy import *" in ENV
# These classes have to be added separately:
ENV["Infinity"] = S.Infinity
ENV["NegativeInfinity"] = S.NegativeInfinity
ENV["NegativeOne"] = S.NegativeOne
ENV["One"] = S.One
ENV["Zero"] = S.Zero

def sT(expr, string):
    """
    sT := sreprTest

    Tests that srepr delivers the expected string and that
    the condition eval(srepr(expr))==expr holds.
    """
    assert srepr(expr) == string
    assert eval(string, ENV) == expr

def test_printmethod():
    class R(oo.__class__):
        def _sympyrepr(self, printer):
            return "foo"
    assert srepr(R()) == "foo"
    class R(Abs):
        def _sympyrepr(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert srepr(R(x)) == "foo(Symbol('x'))"

def test_Add():
    sT(x+y, "Add(Symbol('x'), Symbol('y'))")

def test_Function():
    sT(Function("f")(x), "Function('f')(Symbol('x'))")
    # test unapplied Function
    sT(Function('f'), "Function('f')")

def test_Geometry():
    sT(Point(0,0),  "Point(Zero, Zero)")
    sT(Ellipse(Point(0, 0), 5, 1),  "Ellipse(Point(Zero, Zero), Integer(5), One)")
    # TODO more tests

def test_Infinity():
    sT(oo, "Infinity")

def test_Integer():
    sT(Integer(4), "Integer(4)")

def test_list():
    sT([x, Integer(4)], "[Symbol('x'), Integer(4)]")

def test_Matrix():
    sT(Matrix([[x**+1, 1], [y, x+y]]), "Matrix([[Symbol('x'), One], [Symbol('y'), Add(Symbol('x'), Symbol('y'))]])")

def test_NaN():
    sT(nan, "nan")

def test_NegativeInfinity():
    sT(-oo, "NegativeInfinity")

def test_NegativeOne():
    sT(-Integer(1), "NegativeOne")

def test_One():
    sT(S.One, "One")

def test_Rational():
    sT(Rational(1,3), "Rational(1, 3)")
    sT(Rational(-1,3), "Rational(-1, 3)")

def test_Real():
    sT(Real('1.23', prec=3), "Real('1.22998', prec=3)")
    sT(Real('1.23456789', prec=9), "Real('1.23456788994', prec=9)")
    sT(Real('1.234567890123456789', prec=19), "Real('1.234567890123456789013', prec=19)")
    sT(Real('0.60038617995049726', 15), "Real('0.60038617995049726', prec=15)")

def test_Symbol():
    sT(x, "Symbol('x')")
    sT(y, "Symbol('y')")

def test_tuple():
    sT((x,), "(Symbol('x'),)")
    sT((x,y), "(Symbol('x'), Symbol('y'))")

def test_WildFunction():
    sT(WildFunction('w'), "WildFunction('w')")

def test_Zero():
    sT(S.Zero, "Zero")

def test_settins():
    raises(TypeError, 'srepr(x, method="garbage")')
