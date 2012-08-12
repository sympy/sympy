from sympy.utilities.pytest import raises
from sympy import Symbol, symbols, Function, Integer, Matrix, nan, oo, Abs, \
    Rational, Float, S, WildFunction, ImmutableMatrix
from sympy.geometry import Point, Circle, Ellipse
from sympy.printing import srepr

x, y = symbols('x,y')

# eval(srepr(expr)) == expr has to succeed in the right environment. The right
# environment is the scope of "from sympy import *" for most cases.
ENV = {}
exec "from sympy import *" in ENV

def sT(expr, string):
    """
    sT := sreprTest

    Tests that srepr delivers the expected string and that
    the condition eval(srepr(expr))==expr holds.
    """
    assert srepr(expr) == string
    assert eval(string, ENV) == expr

def test_printmethod():
    class R(Abs):
        def _sympyrepr(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert srepr(R(x)) == "foo(Symbol('x'))"

def test_Add():
    sT(x+y, "Add(Symbol('x'), Symbol('y'))")
    assert srepr(x**2 + 1, order='lex') == "Add(Pow(Symbol('x'), Integer(2)), Integer(1))"
    assert srepr(x**2 + 1, order='old') == "Add(Integer(1), Pow(Symbol('x'), Integer(2)))"

def test_Function():
    sT(Function("f")(x), "Function('f')(Symbol('x'))")
    # test unapplied Function
    sT(Function('f'), "Function('f')")

def test_Geometry():
    sT(Point(0,0),  "Point(Integer(0), Integer(0))")
    sT(Ellipse(Point(0, 0), 5, 1),  "Ellipse(Point(Integer(0), Integer(0)), Integer(5), Integer(1))")
    # TODO more tests

def test_Singletons():
    sT(S.Catalan, 'Catalan')
    sT(S.ComplexInfinity, 'zoo')
    sT(S.EulerGamma, 'EulerGamma')
    sT(S.Exp1, 'E')
    sT(S.GoldenRatio, 'GoldenRatio')
    sT(S.Half, 'Rational(1, 2)')
    sT(S.ImaginaryUnit, 'I')
    sT(S.Infinity, 'oo')
    sT(S.NaN, 'nan')
    sT(S.NegativeInfinity, '-oo')
    sT(S.NegativeOne, 'Integer(-1)')
    sT(S.One, 'Integer(1)')
    sT(S.Pi, 'pi')
    sT(S.Zero, 'Integer(0)')

def test_Integer():
    sT(Integer(4), "Integer(4)")

def test_list():
    sT([x, Integer(4)], "[Symbol('x'), Integer(4)]")

def test_Matrix():
    # Matrix is really MutableMatrix
    for cls, name in [(Matrix, "Matrix"), (ImmutableMatrix, "ImmutableMatrix")]:
        sT(cls([[x**+1, 1], [y, x+y]]),
           "%s([[Symbol('x'), Integer(1)], [Symbol('y'), Add(Symbol('x'), Symbol('y'))]])"%name)

        sT(cls(), "%s([])"%name)

        sT(cls([[x**+1, 1], [y, x+y]]), "%s([[Symbol('x'), Integer(1)], [Symbol('y'), Add(Symbol('x'), Symbol('y'))]])"%name)

def test_Rational():
    sT(Rational(1,3), "Rational(1, 3)")
    sT(Rational(-1,3), "Rational(-1, 3)")

def test_Float():
    sT(Float('1.23', prec=3), "Float('1.22998', prec=3)")
    sT(Float('1.23456789', prec=9), "Float('1.23456788994', prec=9)")
    sT(Float('1.234567890123456789', prec=19), "Float('1.234567890123456789013', prec=19)")
    sT(Float('0.60038617995049726', 15), "Float('0.60038617995049726', prec=15)")

def test_Symbol():
    sT(x, "Symbol('x')")
    sT(y, "Symbol('y')")

def test_tuple():
    sT((x,), "(Symbol('x'),)")
    sT((x,y), "(Symbol('x'), Symbol('y'))")

def test_WildFunction():
    sT(WildFunction('w'), "WildFunction('w')")

def test_settins():
    raises(TypeError, lambda: srepr(x, method="garbage"))

def test_Mul():
    sT(3*x**3*y, "Mul(Integer(3), Pow(Symbol('x'), Integer(3)), Symbol('y'))")
    assert srepr(3*x**3*y, order='old') == "Mul(Integer(3), Symbol('y'), Pow(Symbol('x'), Integer(3)))"
