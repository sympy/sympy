import sympy
from sympy.concrete.products import Product
from sympy.concrete.summations import Sum
from sympy.core.function import Derivative
from sympy.core.numbers import Integer, Rational, Float, oo
from sympy.core.relational import Rel
from sympy.core.symbol import symbols
from sympy.functions import sin
from sympy.integrals.integrals import Integral
from sympy.series.order import Order

from sympy.printing.precedence import precedence, PRECEDENCE

x, y = symbols("x,y")


def test_Add():
    assert precedence(x + y) == PRECEDENCE["Add"]
    assert precedence(x*y + 1) == PRECEDENCE["Add"]


def test_Function():
    assert precedence(sin(x)) == PRECEDENCE["Func"]

def test_Derivative():
    assert precedence(Derivative(x, y)) == PRECEDENCE["Atom"]

def test_Integral():
    assert precedence(Integral(x, y)) == PRECEDENCE["Atom"]


def test_Mul():
    assert precedence(x*y) == PRECEDENCE["Mul"]
    assert precedence(-x*y) == PRECEDENCE["Add"]

import sympy
from sympy import symbols

class F(sympy.Function):
    precedence = 40

    def _sympystr(self, printer):
        arg1 = printer.parenthesize(self.args[0], self.precedence)
        arg2 = printer.parenthesize(self.args[1], self.precedence)
        return f"{arg1} F {arg2}"

def test_custom_function_mul_precedence():
    a, b = symbols("a b")

    f_instance = F(a, b)

    positive_mul = 2 * f_instance
    assert str(positive_mul) == "2*(a F b)", f"Expected '2*(a F b)', got '{str(positive_mul)}'"

    negative_mul = -2 * f_instance
    assert str(negative_mul) == "-2*(a F b)", f"Expected '-2*(a F b)', got '{str(negative_mul)}'"



def test_Number():
    assert precedence(Integer(0)) == PRECEDENCE["Atom"]
    assert precedence(Integer(1)) == PRECEDENCE["Atom"]
    assert precedence(Integer(-1)) == PRECEDENCE["Add"]
    assert precedence(Integer(10)) == PRECEDENCE["Atom"]
    assert precedence(Rational(5, 2)) == PRECEDENCE["Mul"]
    assert precedence(Rational(-5, 2)) == PRECEDENCE["Add"]
    assert precedence(Float(5)) == PRECEDENCE["Atom"]
    assert precedence(Float(-5)) == PRECEDENCE["Add"]
    assert precedence(oo) == PRECEDENCE["Atom"]
    assert precedence(-oo) == PRECEDENCE["Add"]


def test_Order():
    assert precedence(Order(x)) == PRECEDENCE["Atom"]


def test_Pow():
    assert precedence(x**y) == PRECEDENCE["Pow"]
    assert precedence(-x**y) == PRECEDENCE["Add"]
    assert precedence(x**-y) == PRECEDENCE["Pow"]


def test_Product():
    assert precedence(Product(x, (x, y, y + 1))) == PRECEDENCE["Atom"]


def test_Relational():
    assert precedence(Rel(x + y, y, "<")) == PRECEDENCE["Relational"]


def test_Sum():
    assert precedence(Sum(x, (x, y, y + 1))) == PRECEDENCE["Atom"]


def test_Symbol():
    assert precedence(x) == PRECEDENCE["Atom"]


def test_And_Or():
    # precedence relations between logical operators, ...
    assert precedence(x & y) > precedence(x | y)
    assert precedence(~y) > precedence(x & y)
    # ... and with other operators (cfr. other programming languages)
    assert precedence(x + y) > precedence(x | y)
    assert precedence(x + y) > precedence(x & y)
    assert precedence(x*y) > precedence(x | y)
    assert precedence(x*y) > precedence(x & y)
    assert precedence(~y) > precedence(x*y)
    assert precedence(~y) > precedence(x - y)
    # double checks
    assert precedence(x & y) == PRECEDENCE["And"]
    assert precedence(x | y) == PRECEDENCE["Or"]
    assert precedence(~y) == PRECEDENCE["Not"]
