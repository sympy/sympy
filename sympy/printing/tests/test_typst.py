from sympy.concrete.summations import Sum
from sympy.core.function import Function
from sympy.core.mul import Mul
from sympy.core.numbers import (I, Rational, oo)
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import log
from sympy.printing.typst import typst
from sympy.series.limits import Limit

x, y = symbols('x, y')

def test_printmethod():
    class R(Abs):
        def _typst(self, printer):
            return "foo(%s)" % printer._print(self.args[0])
    assert typst(R(x)) == r"foo(x)"

    class R(Abs):
        def _typst(self, printer):
            return "foo"
    assert typst(R(x)) == r"foo"


def test_typst_basic():
    assert typst(1 + x) == r"x + 1"
    assert typst(x**2) == r"x^(2)"
    assert typst(x**(1 + x)) == r"x^(x + 1)"
    assert typst(x**3 + x + 1 + x**2) == r"x^(3) + x^(2) + x + 1"

    assert typst(2*x*y) == r"2 x y"
    assert typst(2*x*y, mul_symbol='dot') == r"2 dot x dot y"
    assert typst(3*x**2*y, mul_symbol='#h(1cm)') == r"3#h(1cm)x^(2)#h(1cm)y"
    assert typst(1.5*3**x, mul_symbol='#h(1cm)') == r"1.5#h(1cm)3^(x)"

    assert typst(x**S.Half**5) == r"root(32, x)"
    assert typst(Mul(S.Half, x**2, -5, evaluate=False)) == r"1/2 x^(2) (-5)"
    assert typst(Mul(S.Half, x**2, 5, evaluate=False)) == r"1/2 x^(2) 5"
    assert typst(Mul(-5, -5, evaluate=False)) == r"(-5) (-5)"
    assert typst(Mul(5, -5, evaluate=False)) == r"5 (-5)"
    assert typst(Mul(S.Half, -5, S.Half, evaluate=False)) == r"1/2 (-5) 1/2"
    assert typst(Mul(5, I, 5, evaluate=False)) == r"5 i 5"
    assert typst(Mul(5, I, -5, evaluate=False)) == r"5 i (-5)"
    assert typst(Mul(Pow(x, 2), S.Half*x + 1)) == r"x^(2) (x/2 + 1)"
    assert typst(Mul(Pow(x, 3), Rational(2, 3)*x + 1)) == r"x^(3) ((2 x)/3 + 1)"
    assert typst(Mul(Pow(x, 11), 2*x + 1)) == r"x^(11) (2 x + 1)"


def test_typst_limits():
    assert typst(Limit(x, x, oo)) == r"lim_(x -> infinity) x"

    # issue 8175
    f = Function('f')
    assert typst(Limit(f(x), x, 0)) == r"lim_(x -> 0^+) f(x)"
    assert typst(Limit(f(x), x, 0, "-")) == \
        r"lim_(x -> 0^-) f(x)"

    # issue #10806
    assert typst(Limit(f(x), x, 0)**2) == \
        r"(lim_(x -> 0^+) f(x))^(2)"
    # bi-directional limit
    assert typst(Limit(f(x), x, 0, dir='+-')) == \
        r"lim_(x -> 0) f(x)"


def test_typst_log():
    assert typst(log(x)) == r"log(x)"
    assert typst(log(x), ln_notation=True) == r"ln(x)"
    assert typst(log(x) + log(y)) == \
        r"log(x) + log(y)"
    assert typst(log(x) + log(y), ln_notation=True) == \
        r"ln(x) + ln(y)"
    assert typst(pow(log(x), x)) == r"log(x)^(x)"
    assert typst(pow(log(x), x), ln_notation=True) == \
        r"ln(x)^(x)"


def test_typst_sum():
    assert typst(Sum(x*y**2, (x, -2, 2), (y, -5, 5))) == \
        r"sum_(-2 <= x <= 2 \ -5 <= y <= 5) x y^(2)"
    assert typst(Sum(x**2, (x, -2, 2))) == \
        r"sum_(x=-2)^(2) x^(2)"
    assert typst(Sum(x**2 + y, (x, -2, 2))) == \
        r"sum_(x=-2)^(2) (x^(2) + y)"
    assert typst(Sum(x**2 + y, (x, -2, 2))**2) == \
        r"(sum_(x=-2)^(2) (x^(2) + y))^(2)"
