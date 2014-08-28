from sympy.core import (S, pi, oo, symbols, Function,
                        Rational, Integer, Tuple)
from sympy.core import EulerGamma, GoldenRatio, Catalan, Lambda
from sympy.functions import Piecewise, sqrt, Abs, ceiling, exp, sin, cos
from sympy.utilities.pytest import raises
from sympy.utilities.lambdify import implemented_function
from sympy.matrices import eye, Matrix
from sympy.integrals import Integral
from sympy.concrete import Sum
from sympy.utilities.pytest import XFAIL

from sympy import octave_code as mcode

x, y, z = symbols('x,y,z')


def test_Integer():
    assert mcode(Integer(67)) == "67"
    assert mcode(Integer(-1)) == "-1"


def test_Rational():
    assert mcode(Rational(3, 7)) == "3/7"
    assert mcode(Rational(18, 9)) == "2"
    assert mcode(Rational(3, -7)) == "-3/7"
    assert mcode(Rational(-3, -7)) == "3/7"
    assert mcode(x + Rational(3, 7)) == "x + 3/7"
    assert mcode(Rational(3, 7)*x) == "(3/7)*x"


def test_Function():
    assert mcode(sin(x) ** cos(x)) == "sin(x)^cos(x)"
    assert mcode(abs(x)) == "abs(x)"
    assert mcode(ceiling(x)) == "ceil(x)"


def test_Pow():
    assert mcode(x**3) == "x^3"
    assert mcode(x**(y**3)) == "x^(y^3)"
    g = implemented_function('g', Lambda(x, 2*x))
    assert mcode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "(3.5*2*x)^(-x + y^x)/(x^2 + y)"
    assert mcode(x**Rational(2, 3)) == 'x^(2/3)'
    assert mcode(x**-1) == '1/x'
    assert mcode(x**S.Half) == 'sqrt(x)'
    assert mcode(sqrt(x)) == "sqrt(x)"
    assert mcode(x**-S.Half) == '1/sqrt(x)'

def test_Mul():
    assert mcode(x*y*z) == "x*y*z"

def test_imag():
    I = S('I')
    assert mcode(I) == "1i"
    assert mcode(5*I) == "5i"
    assert mcode((S(3)/2)*I) == "(3/2)*1i"
    assert mcode(3+4*I) == "3 + 4i"

def test_constants():
    assert mcode(pi) == "pi"
    assert mcode(oo) == "inf"
    assert mcode(-oo) == "-inf"
    assert mcode(S.NegativeInfinity) == "-inf"
    assert mcode(S.NaN) == "NaN"
    assert mcode(S.Exp1) == "exp(1)"
    assert mcode(exp(1)) == "exp(1)"

def test_ccode_constants_other():
    assert mcode(2*GoldenRatio) == "2*(1+sqrt(5))/2"
    assert mcode(2*Catalan) == "2*0.915965594177219011"
    assert mcode(2*EulerGamma) == "2*0.577215664901532866"

def test_ccode_boolean():
    assert mcode(x & y) == "x && y"
    assert mcode(x | y) == "x || y"
    assert mcode(~x) == "~x"
    assert mcode(x & y & z) == "x && y && z"
    assert mcode(x | y | z) == "x || y || z"
    assert mcode((x & y) | z) == "z || x && y"
    assert mcode((x | y) & z) == "z && (x || y)"

def test_Matrices():
    assert mcode(Matrix(1, 1, [10])) == "[10]"
    A = Matrix([[1, sin(x/2), abs(x)],
                [0, 1, EulerGamma],
                [0, exp(1), ceiling(x)]]);
    assert mcode(A) == (
        "[[1, sin((1/2)*x),               abs(x)]; ...\n"
        "[0,            1, 0.577215664901532866]; ...\n"
        "[0,       exp(1),              ceil(x)]]")
    # row and columns
    assert mcode(A[:,0]) == "[[1];  [0];  [0]]"
    assert mcode(A[0,:]) == "[1, sin((1/2)*x), abs(x)]"
    # empty matrices
    assert mcode(Matrix(0, 0, [])) == '[]'
    assert mcode(Matrix(0, 3, [])) == 'zeros(0,3)'

def test_containers():
    assert mcode([1, 2, 3, [4, 5, [6, 7]], 8, [9, 10], 11]) == \
        "{1, 2, 3, {4, 5, {6, 7}}, 8, {9, 10}, 11}"
    assert mcode((1, 2, (3, 4))) == "{1, 2, {3, 4}}"
    assert mcode([1]) == "{1}"
    assert mcode((1,)) == "{1}"
    assert mcode(Tuple(*[1, 2, 3])) == "{1, 2, 3}"
    assert mcode((1, x*y, (3, x**2))) == "{1, x*y, {3, x^2}}"

def test_octave_piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    assert mcode(expr) == ("(x < 1).*(x) + (~(x < 1)).*(x^2)")
    # FIXME: fix up if we get indenting
    assert mcode(expr, assign_to="c") == (
            "if (x < 1)\n"
            "c = x;\n"
            "else\n"
            "c = x^2;\n"
            "end")
    expr = Piecewise((x**2, x < 1), (x**3, x < 2), (x**4, x < 3), (x**5, True))
    assert mcode(expr) == (
            "(x < 1).*(x^2) + (~(x < 1)).*( ...\n"
            "(x < 2).*(x^3) + (~(x < 2)).*( ...\n"
            "(x < 3).*(x^4) + (~(x < 3)).*(x^5)))")
    assert mcode(expr, assign_to='c') == (
            "if (x < 1)\n"
            "c = x^2;\n"
            "elseif (x < 2)\n"
            "c = x^3;\n"
            "elseif (x < 3)\n"
            "c = x^4;\n"
            "else\n"
            "c = x^5;\n"
            "end")
    # Check that Piecewise without a True (default) condition error
    expr = Piecewise((x, x < 1), (x**2, x > 1), (sin(x), x > 0))
    raises(ValueError, lambda: mcode(expr))

@XFAIL
def test_octave_piecewise_times_const():
    # FIXME: self.parenthesize something! but where?
    expr = Piecewise((x, x < 1), (x**2, True))
    assert mcode(2*expr) == ("2*((x < 1).*(x) + (~(x < 1)).*(x^2))")
