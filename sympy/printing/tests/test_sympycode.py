# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)
from sympy import (Symbol, symbols, oo, limit, Rational, Integral, Derivative,
    log, exp, sqrt, pi, Function, sin, Eq, Ge, Le, Gt, Lt, Ne, Abs, conjugate,
    I, Matrix)
from sympy.utilities.pytest import raises, XFAIL
from ..sympycode import sympy_code


x, y = symbols('x,y')
th = Symbol('theta')
ph = Symbol('phi')


def test_sympy_code_basic():
    # Simple numbers/symbols
    assert sympy_code(-Rational(1)/2) == "e = Rational(-1, 2)"
    assert sympy_code(-Rational(13)/22) == "e = Rational(-13, 22)"
    assert sympy_code(oo) == "e = oo"

    # Powers
    assert sympy_code((x**2)) == "x = Symbol(\'x\')\ne = x**2"
    assert sympy_code(1/x) == "x = Symbol('x')\ne = 1/x"
    assert sympy_code(y*x**-2) == "y = Symbol('y')\nx = Symbol('x')\ne = y/x**2"
    assert sympy_code(
        x**Rational(-5, 2)) == "x = Symbol('x')\ne = x**Rational(-5, 2)"

    # Sums of terms
    assert sympy_code((x**2 + x + 1)) in [
        "x = Symbol('x')\ne = 1 + x + x**2",
        "x = Symbol('x')\ne = x + x**2 + 1",
        "x = Symbol('x')\ne = x**2 + x + 1", ]
    assert sympy_code(1 - x) in [
        "x = Symbol('x')\ne = 1 - x",
        "x = Symbol('x')\ne = -x + 1"]
    assert sympy_code(1 - 2*x) in [
        "x = Symbol('x')\ne = 1 - 2*x",
        "x = Symbol('x')\ne = -2*x + 1"]
    assert sympy_code(1 - Rational(3, 2)*y/x) in [
        "y = Symbol('y')\nx = Symbol('x')\ne = 1 - 3/2*y/x",
        "y = Symbol('y')\nx = Symbol('x')\ne = -3/2*y/x + 1",
        "y = Symbol('y')\nx = Symbol('x')\ne = 1 - 3*y/(2*x)"]

    # Multiplication
    assert sympy_code(x/y) == "x = Symbol('x')\ny = Symbol('y')\ne = x/y"
    assert sympy_code(-x/y) == "x = Symbol('x')\ny = Symbol('y')\ne = -x/y"
    assert sympy_code((x + 2)/y) in [
        "y = Symbol('y')\nx = Symbol('x')\ne = 1/y*(2 + x)",
        "y = Symbol('y')\nx = Symbol('x')\ne = 1/y*(x + 2)",
        "x = Symbol('x')\ny = Symbol('y')\ne = 1/y*(2 + x)",
        "x = Symbol('x')\ny = Symbol('y')\ne = (2 + x)/y",
        "x = Symbol('x')\ny = Symbol('y')\ne = (x + 2)/y"]
    assert sympy_code((1 + x)*y) in [
        "y = Symbol('y')\nx = Symbol('x')\ne = y*(1 + x)",
        "y = Symbol('y')\nx = Symbol('x')\ne = y*(x + 1)", ]

    # Check for proper placement of negative sign
    assert sympy_code(-5*x/(x + 10)) == "x = Symbol('x')\ne = -5*x/(x + 10)"
    assert sympy_code(1 - Rational(3, 2)*(x + 1)) in [
        "x = Symbol('x')\ne = Rational(-3, 2)*x + Rational(-1, 2)",
        "x = Symbol('x')\ne = -3*x/2 + Rational(-1, 2)",
        "x = Symbol('x')\ne = -3*x/2 + Rational(-1, 2)"
    ]


def test_sympy_code_keyword_symbol_name_escaping():
    # Check for escaping of keywords
    assert sympy_code(
        5*Symbol("lambda")) == "lambda_ = Symbol('lambda')\ne = 5*lambda_"
    assert (sympy_code(5*Symbol("lambda") + 7*Symbol("lambda_")) ==
            "lambda__ = Symbol('lambda')\nlambda_ = Symbol('lambda_')\ne = 7*lambda_ + 5*lambda__")
    assert (sympy_code(5*Symbol("for") + Function("for_")(8)) ==
            "for__ = Symbol('for')\nfor_ = Function('for_')\ne = 5*for__ + for_(8)")

def test_sympy_code_keyword_function_name_escaping():
    assert sympy_code(
        5*Function("for")(8)) == "for_ = Function('for')\ne = 5*for_(8)"


def test_sympy_code_relational():
    assert sympy_code(Eq(x, y)) == "e = Eq(x, y)"
    assert sympy_code(Ge(x, y)) == "x = Symbol('x')\ny = Symbol('y')\ne = x >= y"
    assert sympy_code(Le(x, y)) == "x = Symbol('x')\ny = Symbol('y')\ne = x <= y"
    assert sympy_code(Gt(x, y)) == "x = Symbol('x')\ny = Symbol('y')\ne = x > y"
    assert sympy_code(Lt(x, y)) == "x = Symbol('x')\ny = Symbol('y')\ne = x < y"
    assert sympy_code(Ne(x/(y + 1), y**2)) in ["e = Ne(x/(1 + y), y**2)", "e = Ne(x/(y + 1), y**2)"]


def test_sympy_code_functions():
    # Simple
    assert sympy_code((2*x + exp(x))) in "x = Symbol('x')\ne = 2*x + exp(x)"
    assert sympy_code(sqrt(2)) == 'e = sqrt(2)'
    assert sympy_code(2**Rational(1, 3)) == 'e = 2**Rational(1, 3)'
    assert sympy_code(sqrt(2 + pi)) == 'e = sqrt(2 + pi)'
    assert sympy_code((2 + pi)**Rational(1, 3)) == 'e = (2 + pi)**Rational(1, 3)'
    assert sympy_code(2**Rational(1, 4)) == 'e = 2**Rational(1, 4)'
    assert sympy_code(Abs(x)) == "x = Symbol('x')\ne = Abs(x)"
    assert sympy_code(
        Abs(x/(x**2 + 1))) in ["x = Symbol('x')\ne = Abs(x/(1 + x**2))",
            "x = Symbol('x')\ne = Abs(x/(x**2 + 1))"]

    # Univariate/Multivariate functions
    f = Function('f')
    assert sympy_code(f(x)) == "x = Symbol('x')\nf = Function('f')\ne = f(x)"
    assert sympy_code(f(x, y)) == "x = Symbol('x')\ny = Symbol('y')\nf = Function('f')\ne = f(x, y)"
    assert sympy_code(f(x/(y + 1), y)) in [
        "x = Symbol('x')\ny = Symbol('y')\nf = Function('f')\ne = f(x/(1 + y), y)",
        "x = Symbol('x')\ny = Symbol('y')\nf = Function('f')\ne = f(x/(y + 1), y)"]

    # Nesting of square roots
    assert sympy_code(sqrt((sqrt(x + 1)) + 1)) in [
        "x = Symbol('x')\ne = sqrt(1 + sqrt(1 + x))",
        "x = Symbol('x')\ne = sqrt(sqrt(x + 1) + 1)"]

    # Nesting of powers
    assert sympy_code((((x + 1)**Rational(1, 3)) + 1)**Rational(1, 3)) in [
        "x = Symbol('x')\ne = (1 + (1 + x)**Rational(1, 3))**Rational(1, 3)",
        "x = Symbol('x')\ne = ((x + 1)**Rational(1, 3) + 1)**Rational(1, 3)"]

    # Function powers
    assert sympy_code(sin(x)**2) == "x = Symbol('x')\ne = sin(x)**2"


@XFAIL
def test_sympy_code_functions_conjugates():
    a, b = map(Symbol, 'ab')
    assert sympy_code( conjugate(a + b*I) ) == '_     _\na - I*b'
    assert sympy_code( conjugate(exp(a + b*I)) ) == ' _     _\n a - I*b\ne       '


def test_sympy_code_derivatives():
    # Simple
    f_1 = Derivative(log(x), x, evaluate=False)
    assert sympy_code(f_1) == "x = Symbol('x')\ne = Derivative(log(x), x)"

    f_2 = Derivative(log(x), x, evaluate=False) + x
    assert sympy_code(f_2) == "x = Symbol('x')\ne = x + Derivative(log(x), x)"

    # Multiple symbols
    f_3 = Derivative(log(x) + x**2, x, y, evaluate=False)
    assert sympy_code(f_3) == \
        "x = Symbol('x')\ny = Symbol('y')\ne = Derivative(x**2 + log(x), x, y)"

    f_4 = Derivative(2*x*y, y, x, evaluate=False) + x**2
    assert sympy_code(f_4) in [
        "x = Symbol('x')\ny = Symbol('y')\ne = x**2 + Derivative(2*x*y, y, x)",
        "x = Symbol('x')\ny = Symbol('y')\ne = Derivative(2*x*y, y, x) + x**2"]


def test_sympy_code_integrals():
    # Simple
    f_1 = Integral(log(x), x)
    assert sympy_code(f_1) == "x = Symbol('x')\ne = Integral(log(x), x)"

    f_2 = Integral(x**2, x)
    assert sympy_code(f_2) == "x = Symbol('x')\ne = Integral(x**2, x)"

    # Double nesting of pow
    f_3 = Integral(x**(2**x), x)
    assert sympy_code(f_3) == "x = Symbol('x')\ne = Integral(x**(2**x), x)"

    # Definite integrals
    f_4 = Integral(x**2, (x, 1, 2))
    assert sympy_code(f_4) == "x = Symbol('x')\ne = Integral(x**2, (x, 1, 2))"

    f_5 = Integral(x**2, (x, Rational(1, 2), 10))
    assert sympy_code(
        f_5) == "x = Symbol('x')\ne = Integral(x**2, (x, Rational(1, 2), 10))"

    # Nested integrals
    f_6 = Integral(x**2*y**2, x, y)
    assert sympy_code(f_6) == "x = Symbol('x')\ny = Symbol('y')\ne = Integral(x**2*y**2, x, y)"


def test_sympy_code_matrix():
    p = sympy_code(Matrix([[x**2+1, 1], [y, x+y]]))
    s = "x = Symbol('x')\ny = Symbol('y')\ne = MutableDenseMatrix([[x**2 + 1, 1], [y, x + y]])"
    assert p == s

def test_sympy_code_limits():
    assert sympy_code(limit(x, x, oo)) == 'e = oo'
    assert sympy_code(limit(x**2, x, 0)) == 'e = 0'


def test_settings():
    raises(TypeError, lambda: sympy_code(x, method="garbage"))
