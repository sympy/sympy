# -*- coding: utf-8 -*-

from sympy import Symbol, symbols, oo, limit, Rational, Integral, Derivative, \
    log, exp, sqrt, pi, Function, sin, Eq, Le, Gt, Ne, raises

from sympy.printing.python import python

x, y = symbols('xy')
th  = Symbol('theta')
ph  = Symbol('phi')

def test_python_basic():
    # Simple numbers/symbols
    assert python(-Rational(1)/2) == "e = Rational(-1, 2)"
    assert python(-Rational(13)/22) == "e = Rational(-13, 22)"
    assert python(oo) == "e = oo"

    # Powers
    assert python((x**2)) == "x = Symbol(\'x\')\ne = x**2"
    assert python(1/x) == "x = Symbol('x')\ne = 1/x"
    assert python(y*x**-2) == "y = Symbol('y')\nx = Symbol('x')\ne = y/x**2"
    assert python(x**Rational(-5, 2)) == "x = Symbol('x')\ne = x**(Rational(-5, 2))"

    # Sums of terms
    assert python((x**2 + x + 1)) in [
            "x = Symbol('x')\ne = 1 + x + x**2",
            "x = Symbol('x')\ne = x + x**2 + 1",
            "x = Symbol('x')\ne = x**2 + x + 1",]
    assert python(1-x) in [
            "x = Symbol('x')\ne = 1 - x",
            "x = Symbol('x')\ne = -x + 1"]
    assert python(1-2*x) in [
            "x = Symbol('x')\ne = 1 - 2*x",
            "x = Symbol('x')\ne = -2*x + 1"]
    assert python(1-Rational(3,2)*y/x) in [
            "y = Symbol('y')\nx = Symbol('x')\ne = 1 - 3/2*y/x",
            "y = Symbol('y')\nx = Symbol('x')\ne = -3/2*y/x + 1",
            "y = Symbol('y')\nx = Symbol('x')\ne = 1 - 3*y/(2*x)"]

    # Multiplication
    assert python(x/y) == "x = Symbol('x')\ny = Symbol('y')\ne = x/y"
    assert python(-x/y) == "x = Symbol('x')\ny = Symbol('y')\ne = -x/y"
    assert python((x+2)/y) in [
            "y = Symbol('y')\nx = Symbol('x')\ne = 1/y*(2 + x)",
            "y = Symbol('y')\nx = Symbol('x')\ne = 1/y*(x + 2)",
            "x = Symbol('x')\ny = Symbol('y')\ne = 1/y*(2 + x)",
            "x = Symbol('x')\ny = Symbol('y')\ne = (2 + x)/y"]
    assert python((1+x)*y) in [
            "y = Symbol('y')\nx = Symbol('x')\ne = y*(1 + x)",
            "y = Symbol('y')\nx = Symbol('x')\ne = y*(x + 1)",]

    # Check for proper placement of negative sign
    assert python(-5*x/(x+10)) == "x = Symbol('x')\ne = -5*x/(10 + x)"
    assert python(1 - Rational(3,2)*(x+1)) in [
            "x = Symbol('x')\ne = Rational(-1, 2) + Rational(-3, 2)*x",
            "x = Symbol('x')\ne = Rational(-1, 2) - 3*x/2",
            "x = Symbol('x')\ne = Rational(-1, 2) - 3*x/2"
            ]

def test_python_relational():
    assert python(Eq(x, y)) == "x = Symbol('x')\ny = Symbol('y')\ne = x == y"
    assert python(Le(x, y)) == "x = Symbol('x')\ny = Symbol('y')\ne = x <= y"
    assert python(Gt(x, y)) == "y = Symbol('y')\nx = Symbol('x')\ne = y < x"
    assert python(Ne(x/(y+1), y**2)) in [
            "x = Symbol('x')\ny = Symbol('y')\ne = x/(1 + y) != y**2",
            "x = Symbol('x')\ny = Symbol('y')\ne = x/(y + 1) != y**2"]

def test_python_functions():
    # Simple
    assert python((2*x + exp(x))) in "x = Symbol('x')\ne = 2*x + exp(x)"
    assert python(sqrt(2)) == 'e = 2**(Half(1, 2))'
    assert python(sqrt(2+pi)) == 'e = (2 + pi)**(Half(1, 2))'
    assert python(abs(x)) == "x = Symbol('x')\ne = abs(x)"
    assert python(abs(x/(x**2+1))) in ["x = Symbol('x')\ne = abs(x/(1 + x**2))",
            "x = Symbol('x')\ne = abs(x/(x**2 + 1))"]

    # Univariate/Multivariate functions
    f = Function('f')
    assert python(f(x)) == "x = Symbol('x')\nf = Function('f')\ne = f(x)"
    assert python(f(x, y)) == "x = Symbol('x')\ny = Symbol('y')\nf = Function('f')\ne = f(x, y)"
    assert python(f(x/(y+1), y)) in [
        "x = Symbol('x')\ny = Symbol('y')\nf = Function('f')\ne = f(x/(1 + y), y)",
        "x = Symbol('x')\ny = Symbol('y')\nf = Function('f')\ne = f(x/(y + 1), y)"]

    # Nesting of square roots
    assert python(sqrt((sqrt(x+1))+1)) in [
            "x = Symbol('x')\ne = (1 + (1 + x)**(Half(1, 2)))**(Half(1, 2))",
            "x = Symbol('x')\ne = ((x + 1)**(Half(1, 2)) + 1)**(Half(1, 2))"]
    # Function powers
    assert python(sin(x)**2) == "x = Symbol('x')\ne = sin(x)**2"

    # Conjugates
    a, b = map(Symbol, 'ab')
    #assert python( conjugate(a+b*I) ) == '_     _\na - I*b'
    #assert python( conjugate(exp(a+b*I)) ) == ' _     _\n a - I*b\ne       '

def test_python_derivatives():
    # Simple
    f_1 = Derivative(log(x), x, evaluate=False)
    assert python(f_1) == "x = Symbol('x')\ne = D(log(x), x)"

    f_2 = Derivative(log(x), x, evaluate=False) + x
    assert python(f_2) == "x = Symbol('x')\ne = x + D(log(x), x)"

    # Multiple symbols
    f_3 = Derivative(log(x) + x**2, x, y, evaluate=False)
    #assert python(f_3) ==

    f_4 = Derivative(2*x*y, y, x, evaluate=False) + x**2
    assert python(f_4) in [
            "x = Symbol('x')\ny = Symbol('y')\ne = x**2 + D(2*x*y, y, x)",
            "x = Symbol('x')\ny = Symbol('y')\ne = D(2*x*y, y, x) + x**2"]

def test_python_integrals():
    # Simple
    f_1 = Integral(log(x), x)
    assert python(f_1) == "x = Symbol('x')\ne = Integral(log(x), x)"

    f_2 = Integral(x**2, x)
    assert python(f_2) == "x = Symbol('x')\ne = Integral(x**2, x)"

    # Double nesting of pow
    f_3 = Integral(x**(2**x), x)
    assert python(f_3) == "x = Symbol('x')\ne = Integral(x**(2**x), x)"

    # Definite integrals
    f_4 = Integral(x**2, (x,1,2))
    assert python(f_4) == "x = Symbol('x')\ne = Integral(x**2, (x, 1, 2))"

    f_5 = Integral(x**2, (x,Rational(1,2),10))
    assert python(f_5) == "x = Symbol('x')\ne = Integral(x**2, (x, Half(1, 2), 10))"

    # Nested integrals
    f_6 = Integral(x**2*y**2, x,y)
    assert python(f_6) == "x = Symbol('x')\ny = Symbol('y')\ne = Integral(x**2*y**2, x, y)"


# Not implemented yet
#def test_python_matrix():
#    p = python(Matrix([[x**2+1, 1], [y, x+y]]))
#    s = ''
#    assert p == s

def test_python_limits():
    assert python(limit(x, x, oo)) == 'e = oo'
    assert python(limit(x**2, x, 0)) == 'e = 0'

def test_settings():
    raises(TypeError, 'python(x, method="garbage")')
