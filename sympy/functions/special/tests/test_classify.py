from sympy import Rational, Symbol, gamma, sqrt, exp, cos, sin, acosh, E, pi
from sympy.functions.special.classify import \
    LINEAR, POLYNOMIAL, RATIONAL, EXP, LOG, TRIGONOMETRIC, NONELEMENTARY, ROOT,\
    INVHYPERBOLIC
from sympy.functions.special.classify import classify

x = Symbol('x')

def test_nonfunction():
    assert classify(pi, x) == []

def test_linear():
    assert classify(3*x + 1, x) == [LINEAR]
    assert classify((4*x+5)*(1+pi), x) == [LINEAR]
    assert classify((E*(pi + 8*x)+5)*(1+pi), x) == [LINEAR]

def test_polynomal():
    assert classify(3*x + x**2, x) == [POLYNOMIAL]
    assert classify((4*x+5)*(2*x+1), x) == [POLYNOMIAL]

def test_rational():
    assert classify(1/x, x) == [RATIONAL]
    assert classify(1/(2*x+1), x) == [RATIONAL]
    assert classify(1/(2*x+x**2), x) == [RATIONAL]
    assert classify((1+x)/(1-x), x) == [RATIONAL]
    assert classify((4+x)/(2*x+x**2), x) == [RATIONAL]
    assert classify(x**2 + (4+x)/(2*x+x**2), x) == [RATIONAL]
    assert classify(x**2*(4+x)/(2*x+x**2), x) == [RATIONAL]

def test_mixed():
    assert classify(sqrt(3*x+1), x) == [ROOT, LINEAR]
    assert classify(x**Rational(3,2), x) == [ROOT, LINEAR]
    assert classify(x**pi, x) == [EXP, LINEAR, LOG, LINEAR]
    assert classify(2*x + exp(3*x+4) + exp(2*x+1), x) == [LINEAR, EXP, LINEAR]
    assert classify(2*x + (exp(3*x) + exp(x+1))/(5*x+1), x) == [RATIONAL, EXP, LINEAR]
    assert classify(2*x + (cos(x+4)**4 + sin(2*x+1)**2), x) == [POLYNOMIAL, TRIGONOMETRIC, LINEAR]
    assert classify(2*x + (cos(x+4)**4 + sin(2*x+1)**2) / (5*x + 1), x) == [RATIONAL, TRIGONOMETRIC, LINEAR]
    assert classify(2**x, x) == [EXP, LINEAR]
    assert classify(x**x, x) == [EXP, LOG, LINEAR]
    assert classify((2*x+1)**x, x) == [EXP, LOG, LINEAR]
    assert classify((acosh(x) + 1)**x, x) == [EXP, LOG, LINEAR, INVHYPERBOLIC, LINEAR]
    assert classify(exp(exp(x)), x) == [EXP, EXP, LINEAR]
    assert classify(exp(exp(1/x)), x) == [EXP, EXP, RATIONAL]
    assert classify(exp(exp(x)) + exp(exp(1/x)), x) == [LINEAR, EXP, EXP, RATIONAL]
    assert classify(exp(gamma(x)**2), x) == [EXP, POLYNOMIAL, NONELEMENTARY, LINEAR]
