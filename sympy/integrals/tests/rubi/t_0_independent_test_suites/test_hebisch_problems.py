"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Hebisch Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

def test_integrate_0_Independent_test_suites_Hebisch_Problems_1():
    f = (x**6 - x**5 + x**4 - x**3 + 1)*exp(x)
    F = x**6*exp(x) - 7*x**5*exp(x) + 36*x**4*exp(x) - 145*x**3*exp(x) + 435*x**2*exp(x) - 870*x*exp(x) + 871*exp(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hebisch_Problems_2():
    f = (2 - x**2)*exp(x/(x**2 + 2))/(x**3 + 2*x)
    F = sympy.Function('ExpIntegralEi')((x * ((Integer(2) + (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hebisch_Problems_3():
    f = (2*x**4 - x**3 + 3*x**2 + 2*x + 2)*exp(x/(x**2 + 2))/(x**3 + 2*x)
    F = ((sympy.E)**((x * ((Integer(2) + (x)**(Integer(2))))**(Integer(-1)))) * (Integer(2) + (x)**(Integer(2)))) + sympy.Function('ExpIntegralEi')((x * ((Integer(2) + (x)**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hebisch_Problems_4():
    f = (exp(x) + 1)*exp(x + exp(x))/(x + exp(x))
    F = sympy.Function('ExpIntegralEi')(((sympy.E)**(x) + x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hebisch_Problems_5():
    f = (x**3 - x**2 - 3*x + 1)*exp(1/(x**2 - 1))/(x**3 - x**2 - x + 1)
    F = (x + 1)*exp(1/(x**2 - 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hebisch_Problems_6():
    f = (log(x)**2 - 1)*exp(1 + 1/log(x))/log(x)**2
    F = x*exp(1 + 1/log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Hebisch_Problems_7():
    f = ((x + 1)*log(x)**2 - 1)*exp(x + 1/log(x))/log(x)**2
    F = x*exp(x + 1/log(x))
    assert integrate(f, x) == F

