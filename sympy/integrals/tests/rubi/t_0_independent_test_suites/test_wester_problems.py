"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Wester Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, m = symbols('a b m')

def test_integrate_0_Independent_test_suites_Wester_Problems_1():
    f = (3*x - 5)**2/(2*x - 1)**(sympy.S(7)/2)
    F = -9/(4*sqrt(2*x - 1)) + 7/(2*(2*x - 1)**(sympy.S(3)/2)) - 49/(20*(2*x - 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Wester_Problems_2():
    f = 1/(2*exp(m*x) - 5*exp(-m*x))
    F = -sqrt(10)*atanh(sqrt(10)*exp(m*x)/5)/(10*m)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Wester_Problems_3():
    f = 1/(a + b*cos(x))
    F = 2*atan(sqrt(a - b)*tan(x/2)/sqrt(a + b))/(sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Wester_Problems_4():
    f = 1/(4*sin(x) + 3*cos(x) + 3)
    F = log(4*tan(x/2) + 3)/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Wester_Problems_5():
    f = 1/(4*sin(x) + 3*cos(x) + 4)
    F = -log(3*cot(x/2 + pi/4) + 4)/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Wester_Problems_6():
    f = 1/(4*sin(x) + 3*cos(x) + 6)
    F = sqrt(11)*x/11 + 2*sqrt(11)*atan((-3*sin(x) + 4*cos(x))/(4*sin(x) + 3*cos(x) + sqrt(11) + 6))/11
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Wester_Problems_7():
    f = log((-a**2 + x**2)**2)/2
    F = 2*a*atanh(x/a) + x*log((-a**2 + x**2)**2)/2 - 2*x
    assert integrate(f, x) == F

