"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Jeffrey Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

p, q, r = symbols('p q r')

def test_integrate_0_Independent_test_suites_Jeffrey_Problems_1():
    f = 3/(5 - 4*cos(x))
    F = x + 2*atan(sin(x)/(2 - cos(x)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_2():
    f = (2*sin(x) + cos(x) + 1)/(-2*sin(x)*cos(x) + 2*sin(x) + cos(x)**2 + 3)
    F = -atan((-sin(x) + 2*cos(x))/(sin(x) + 2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_3():
    f = (5*sin(x) + cos(x) + 2)/(-2*sin(x)**2 + sin(x)*cos(x) - 2*sin(x) + 4*cos(x))
    F = -log(sin(x) - 3*cos(x) + 1) + log(sin(x) + cos(x) + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_4():
    f = (2*sin(x) + 7*cos(x) + 3)/(-sin(x)*cos(x) - 5*sin(x) + 3*cos(x)**2 + 4*cos(x) + 1)
    F = -log(-2*sin(x) + cos(x) + 1) + log(sin(x) + cos(x) + 3)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_5():
    f = (5*cos(x)**2 + 4*cos(x) - 1)/(4*cos(x)**3 - 3*cos(x)**2 - 4*cos(x) - 1)
    F = x - 2*atan((7*sin(x)*cos(x) + 3*sin(x))/(5*cos(x)**2 + 2*cos(x) + 1)) - 2*atan(sin(x)/(cos(x) + 3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_6():
    f = (7*cos(x)**2 + 2*cos(x) - 5)/(4*cos(x)**3 - 9*cos(x)**2 + 2*cos(x) - 1)
    F = x - 2*atan(2*sin(x)*cos(x)/(2*cos(x)**2 - cos(x) + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_7():
    f = 3/(4*sin(x) + 5)
    F = x + 2*atan(cos(x)/(sin(x) + 2))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_8():
    f = 2/(cos(x)**2 + 1)
    F = sqrt(2)*x - sqrt(2)*atan(sin(x)*cos(x)/(cos(x)**2 + 1 + sqrt(2)))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Jeffrey_Problems_9():
    f = 1/(p + q*cos(x) + r*sin(x))
    F = 2*atan((r + (p - q)*tan(x/2))/sqrt(p**2 - q**2 - r**2))/sqrt(p**2 - q**2 - r**2)
    assert integrate(f, x) == F

