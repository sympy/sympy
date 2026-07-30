"""Generated from MathematicaSyntaxTestSuite.

Source: 0 Independent test suites/Bronstein Problems.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

def test_integrate_0_Independent_test_suites_Bronstein_Problems_1():
    f = sqrt(x**8 + 1)*(2*x**8 + 1)/(x**17 + 2*x**9 + x)
    F = -atanh(sqrt(x**8 + 1))/4 - 1/(4*sqrt(x**8 + 1))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_2():
    f = 1/(x**2 + 1)
    F = atan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_3():
    f = 1/(x*sqrt(x**8 + 1))
    F = -atanh(sqrt(x**8 + 1))/4
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_4():
    f = x/sqrt(1 - x**3)
    F = 2*sqrt(1 - x**3)/(-x + 1 + sqrt(3)) - 3**(sympy.S(1)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*sqrt(2 - sqrt(3))*elliptic_e(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3)) + 2*sqrt(2)*3**(sympy.S(3)/4)*sqrt((x**2 + x + 1)/(-x + 1 + sqrt(3))**2)*(1 - x)*elliptic_f(asin((-x - sqrt(3) + 1)/(-x + 1 + sqrt(3))), -7 - 4*sqrt(3))/(3*sqrt((1 - x)/(-x + 1 + sqrt(3))**2)*sqrt(1 - x**3))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_5():
    f = 1/(x*sqrt(1 - x**3))
    F = -2*atanh(sqrt(1 - x**3))/3
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_6():
    f = x/sqrt(x**4 + 10*x**2 - 96*x - 71)
    F = log(x**8 + 20*x**6 - 128*x**5 + 54*x**4 - 1408*x**3 + 3124*x**2 + sqrt(x**4 + 10*x**2 - 96*x - 71)*(x**6 + 15*x**4 - 80*x**3 + 27*x**2 - 528*x + 781) + 10001)/8
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_7():
    f = (x - tan(x))/tan(x)**2
    F = -x**2/2 - x*cot(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_8():
    f = x*tan(x) + tan(x)**2 + 1
    F = ((sympy.I * (x)**(Integer(2))) * (Integer(2))**(Integer(-1))) + (Integer(-1) * (x * sympy.log((Integer(1) + (sympy.E)**((Integer(2) * sympy.I * x)))))) + ((Integer(2))**(Integer(-1)) * sympy.I * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * (sympy.E)**((Integer(2) * sympy.I * x))))) + sympy.tan(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_9():
    f = sin(x)/x
    F = sympy.Function('SinIntegral')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_10():
    f = (5*x**2 + 3*(x + exp(x))**(sympy.S(1)/3) + (2*x**2 + 3*x)*exp(x))/(x*(x + exp(x))**(sympy.S(1)/3))
    F = 3*x*(x + exp(x))**(sympy.S(2)/3) + 3*log(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_11():
    f = (1 + 1/x)/(x + log(x))**(sympy.S(3)/2) + 1/x
    F = log(x) - 2/sqrt(x + log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_12():
    f = (x**2 + 2*x*log(x) + (x + 1)*sqrt(x + log(x)) + log(x)**2)/(x**3 + 2*x**2*log(x) + x*log(x)**2)
    F = log(x) - 2/sqrt(x + log(x))
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_13():
    f = (-x**2 + 2*log(x)**2 - log(x))/(-x**2*log(x) + log(x)**3)
    F = ((Integer(-1) * (Integer(2))**(Integer(-1))) * sympy.log((x + (Integer(-1) * sympy.log(x))))) + ((Integer(2))**(Integer(-1)) * sympy.log((x + sympy.log(x)))) + sympy.Function('LogIntegral')(x)
    assert integrate(f, x) == F


def test_integrate_0_Independent_test_suites_Bronstein_Problems_14():
    f = (x**4 - 3*x**2 + 6)/(x**6 - 5*x**4 + 5*x**2 + 4)
    F = atan(x*(x**4 - 3*x**2 + 1)/2) + atan(2*x - sqrt(3)) + atan(2*x + sqrt(3))
    assert integrate(f, x) == F

