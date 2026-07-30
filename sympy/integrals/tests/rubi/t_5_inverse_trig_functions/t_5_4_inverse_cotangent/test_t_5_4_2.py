"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.4 Inverse cotangent/5.4.2 Exponentials of inverse cotangent.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, c, n = symbols('a c n')

def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_1():
    f = exp(acot(x))
    F = ((x - I)/x)**(1 + I/2)*((x + I)/x)**(-1 - I/2)*(sympy.S(4)/5 + 8*I/5)*hyper((2, 1 + I/2), (2 + I/2,), (1 - I/x)/(1 + I/x))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_2():
    f = exp(acot(x))/(a*x**2 + a)
    F = -exp(acot(x))/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_3():
    f = exp(acot(x))/(a*x**2 + a)**2
    F = -(1 - 2*x)*exp(acot(x))/(5*a**2*(x**2 + 1)) - 2*exp(acot(x))/(5*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_4():
    f = exp(acot(x))/(a*x**2 + a)**3
    F = -(1 - 4*x)*exp(acot(x))/(17*a**3*(x**2 + 1)**2) - 12*(1 - 2*x)*exp(acot(x))/(85*a**3*(x**2 + 1)) - 24*exp(acot(x))/(85*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_5():
    f = exp(acot(x))/(a*x**2 + a)**(sympy.S(3)/2)
    F = -(1 - x)*exp(acot(x))/(2*a*sqrt(a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_6():
    f = exp(acot(x))/(a*x**2 + a)**(sympy.S(5)/2)
    F = -(1 - 3*x)*exp(acot(x))/(10*a*(a*x**2 + a)**(sympy.S(3)/2)) - 3*(1 - x)*exp(acot(x))/(10*a**2*sqrt(a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_7():
    f = exp(acot(x))/(a*x**2 + a)**(sympy.S(7)/2)
    F = -(1 - 5*x)*exp(acot(x))/(26*a*(a*x**2 + a)**(sympy.S(5)/2)) - (1 - 3*x)*exp(acot(x))/(13*a**2*(a*x**2 + a)**(sympy.S(3)/2)) - 3*(1 - x)*exp(acot(x))/(13*a**3*sqrt(a*x**2 + a))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_8():
    f = exp(n*acot(a*x))/(a**2*c*x**2 + c)**(sympy.S(1)/3)
    F = 3*x*((a - I/x)/(a + I/x))**(-I*n/2 + sympy.S(1)/3)*(1 + 1/(a**2*x**2))**(sympy.S(1)/3)*(1 - I/(a*x))**(I*n/2 + sympy.S(-1)/3)*(1 + I/(a*x))**(-I*n/2 + sympy.S(2)/3)*hyper((sympy.S(-1)/3, -I*n/2 + sympy.S(1)/3), (sympy.S(2)/3,), 2*I/(x*(a + I/x)))/(a**2*c*x**2 + c)**(sympy.S(1)/3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_9():
    f = exp(n*acot(a*x))/(a**2*c*x**2 + c)**(sympy.S(2)/3)
    F = -3*x*((a - I/x)/(a + I/x))**(-I*n/2 + sympy.S(2)/3)*(1 + 1/(a**2*x**2))**(sympy.S(2)/3)*(1 - I/(a*x))**(I*n/2 + sympy.S(-2)/3)*(1 + I/(a*x))**(-I*n/2 + sympy.S(1)/3)*hyper((sympy.S(1)/3, -I*n/2 + sympy.S(2)/3), (sympy.S(4)/3,), 2*I/(x*(a + I/x)))/(a**2*c*x**2 + c)**(sympy.S(2)/3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_10():
    f = exp(n*acot(a*x))/(a**2*c*x**2 + c)**(sympy.S(4)/3)
    F = -6*x*((a - I/x)/(a + I/x))**(-I*n/2 + sympy.S(1)/3)*(1 + 1/(a**2*x**2))**(sympy.S(1)/3)*(1 - I/(a*x))**(I*n/2 + sympy.S(-1)/3)*(1 + I/(a*x))**(-I*n/2 + sympy.S(2)/3)*hyper((sympy.S(-1)/3, -I*n/2 + sympy.S(1)/3), (sympy.S(2)/3,), 2*I/(x*(a + I/x)))/(c*(9*n**2 + 4)*(a**2*c*x**2 + c)**(sympy.S(1)/3)) - 3*(-2*a*x + 3*n)*exp(n*acot(a*x))/(a*c*(9*n**2 + 4)*(a**2*c*x**2 + c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_11():
    f = exp(n*acot(a*x))/(a**2*c*x**2 + c)**(sympy.S(5)/3)
    F = -12*x*((a - I/x)/(a + I/x))**(-I*n/2 + sympy.S(2)/3)*(1 + 1/(a**2*x**2))**(sympy.S(2)/3)*(1 - I/(a*x))**(I*n/2 + sympy.S(-2)/3)*(1 + I/(a*x))**(-I*n/2 + sympy.S(1)/3)*hyper((sympy.S(1)/3, -I*n/2 + sympy.S(2)/3), (sympy.S(4)/3,), 2*I/(x*(a + I/x)))/(c*(9*n**2 + 16)*(a**2*c*x**2 + c)**(sympy.S(2)/3)) - 3*(-4*a*x + 3*n)*exp(n*acot(a*x))/(a*c*(9*n**2 + 16)*(a**2*c*x**2 + c)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_4_Inverse_cotangent_5_4_2_Exponentials_of_inverse_cotangent_12():
    f = exp(n*acot(a*x))/(a**2*c*x**2 + c)**(sympy.S(7)/3)
    F = -240*x*((a - I/x)/(a + I/x))**(-I*n/2 + sympy.S(1)/3)*(1 + 1/(a**2*x**2))**(sympy.S(1)/3)*(1 - I/(a*x))**(I*n/2 + sympy.S(-1)/3)*(1 + I/(a*x))**(-I*n/2 + sympy.S(2)/3)*hyper((sympy.S(-1)/3, -I*n/2 + sympy.S(1)/3), (sympy.S(2)/3,), 2*I/(x*(a + I/x)))/(c**2*(9*n**2 + 4)*(9*n**2 + 64)*(a**2*c*x**2 + c)**(sympy.S(1)/3)) - 3*(-8*a*x + 3*n)*exp(n*acot(a*x))/(a*c*(9*n**2 + 64)*(a**2*c*x**2 + c)**(sympy.S(4)/3)) - 120*(-2*a*x + 3*n)*exp(n*acot(a*x))/(a*c**2*(9*n**2 + 4)*(9*n**2 + 64)*(a**2*c*x**2 + c)**(sympy.S(1)/3))
    assert integrate(f, x) == F

