"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.2 Hyperbolic cosine/6.2.5 Hyperbolic cosine functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, C, F, a, b, c, d, e, f, m, n = symbols('A B C F a b c d e f m n')

def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_1():
    f = cosh(a + b*x)
    F = sinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_2():
    f = cosh(a + b*x)**2
    F = x/2 + sinh(a + b*x)*cosh(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_3():
    f = cosh(a + b*x)**3
    F = sinh(a + b*x)**3/(3*b) + sinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_4():
    f = cosh(a + b*x)**4
    F = 3*x/8 + sinh(a + b*x)*cosh(a + b*x)**3/(4*b) + 3*sinh(a + b*x)*cosh(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_5():
    f = cosh(a + b*x)**5
    F = sinh(a + b*x)**5/(5*b) + 2*sinh(a + b*x)**3/(3*b) + sinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_6():
    f = cosh(a + b*x)**6
    F = 5*x/16 + sinh(a + b*x)*cosh(a + b*x)**5/(6*b) + 5*sinh(a + b*x)*cosh(a + b*x)**3/(24*b) + 5*sinh(a + b*x)*cosh(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_7():
    f = cosh(a + b*x)**(sympy.S(7)/2)
    F = 2*sinh(a + b*x)*cosh(a + b*x)**(sympy.S(5)/2)/(7*b) + 10*sinh(a + b*x)*sqrt(cosh(a + b*x))/(21*b) - 10*I*elliptic_f(I*(a + b*x)/2, 2)/(21*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_8():
    f = cosh(a + b*x)**(sympy.S(5)/2)
    F = 2*sinh(a + b*x)*cosh(a + b*x)**(sympy.S(3)/2)/(5*b) - 6*I*elliptic_e(I*(a + b*x)/2, 2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_9():
    f = cosh(a + b*x)**(sympy.S(3)/2)
    F = 2*sinh(a + b*x)*sqrt(cosh(a + b*x))/(3*b) - 2*I*elliptic_f(I*(a + b*x)/2, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_10():
    f = sqrt(cosh(a + b*x))
    F = -2*I*elliptic_e(I*(a + b*x)/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_11():
    f = 1/sqrt(cosh(a + b*x))
    F = -2*I*elliptic_f(I*(a + b*x)/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_12():
    f = cosh(a + b*x)**(sympy.S(-3)/2)
    F = 2*sinh(a + b*x)/(b*sqrt(cosh(a + b*x))) + 2*I*elliptic_e(I*(a + b*x)/2, 2)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_13():
    f = cosh(a + b*x)**(sympy.S(-5)/2)
    F = 2*sinh(a + b*x)/(3*b*cosh(a + b*x)**(sympy.S(3)/2)) - 2*I*elliptic_f(I*(a + b*x)/2, 2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_14():
    f = cosh(a + b*x)**(sympy.S(-7)/2)
    F = 6*sinh(a + b*x)/(5*b*sqrt(cosh(a + b*x))) + 2*sinh(a + b*x)/(5*b*cosh(a + b*x)**(sympy.S(5)/2)) + 6*I*elliptic_e(I*(a + b*x)/2, 2)/(5*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_15():
    f = (a*cosh(x))**(sympy.S(7)/2)
    F = -10*I*a**4*sqrt(cosh(x))*elliptic_f(I*x/2, 2)/(21*sqrt(a*cosh(x))) + 10*a**3*sqrt(a*cosh(x))*sinh(x)/21 + 2*a*(a*cosh(x))**(sympy.S(5)/2)*sinh(x)/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_16():
    f = (a*cosh(x))**(sympy.S(5)/2)
    F = -6*I*a**2*sqrt(a*cosh(x))*elliptic_e(I*x/2, 2)/(5*sqrt(cosh(x))) + 2*a*(a*cosh(x))**(sympy.S(3)/2)*sinh(x)/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_17():
    f = (a*cosh(x))**(sympy.S(3)/2)
    F = -2*I*a**2*sqrt(cosh(x))*elliptic_f(I*x/2, 2)/(3*sqrt(a*cosh(x))) + 2*a*sqrt(a*cosh(x))*sinh(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_18():
    f = sqrt(a*cosh(x))
    F = -2*I*sqrt(a*cosh(x))*elliptic_e(I*x/2, 2)/sqrt(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_19():
    f = 1/sqrt(a*cosh(x))
    F = -2*I*sqrt(cosh(x))*elliptic_f(I*x/2, 2)/sqrt(a*cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_20():
    f = (a*cosh(x))**(sympy.S(-3)/2)
    F = 2*sinh(x)/(a*sqrt(a*cosh(x))) + 2*I*sqrt(a*cosh(x))*elliptic_e(I*x/2, 2)/(a**2*sqrt(cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_21():
    f = (a*cosh(x))**(sympy.S(-5)/2)
    F = 2*sinh(x)/(3*a*(a*cosh(x))**(sympy.S(3)/2)) - 2*I*sqrt(cosh(x))*elliptic_f(I*x/2, 2)/(3*a**2*sqrt(a*cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_22():
    f = (a*cosh(x))**(sympy.S(-7)/2)
    F = 2*sinh(x)/(5*a*(a*cosh(x))**(sympy.S(5)/2)) + 6*sinh(x)/(5*a**3*sqrt(a*cosh(x))) + 6*I*sqrt(a*cosh(x))*elliptic_e(I*x/2, 2)/(5*a**4*sqrt(cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_23():
    f = (b*cosh(c + d*x))**n
    F = -(b*cosh(c + d*x))**(n + 1)*sinh(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cosh(c + d*x)**2)/(b*d*sqrt(-sinh(c + d*x)**2)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_24():
    f = cosh(x)**4/(a*cosh(x) + a)
    F = -sinh(x)*cosh(x)**3/(a*cosh(x) + a) - 3*x/(2*a) + 4*sinh(x)**3/(3*a) - 3*sinh(x)*cosh(x)/(2*a) + 4*sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_25():
    f = cosh(x)**3/(a*cosh(x) + a)
    F = -sinh(x)*cosh(x)**2/(a*cosh(x) + a) + 3*x/(2*a) + 3*sinh(x)*cosh(x)/(2*a) - 2*sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_26():
    f = cosh(x)**2/(a*cosh(x) + a)
    F = -x/a + sinh(x)/a + sinh(x)/(a*(cosh(x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_27():
    f = cosh(x)/(a*cosh(x) + a)
    F = -sinh(x)/(a*cosh(x) + a) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_28():
    f = sech(x)/(a*cosh(x) + a)
    F = -sinh(x)/(a*cosh(x) + a) + atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_29():
    f = sech(x)**2/(a*cosh(x) + a)
    F = -tanh(x)/(a*cosh(x) + a) + 2*tanh(x)/a - atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_30():
    f = sech(x)**3/(a*cosh(x) + a)
    F = -tanh(x)*sech(x)/(a*cosh(x) + a) + 3*tanh(x)*sech(x)/(2*a) - 2*tanh(x)/a + 3*atan(sinh(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_31():
    f = sech(x)**4/(a*cosh(x) + a)
    F = -tanh(x)*sech(x)**2/(a*cosh(x) + a) - 4*tanh(x)**3/(3*a) - 3*tanh(x)*sech(x)/(2*a) + 4*tanh(x)/a - 3*atan(sinh(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_32():
    f = 1/(cosh(c + d*x) + 1)
    F = sinh(c + d*x)/(d*(cosh(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_33():
    f = (cosh(c + d*x) + 1)**(-2)
    F = sinh(c + d*x)/(3*d*(cosh(c + d*x) + 1)) + sinh(c + d*x)/(3*d*(cosh(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_34():
    f = (cosh(c + d*x) + 1)**(-3)
    F = 2*sinh(c + d*x)/(15*d*(cosh(c + d*x) + 1)) + 2*sinh(c + d*x)/(15*d*(cosh(c + d*x) + 1)**2) + sinh(c + d*x)/(5*d*(cosh(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_35():
    f = (cosh(c + d*x) + 1)**(-4)
    F = 2*sinh(c + d*x)/(35*d*(cosh(c + d*x) + 1)) + 2*sinh(c + d*x)/(35*d*(cosh(c + d*x) + 1)**2) + 3*sinh(c + d*x)/(35*d*(cosh(c + d*x) + 1)**3) + sinh(c + d*x)/(7*d*(cosh(c + d*x) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_36():
    f = 1/(1 - cosh(c + d*x))
    F = -sinh(c + d*x)/(d*(1 - cosh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_37():
    f = (1 - cosh(c + d*x))**(-2)
    F = -sinh(c + d*x)/(3*d*(1 - cosh(c + d*x))) - sinh(c + d*x)/(3*d*(1 - cosh(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_38():
    f = (1 - cosh(c + d*x))**(-3)
    F = -2*sinh(c + d*x)/(15*d*(1 - cosh(c + d*x))) - 2*sinh(c + d*x)/(15*d*(1 - cosh(c + d*x))**2) - sinh(c + d*x)/(5*d*(1 - cosh(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_39():
    f = (1 - cosh(c + d*x))**(-4)
    F = -2*sinh(c + d*x)/(35*d*(1 - cosh(c + d*x))) - 2*sinh(c + d*x)/(35*d*(1 - cosh(c + d*x))**2) - 3*sinh(c + d*x)/(35*d*(1 - cosh(c + d*x))**3) - sinh(c + d*x)/(7*d*(1 - cosh(c + d*x))**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_40():
    f = cosh(x)/sqrt(a*cosh(x) + a)
    F = 2*sinh(x)/sqrt(a*cosh(x) + a) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(a*cosh(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_41():
    f = cosh(x)/sqrt(-a*cosh(x) + a)
    F = 2*sinh(x)/sqrt(-a*cosh(x) + a) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(-a*cosh(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_42():
    f = (a*cosh(c + d*x) + a)**(sympy.S(5)/2)
    F = 64*a**3*sinh(c + d*x)/(15*d*sqrt(a*cosh(c + d*x) + a)) + 16*a**2*sqrt(a*cosh(c + d*x) + a)*sinh(c + d*x)/(15*d) + 2*a*(a*cosh(c + d*x) + a)**(sympy.S(3)/2)*sinh(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_43():
    f = (a*cosh(c + d*x) + a)**(sympy.S(3)/2)
    F = 8*a**2*sinh(c + d*x)/(3*d*sqrt(a*cosh(c + d*x) + a)) + 2*a*sqrt(a*cosh(c + d*x) + a)*sinh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_44():
    f = sqrt(a*cosh(c + d*x) + a)
    F = 2*a*sinh(c + d*x)/(d*sqrt(a*cosh(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_45():
    f = 1/sqrt(a*cosh(c + d*x) + a)
    F = sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(c + d*x)/(2*sqrt(a*cosh(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_46():
    f = (a*cosh(c + d*x) + a)**(sympy.S(-3)/2)
    F = sinh(c + d*x)/(2*d*(a*cosh(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(c + d*x)/(2*sqrt(a*cosh(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_47():
    f = (a*cosh(c + d*x) + a)**(sympy.S(-5)/2)
    F = sinh(c + d*x)/(4*d*(a*cosh(c + d*x) + a)**(sympy.S(5)/2)) + 3*sinh(c + d*x)/(16*a*d*(a*cosh(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(c + d*x)/(2*sqrt(a*cosh(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_48():
    f = (-a*cosh(c + d*x) + a)**(sympy.S(5)/2)
    F = -64*a**3*sinh(c + d*x)/(15*d*sqrt(-a*cosh(c + d*x) + a)) - 16*a**2*sqrt(-a*cosh(c + d*x) + a)*sinh(c + d*x)/(15*d) - 2*a*(-a*cosh(c + d*x) + a)**(sympy.S(3)/2)*sinh(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_49():
    f = (-a*cosh(c + d*x) + a)**(sympy.S(3)/2)
    F = -8*a**2*sinh(c + d*x)/(3*d*sqrt(-a*cosh(c + d*x) + a)) - 2*a*sqrt(-a*cosh(c + d*x) + a)*sinh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_50():
    f = sqrt(-a*cosh(c + d*x) + a)
    F = -2*a*sinh(c + d*x)/(d*sqrt(-a*cosh(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_51():
    f = 1/sqrt(-a*cosh(c + d*x) + a)
    F = -sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(c + d*x)/(2*sqrt(-a*cosh(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_52():
    f = (-a*cosh(c + d*x) + a)**(sympy.S(-3)/2)
    F = -sinh(c + d*x)/(2*d*(-a*cosh(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(c + d*x)/(2*sqrt(-a*cosh(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_53():
    f = (-a*cosh(c + d*x) + a)**(sympy.S(-5)/2)
    F = -sinh(c + d*x)/(4*d*(-a*cosh(c + d*x) + a)**(sympy.S(5)/2)) - 3*sinh(c + d*x)/(16*a*d*(-a*cosh(c + d*x) + a)**(sympy.S(3)/2)) - 3*sqrt(2)*atan(sqrt(2)*sqrt(a)*sinh(c + d*x)/(2*sqrt(-a*cosh(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_54():
    f = cosh(x)**4/(a + b*cosh(x))
    F = 2*a**4*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b**4*sqrt(a - b)*sqrt(a + b)) - a*sinh(x)*cosh(x)/(2*b**2) - a*x*(2*a**2 + b**2)/(2*b**4) + sinh(x)*cosh(x)**2/(3*b) + (3*a**2 + 2*b**2)*sinh(x)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_55():
    f = cosh(x)**3/(a + b*cosh(x))
    F = -2*a**3*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b**3*sqrt(a - b)*sqrt(a + b)) - a*sinh(x)/b**2 + sinh(x)*cosh(x)/(2*b) + x*(2*a**2 + b**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_56():
    f = cosh(x)**2/(a + b*cosh(x))
    F = 2*a**2*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b**2*sqrt(a - b)*sqrt(a + b)) - a*x/b**2 + sinh(x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_57():
    f = cosh(x)/(a + b*cosh(x))
    F = -2*a*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b*sqrt(a - b)*sqrt(a + b)) + x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_58():
    f = sech(x)/(a + b*cosh(x))
    F = -2*b*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*sqrt(a - b)*sqrt(a + b)) + atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_59():
    f = sech(x)**2/(a + b*cosh(x))
    F = tanh(x)/a + 2*b**2*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**2*sqrt(a - b)*sqrt(a + b)) - b*atan(sinh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_60():
    f = sech(x)**3/(a + b*cosh(x))
    F = tanh(x)*sech(x)/(2*a) - b*tanh(x)/a**2 - 2*b**3*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**3*sqrt(a - b)*sqrt(a + b)) + (a**2 + 2*b**2)*atan(sinh(x))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_61():
    f = sech(x)**4/(a + b*cosh(x))
    F = tanh(x)*sech(x)**2/(3*a) - b*tanh(x)*sech(x)/(2*a**2) + (2*a**2 + 3*b**2)*tanh(x)/(3*a**3) + 2*b**4*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**4*sqrt(a - b)*sqrt(a + b)) - b*(a**2 + 2*b**2)*atan(sinh(x))/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_62():
    f = (a + b*cosh(c + d*x))**5
    F = 7*a*b**2*(22*a**2 + 23*b**2)*sinh(c + d*x)*cosh(c + d*x)/(120*d) + 9*a*b*(a + b*cosh(c + d*x))**3*sinh(c + d*x)/(20*d) + a*x*(8*a**4 + 40*a**2*b**2 + 15*b**4)/8 + b*(a + b*cosh(c + d*x))**4*sinh(c + d*x)/(5*d) + b*(a + b*cosh(c + d*x))**2*(47*a**2 + 16*b**2)*sinh(c + d*x)/(60*d) + b*(107*a**4 + 192*a**2*b**2 + 16*b**4)*sinh(c + d*x)/(30*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_63():
    f = (a + b*cosh(c + d*x))**4
    F = 7*a*b*(a + b*cosh(c + d*x))**2*sinh(c + d*x)/(12*d) + a*b*(19*a**2 + 16*b**2)*sinh(c + d*x)/(6*d) + b**2*(26*a**2 + 9*b**2)*sinh(c + d*x)*cosh(c + d*x)/(24*d) + b*(a + b*cosh(c + d*x))**3*sinh(c + d*x)/(4*d) + x*(a**4 + 3*a**2*b**2 + 3*b**4/8)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_64():
    f = (a + b*cosh(c + d*x))**3
    F = 5*a*b**2*sinh(c + d*x)*cosh(c + d*x)/(6*d) + a*x*(2*a**2 + 3*b**2)/2 + b*(a + b*cosh(c + d*x))**2*sinh(c + d*x)/(3*d) + 2*b*(4*a**2 + b**2)*sinh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_65():
    f = (a + b*cosh(c + d*x))**2
    F = 2*a*b*sinh(c + d*x)/d + b**2*sinh(c + d*x)*cosh(c + d*x)/(2*d) + x*(a**2 + b**2/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_66():
    f = a + b*cosh(c + d*x)
    F = a*x + b*sinh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_67():
    f = 1/(a + b*cosh(c + d*x))
    F = 2*atanh(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_68():
    f = (a + b*cosh(c + d*x))**(-2)
    F = 2*a*atanh(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - b*sinh(c + d*x)/(d*(a + b*cosh(c + d*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_69():
    f = (a + b*cosh(c + d*x))**(-3)
    F = -3*a*b*sinh(c + d*x)/(2*d*(a + b*cosh(c + d*x))*(a**2 - b**2)**2) - b*sinh(c + d*x)/(d*(a + b*cosh(c + d*x))**2*(2*a**2 - 2*b**2)) + (2*a**2 + b**2)*atanh(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_70():
    f = (a + b*cosh(c + d*x))**(-4)
    F = -5*a*b*sinh(c + d*x)/(6*d*(a + b*cosh(c + d*x))**2*(a**2 - b**2)**2) + a*(2*a**2 + 3*b**2)*atanh(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - b*(11*a**2 + 4*b**2)*sinh(c + d*x)/(6*d*(a + b*cosh(c + d*x))*(a**2 - b**2)**3) - b*sinh(c + d*x)/(d*(a + b*cosh(c + d*x))**3*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_71():
    f = 1/(5*cosh(c + d*x) + 3)
    F = atan(tanh(c/2 + d*x/2)/2)/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_72():
    f = (5*cosh(c + d*x) + 3)**(-2)
    F = -3*atan(tanh(c/2 + d*x/2)/2)/(32*d) + 5*sinh(c + d*x)/(16*d*(5*cosh(c + d*x) + 3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_73():
    f = (5*cosh(c + d*x) + 3)**(-3)
    F = 43*atan(tanh(c/2 + d*x/2)/2)/(1024*d) - 45*sinh(c + d*x)/(512*d*(5*cosh(c + d*x) + 3)) + 5*sinh(c + d*x)/(32*d*(5*cosh(c + d*x) + 3)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_74():
    f = (5*cosh(c + d*x) + 3)**(-4)
    F = -279*atan(tanh(c/2 + d*x/2)/2)/(16384*d) + 995*sinh(c + d*x)/(24576*d*(5*cosh(c + d*x) + 3)) - 25*sinh(c + d*x)/(512*d*(5*cosh(c + d*x) + 3)**2) + 5*sinh(c + d*x)/(48*d*(5*cosh(c + d*x) + 3)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_75():
    f = 1/(3*cosh(c + d*x) + 5)
    F = x/4 - atanh(sinh(c + d*x)/(cosh(c + d*x) + 3))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_76():
    f = (3*cosh(c + d*x) + 5)**(-2)
    F = 5*x/64 - 5*atanh(sinh(c + d*x)/(cosh(c + d*x) + 3))/(32*d) - 3*sinh(c + d*x)/(16*d*(3*cosh(c + d*x) + 5))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_77():
    f = (3*cosh(c + d*x) + 5)**(-3)
    F = 59*x/2048 - 59*atanh(sinh(c + d*x)/(cosh(c + d*x) + 3))/(1024*d) - 45*sinh(c + d*x)/(512*d*(3*cosh(c + d*x) + 5)) - 3*sinh(c + d*x)/(32*d*(3*cosh(c + d*x) + 5)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_78():
    f = (3*cosh(c + d*x) + 5)**(-4)
    F = 385*x/32768 - 385*atanh(sinh(c + d*x)/(cosh(c + d*x) + 3))/(16384*d) - 311*sinh(c + d*x)/(8192*d*(3*cosh(c + d*x) + 5)) - 25*sinh(c + d*x)/(512*d*(3*cosh(c + d*x) + 5)**2) - sinh(c + d*x)/(16*d*(3*cosh(c + d*x) + 5)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_79():
    f = (a + b*cosh(x))**(sympy.S(5)/2)
    F = 16*a*b*sqrt(a + b*cosh(x))*sinh(x)/15 + 16*I*a*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)*elliptic_f(I*x/2, 2*b/(a + b))/(15*sqrt(a + b*cosh(x))) + 2*b*(a + b*cosh(x))**(sympy.S(3)/2)*sinh(x)/5 - 2*I*sqrt(a + b*cosh(x))*(23*a**2 + 9*b**2)*elliptic_e(I*x/2, 2*b/(a + b))/(15*sqrt((a + b*cosh(x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_80():
    f = (a + b*cosh(x))**(sympy.S(3)/2)
    F = -8*I*a*sqrt(a + b*cosh(x))*elliptic_e(I*x/2, 2*b/(a + b))/(3*sqrt((a + b*cosh(x))/(a + b))) + 2*b*sqrt(a + b*cosh(x))*sinh(x)/3 + 2*I*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)*elliptic_f(I*x/2, 2*b/(a + b))/(3*sqrt(a + b*cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_81():
    f = sqrt(a + b*cosh(c + d*x))
    F = -2*I*sqrt(a + b*cosh(c + d*x))*elliptic_e(I*(c + d*x)/2, 2*b/(a + b))/(d*sqrt((a + b*cosh(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_82():
    f = 1/sqrt(a + b*cosh(x))
    F = -2*I*sqrt((a + b*cosh(x))/(a + b))*elliptic_f(I*x/2, 2*b/(a + b))/sqrt(a + b*cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_83():
    f = (a + b*cosh(x))**(sympy.S(-3)/2)
    F = -2*b*sinh(x)/(sqrt(a + b*cosh(x))*(a**2 - b**2)) - 2*I*sqrt(a + b*cosh(x))*elliptic_e(I*x/2, 2*b/(a + b))/(sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_84():
    f = (a + b*cosh(x))**(sympy.S(-5)/2)
    F = -8*a*b*sinh(x)/(3*sqrt(a + b*cosh(x))*(a**2 - b**2)**2) - 8*I*a*sqrt(a + b*cosh(x))*elliptic_e(I*x/2, 2*b/(a + b))/(3*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)**2) - 2*b*sinh(x)/((a + b*cosh(x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + 2*I*sqrt((a + b*cosh(x))/(a + b))*elliptic_f(I*x/2, 2*b/(a + b))/(sqrt(a + b*cosh(x))*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_85():
    f = (a + b*cosh(x))**(sympy.S(-7)/2)
    F = -16*a*b*sinh(x)/(15*(a + b*cosh(x))**(sympy.S(3)/2)*(a**2 - b**2)**2) + 16*I*a*sqrt((a + b*cosh(x))/(a + b))*elliptic_f(I*x/2, 2*b/(a + b))/(15*sqrt(a + b*cosh(x))*(a**2 - b**2)**2) - 2*b*(23*a**2 + 9*b**2)*sinh(x)/(15*sqrt(a + b*cosh(x))*(a**2 - b**2)**3) - 2*b*sinh(x)/((a + b*cosh(x))**(sympy.S(5)/2)*(5*a**2 - 5*b**2)) - 2*I*sqrt(a + b*cosh(x))*(23*a**2 + 9*b**2)*elliptic_e(I*x/2, 2*b/(a + b))/(15*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_86():
    f = cosh(x)/sqrt(a + b*cosh(x))
    F = 2*I*a*sqrt((a + b*cosh(x))/(a + b))*elliptic_f(I*x/2, 2*b/(a + b))/(b*sqrt(a + b*cosh(x))) - 2*I*sqrt(a + b*cosh(x))*elliptic_e(I*x/2, 2*b/(a + b))/(b*sqrt((a + b*cosh(x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_87():
    f = (A + B*cosh(x))*(a*cosh(x) + a)**(sympy.S(5)/2)
    F = 2*B*(a*cosh(x) + a)**(sympy.S(5)/2)*sinh(x)/7 + 64*a**3*(7*A + 5*B)*sinh(x)/(105*sqrt(a*cosh(x) + a)) + 16*a**2*(7*A + 5*B)*sqrt(a*cosh(x) + a)*sinh(x)/105 + 2*a*(7*A + 5*B)*(a*cosh(x) + a)**(sympy.S(3)/2)*sinh(x)/35
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_88():
    f = (A + B*cosh(x))*(a*cosh(x) + a)**(sympy.S(3)/2)
    F = 2*B*(a*cosh(x) + a)**(sympy.S(3)/2)*sinh(x)/5 + 8*a**2*(5*A + 3*B)*sinh(x)/(15*sqrt(a*cosh(x) + a)) + 2*a*(5*A + 3*B)*sqrt(a*cosh(x) + a)*sinh(x)/15
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_89():
    f = (A + B*cosh(x))*sqrt(a*cosh(x) + a)
    F = 2*B*sqrt(a*cosh(x) + a)*sinh(x)/3 + 2*a*(3*A + B)*sinh(x)/(3*sqrt(a*cosh(x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_90():
    f = (A + B*cosh(x))*(-a*cosh(x) + a)**(sympy.S(5)/2)
    F = 2*B*(-a*cosh(x) + a)**(sympy.S(5)/2)*sinh(x)/7 - 64*a**3*(7*A - 5*B)*sinh(x)/(105*sqrt(-a*cosh(x) + a)) - 16*a**2*(7*A - 5*B)*sqrt(-a*cosh(x) + a)*sinh(x)/105 - 2*a*(7*A - 5*B)*(-a*cosh(x) + a)**(sympy.S(3)/2)*sinh(x)/35
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_91():
    f = (A + B*cosh(x))*(-a*cosh(x) + a)**(sympy.S(3)/2)
    F = 2*B*(-a*cosh(x) + a)**(sympy.S(3)/2)*sinh(x)/5 - 8*a**2*(5*A - 3*B)*sinh(x)/(15*sqrt(-a*cosh(x) + a)) - 2*a*(5*A - 3*B)*sqrt(-a*cosh(x) + a)*sinh(x)/15
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_92():
    f = (A + B*cosh(x))*sqrt(-a*cosh(x) + a)
    F = 2*B*sqrt(-a*cosh(x) + a)*sinh(x)/3 - 2*a*(3*A - B)*sinh(x)/(3*sqrt(-a*cosh(x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_93():
    f = (A + B*cosh(x))/(cosh(x) + 1)
    F = B*x + (A - B)*sinh(x)/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_94():
    f = (A + B*cosh(x))/(cosh(x) + 1)**2
    F = (A - B)*sinh(x)/(3*(cosh(x) + 1)**2) + (A + 2*B)*sinh(x)/(3*cosh(x) + 3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_95():
    f = (A + B*cosh(x))/(cosh(x) + 1)**3
    F = (A - B)*sinh(x)/(5*(cosh(x) + 1)**3) + (2*A + 3*B)*sinh(x)/(15*cosh(x) + 15) + (2*A + 3*B)*sinh(x)/(15*(cosh(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_96():
    f = (A + B*cosh(x))/(cosh(x) + 1)**4
    F = (A - B)*sinh(x)/(7*(cosh(x) + 1)**4) + (3*A + 4*B)*sinh(x)/(35*(cosh(x) + 1)**3) + (6*A + 8*B)*sinh(x)/(105*cosh(x) + 105) + (6*A + 8*B)*sinh(x)/(105*(cosh(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_97():
    f = (A + B*cosh(x))/(1 - cosh(x))
    F = -B*x - (A + B)*sinh(x)/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_98():
    f = (A + B*cosh(x))/(1 - cosh(x))**2
    F = -(A - 2*B)*sinh(x)/(3 - 3*cosh(x)) - (A + B)*sinh(x)/(3*(1 - cosh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_99():
    f = (A + B*cosh(x))/(1 - cosh(x))**3
    F = -(2*A - 3*B)*sinh(x)/(15 - 15*cosh(x)) - (2*A - 3*B)*sinh(x)/(15*(1 - cosh(x))**2) - (A + B)*sinh(x)/(5*(1 - cosh(x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_100():
    f = (A + B*cosh(x))/(1 - cosh(x))**4
    F = -(6*A - 8*B)*sinh(x)/(105 - 105*cosh(x)) - (6*A - 8*B)*sinh(x)/(105*(1 - cosh(x))**2) - (3*A - 4*B)*sinh(x)/(35*(1 - cosh(x))**3) - (A + B)*sinh(x)/(7*(1 - cosh(x))**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_101():
    f = (A + B*cosh(x))/sqrt(a*cosh(x) + a)
    F = 2*B*sinh(x)/sqrt(a*cosh(x) + a) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(a*cosh(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_102():
    f = (A + B*cosh(x))/(a*cosh(x) + a)**(sympy.S(3)/2)
    F = (A - B)*sinh(x)/(2*(a*cosh(x) + a)**(sympy.S(3)/2)) + sqrt(2)*(A + 3*B)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(a*cosh(x) + a)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_103():
    f = (A + B*cosh(x))/(a*cosh(x) + a)**(sympy.S(5)/2)
    F = (A - B)*sinh(x)/(4*(a*cosh(x) + a)**(sympy.S(5)/2)) + (3*A + 5*B)*sinh(x)/(16*a*(a*cosh(x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A + 5*B)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(a*cosh(x) + a)))/(32*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_104():
    f = (A + B*cosh(x))/sqrt(-a*cosh(x) + a)
    F = 2*B*sinh(x)/sqrt(-a*cosh(x) + a) - sqrt(2)*(A + B)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(-a*cosh(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_105():
    f = (A + B*cosh(x))/(-a*cosh(x) + a)**(sympy.S(3)/2)
    F = -(A + B)*sinh(x)/(2*(-a*cosh(x) + a)**(sympy.S(3)/2)) - sqrt(2)*(A - 3*B)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(-a*cosh(x) + a)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_106():
    f = (A + B*cosh(x))/(-a*cosh(x) + a)**(sympy.S(5)/2)
    F = -(A + B)*sinh(x)/(4*(-a*cosh(x) + a)**(sympy.S(5)/2)) - (3*A - 5*B)*sinh(x)/(16*a*(-a*cosh(x) + a)**(sympy.S(3)/2)) - sqrt(2)*(3*A - 5*B)*atan(sqrt(2)*sqrt(a)*sinh(x)/(2*sqrt(-a*cosh(x) + a)))/(32*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_107():
    f = (A + B*cosh(x))*(a + b*cosh(x))**(sympy.S(5)/2)
    F = 2*B*(a + b*cosh(x))**(sympy.S(5)/2)*sinh(x)/7 + (a + b*cosh(x))**(sympy.S(3)/2)*(2*A*b/5 + 2*B*a/7)*sinh(x) + sqrt(a + b*cosh(x))*(16*A*a*b/15 + 2*B*a**2/7 + 10*B*b**2/21)*sinh(x) + 2*I*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)*(56*A*a*b + 15*B*a**2 + 25*B*b**2)*elliptic_f(I*x/2, 2*b/(a + b))/(105*b*sqrt(a + b*cosh(x))) - 2*I*sqrt(a + b*cosh(x))*(161*A*a**2*b + 63*A*b**3 + 15*B*a**3 + 145*B*a*b**2)*elliptic_e(I*x/2, 2*b/(a + b))/(105*b*sqrt((a + b*cosh(x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_108():
    f = (A + B*cosh(x))*(a + b*cosh(x))**(sympy.S(3)/2)
    F = 2*B*(a + b*cosh(x))**(sympy.S(3)/2)*sinh(x)/5 + sqrt(a + b*cosh(x))*(2*A*b/3 + 2*B*a/5)*sinh(x) + 2*I*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)*(5*A*b + 3*B*a)*elliptic_f(I*x/2, 2*b/(a + b))/(15*b*sqrt(a + b*cosh(x))) - 2*I*sqrt(a + b*cosh(x))*(20*A*a*b + 3*B*a**2 + 9*B*b**2)*elliptic_e(I*x/2, 2*b/(a + b))/(15*b*sqrt((a + b*cosh(x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_109():
    f = (A + B*cosh(x))*sqrt(a + b*cosh(x))
    F = 2*B*sqrt(a + b*cosh(x))*sinh(x)/3 + 2*I*B*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)*elliptic_f(I*x/2, 2*b/(a + b))/(3*b*sqrt(a + b*cosh(x))) - 2*I*sqrt(a + b*cosh(x))*(3*A*b + B*a)*elliptic_e(I*x/2, 2*b/(a + b))/(3*b*sqrt((a + b*cosh(x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_110():
    f = (A + B*cosh(x))/(a + b*cosh(x))
    F = B*x/b + (2*A*b - 2*B*a)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_111():
    f = (A + B*cosh(x))/(a + b*cosh(x))**2
    F = -(A*b - B*a)*sinh(x)/((a + b*cosh(x))*(a**2 - b**2)) + (2*A*a - 2*B*b)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_112():
    f = (A + B*cosh(x))/(a + b*cosh(x))**3
    F = -(3*A*a*b - B*a**2 - 2*B*b**2)*sinh(x)/(2*(a + b*cosh(x))*(a**2 - b**2)**2) - (A*b - B*a)*sinh(x)/((a + b*cosh(x))**2*(2*a**2 - 2*b**2)) + (2*A*a**2 + A*b**2 - 3*B*a*b)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_113():
    f = (A + B*cosh(x))/(a + b*cosh(x))**4
    F = -(11*A*a**2*b + 4*A*b**3 - 2*B*a**3 - 13*B*a*b**2)*sinh(x)/(6*(a + b*cosh(x))*(a**2 - b**2)**3) - (5*A*a*b - 2*B*a**2 - 3*B*b**2)*sinh(x)/(6*(a + b*cosh(x))**2*(a**2 - b**2)**2) - (A*b - B*a)*sinh(x)/((a + b*cosh(x))**3*(3*a**2 - 3*b**2)) + (2*A*a**3 + 3*A*a*b**2 - 4*B*a**2*b - B*b**3)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_114():
    f = (B*cosh(x) + B*b/a)/(a + b*cosh(x))
    F = B*x/b - 2*B*sqrt(a - b)*sqrt(a + b)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_115():
    f = (B*a/b + B*cosh(x))/(a + b*cosh(x))
    F = B*x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_116():
    f = (a + b*cosh(x))/(a*cosh(x) + b)**2
    F = sinh(x)/(a*cosh(x) + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_117():
    f = (cosh(x) + 3)/(2 - cosh(x))
    F = -x + 5*sqrt(3)*x/3 + 10*sqrt(3)*atanh(sinh(x)/(-cosh(x) + sqrt(3) + 2))/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_118():
    f = (A + B*cosh(x))/sqrt(a + b*cosh(x))
    F = -2*I*B*sqrt(a + b*cosh(x))*elliptic_e(I*x/2, 2*b/(a + b))/(b*sqrt((a + b*cosh(x))/(a + b))) - 2*I*sqrt((a + b*cosh(x))/(a + b))*(A*b - B*a)*elliptic_f(I*x/2, 2*b/(a + b))/(b*sqrt(a + b*cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_119():
    f = (A + B*cosh(x))/(a + b*cosh(x))**(sympy.S(3)/2)
    F = -2*I*B*sqrt((a + b*cosh(x))/(a + b))*elliptic_f(I*x/2, 2*b/(a + b))/(b*sqrt(a + b*cosh(x))) - (2*A*b - 2*B*a)*sinh(x)/(sqrt(a + b*cosh(x))*(a**2 - b**2)) - 2*I*sqrt(a + b*cosh(x))*(A*b - B*a)*elliptic_e(I*x/2, 2*b/(a + b))/(b*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_120():
    f = (A + B*cosh(x))/(a + b*cosh(x))**(sympy.S(5)/2)
    F = -(8*A*a*b - 2*B*a**2 - 6*B*b**2)*sinh(x)/(3*sqrt(a + b*cosh(x))*(a**2 - b**2)**2) - (2*A*b - 2*B*a)*sinh(x)/((a + b*cosh(x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + 2*I*sqrt((a + b*cosh(x))/(a + b))*(A*b - B*a)*elliptic_f(I*x/2, 2*b/(a + b))/(3*b*sqrt(a + b*cosh(x))*(a**2 - b**2)) - 2*I*sqrt(a + b*cosh(x))*(4*A*a*b - B*a**2 - 3*B*b**2)*elliptic_e(I*x/2, 2*b/(a + b))/(3*b*sqrt((a + b*cosh(x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_121():
    f = (a*cosh(x)**2)**(sympy.S(7)/2)
    F = 16*a**3*sqrt(a*cosh(x)**2)*tanh(x)/35 + 8*a**2*(a*cosh(x)**2)**(sympy.S(3)/2)*tanh(x)/35 + 6*a*(a*cosh(x)**2)**(sympy.S(5)/2)*tanh(x)/35 + (a*cosh(x)**2)**(sympy.S(7)/2)*tanh(x)/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_122():
    f = (a*cosh(x)**2)**(sympy.S(5)/2)
    F = 8*a**2*sqrt(a*cosh(x)**2)*tanh(x)/15 + 4*a*(a*cosh(x)**2)**(sympy.S(3)/2)*tanh(x)/15 + (a*cosh(x)**2)**(sympy.S(5)/2)*tanh(x)/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_123():
    f = (a*cosh(x)**2)**(sympy.S(3)/2)
    F = 2*a*sqrt(a*cosh(x)**2)*tanh(x)/3 + (a*cosh(x)**2)**(sympy.S(3)/2)*tanh(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_124():
    f = sqrt(a*cosh(x)**2)
    F = sqrt(a*cosh(x)**2)*tanh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_125():
    f = 1/sqrt(a*cosh(x)**2)
    F = cosh(x)*atan(sinh(x))/sqrt(a*cosh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_126():
    f = (a*cosh(x)**2)**(sympy.S(-3)/2)
    F = cosh(x)*atan(sinh(x))/(2*a*sqrt(a*cosh(x)**2)) + tanh(x)/(2*a*sqrt(a*cosh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_127():
    f = (a*cosh(x)**2)**(sympy.S(-5)/2)
    F = tanh(x)/(4*a*(a*cosh(x)**2)**(sympy.S(3)/2)) + 3*cosh(x)*atan(sinh(x))/(8*a**2*sqrt(a*cosh(x)**2)) + 3*tanh(x)/(8*a**2*sqrt(a*cosh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_128():
    f = (a*cosh(x)**3)**(sympy.S(5)/2)
    F = 2*a**2*sqrt(a*cosh(x)**3)*sinh(x)*cosh(x)**5/15 + 26*a**2*sqrt(a*cosh(x)**3)*sinh(x)*cosh(x)**3/165 + 78*a**2*sqrt(a*cosh(x)**3)*sinh(x)*cosh(x)/385 + 26*a**2*sqrt(a*cosh(x)**3)*tanh(x)/77 - 26*I*a**2*sqrt(a*cosh(x)**3)*elliptic_f(I*x/2, 2)/(77*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_129():
    f = (a*cosh(x)**3)**(sympy.S(3)/2)
    F = 2*a*sqrt(a*cosh(x)**3)*sinh(x)*cosh(x)**2/9 + 14*a*sqrt(a*cosh(x)**3)*sinh(x)/45 - 14*I*a*sqrt(a*cosh(x)**3)*elliptic_e(I*x/2, 2)/(15*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_130():
    f = sqrt(a*cosh(x)**3)
    F = 2*sqrt(a*cosh(x)**3)*tanh(x)/3 - 2*I*sqrt(a*cosh(x)**3)*elliptic_f(I*x/2, 2)/(3*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_131():
    f = 1/sqrt(a*cosh(x)**3)
    F = 2*sinh(x)*cosh(x)/sqrt(a*cosh(x)**3) + 2*I*cosh(x)**(sympy.S(3)/2)*elliptic_e(I*x/2, 2)/sqrt(a*cosh(x)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_132():
    f = (a*cosh(x)**3)**(sympy.S(-3)/2)
    F = 10*sinh(x)/(21*a*sqrt(a*cosh(x)**3)) - 10*I*cosh(x)**(sympy.S(3)/2)*elliptic_f(I*x/2, 2)/(21*a*sqrt(a*cosh(x)**3)) + 2*tanh(x)*sech(x)/(7*a*sqrt(a*cosh(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_133():
    f = (a*cosh(x)**3)**(sympy.S(-5)/2)
    F = 154*sinh(x)*cosh(x)/(195*a**2*sqrt(a*cosh(x)**3)) + 154*I*cosh(x)**(sympy.S(3)/2)*elliptic_e(I*x/2, 2)/(195*a**2*sqrt(a*cosh(x)**3)) + 2*tanh(x)*sech(x)**4/(13*a**2*sqrt(a*cosh(x)**3)) + 22*tanh(x)*sech(x)**2/(117*a**2*sqrt(a*cosh(x)**3)) + 154*tanh(x)/(585*a**2*sqrt(a*cosh(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_134():
    f = (a*cosh(x)**4)**(sympy.S(5)/2)
    F = 63*a**2*x*sqrt(a*cosh(x)**4)*sech(x)**2/256 + a**2*sqrt(a*cosh(x)**4)*sinh(x)*cosh(x)**7/10 + 9*a**2*sqrt(a*cosh(x)**4)*sinh(x)*cosh(x)**5/80 + 21*a**2*sqrt(a*cosh(x)**4)*sinh(x)*cosh(x)**3/160 + 21*a**2*sqrt(a*cosh(x)**4)*sinh(x)*cosh(x)/128 + 63*a**2*sqrt(a*cosh(x)**4)*tanh(x)/256
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_135():
    f = (a*cosh(x)**4)**(sympy.S(3)/2)
    F = 5*a*x*sqrt(a*cosh(x)**4)*sech(x)**2/16 + a*sqrt(a*cosh(x)**4)*sinh(x)*cosh(x)**3/6 + 5*a*sqrt(a*cosh(x)**4)*sinh(x)*cosh(x)/24 + 5*a*sqrt(a*cosh(x)**4)*tanh(x)/16
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_136():
    f = sqrt(a*cosh(x)**4)
    F = x*sqrt(a*cosh(x)**4)*sech(x)**2/2 + sqrt(a*cosh(x)**4)*tanh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_137():
    f = 1/sqrt(a*cosh(x)**4)
    F = sinh(x)*cosh(x)/sqrt(a*cosh(x)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_138():
    f = (a*cosh(x)**4)**(sympy.S(-3)/2)
    F = sinh(x)**2*tanh(x)**3/(5*a*sqrt(a*cosh(x)**4)) - 2*sinh(x)**2*tanh(x)/(3*a*sqrt(a*cosh(x)**4)) + sinh(x)*cosh(x)/(a*sqrt(a*cosh(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_139():
    f = (a*cosh(x)**4)**(sympy.S(-5)/2)
    F = sinh(x)**2*tanh(x)**7/(9*a**2*sqrt(a*cosh(x)**4)) - 4*sinh(x)**2*tanh(x)**5/(7*a**2*sqrt(a*cosh(x)**4)) + 6*sinh(x)**2*tanh(x)**3/(5*a**2*sqrt(a*cosh(x)**4)) - 4*sinh(x)**2*tanh(x)/(3*a**2*sqrt(a*cosh(x)**4)) + sinh(x)*cosh(x)/(a**2*sqrt(a*cosh(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_140():
    f = sinh(x)/(cosh(x) + 1)**2
    F = -1/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_141():
    f = sinh(x)/(1 - cosh(x))**2
    F = 1/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_142():
    f = sinh(x)**2/(cosh(x) + 1)**2
    F = x - 2*sinh(x)/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_143():
    f = sinh(x)**2/(1 - cosh(x))**2
    F = x + 2*sinh(x)/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_144():
    f = sinh(x)**3/(cosh(x) + 1)**2
    F = -2*log(cosh(x) + 1) + cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_145():
    f = sinh(x)**3/(1 - cosh(x))**2
    F = 2*log(1 - cosh(x)) + cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_146():
    f = sinh(x)/(cosh(x) + 1)**3
    F = -1/(2*(cosh(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_147():
    f = sinh(x)/(1 - cosh(x))**3
    F = 1/(2*(1 - cosh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_148():
    f = sinh(x)**2/(cosh(x) + 1)**3
    F = sinh(x)**3/(3*(cosh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_149():
    f = sinh(x)**2/(1 - cosh(x))**3
    F = -sinh(x)**3/(3*(1 - cosh(x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_150():
    f = sinh(x)**3/(cosh(x) + 1)**3
    F = log(cosh(x) + 1) + 2/(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_151():
    f = sinh(x)**3/(1 - cosh(x))**3
    F = -log(1 - cosh(x)) - 2/(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_152():
    f = sinh(x)**8/(a*cosh(x) + a)
    F = 5*x/(16*a) + sinh(x)**7/(7*a) - sinh(x)**5*cosh(x)/(6*a) + 5*sinh(x)**3*cosh(x)/(24*a) - 5*sinh(x)*cosh(x)/(16*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_153():
    f = sinh(x)**7/(a*cosh(x) + a)
    F = (-a*cosh(x) + a)**4/a**5 - 4*(-a*cosh(x) + a)**5/(5*a**6) + (-a*cosh(x) + a)**6/(6*a**7)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_154():
    f = sinh(x)**6/(a*cosh(x) + a)
    F = -3*x/(8*a) + sinh(x)**5/(5*a) - sinh(x)**3*cosh(x)/(4*a) + 3*sinh(x)*cosh(x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_155():
    f = sinh(x)**5/(a*cosh(x) + a)
    F = -2*(-a*cosh(x) + a)**3/(3*a**4) + (-a*cosh(x) + a)**4/(4*a**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_156():
    f = sinh(x)**4/(a*cosh(x) + a)
    F = x/(2*a) + sinh(x)**3/(3*a) - sinh(x)*cosh(x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_157():
    f = sinh(x)**3/(a*cosh(x) + a)
    F = cosh(x)**2/(2*a) - cosh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_158():
    f = sinh(x)**2/(a*cosh(x) + a)
    F = -x/a + sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_159():
    f = sinh(x)/(a*cosh(x) + a)
    F = log(cosh(x) + 1)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_160():
    f = csch(x)/(a*cosh(x) + a)
    F = 1/(2*a*cosh(x) + 2*a) - atanh(cosh(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_161():
    f = csch(x)**2/(a*cosh(x) + a)
    F = csch(x)/(3*a*cosh(x) + 3*a) - 2*coth(x)/(3*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_162():
    f = csch(x)**3/(a*cosh(x) + a)
    F = -a/(8*(a*cosh(x) + a)**2) - 1/(4*a*cosh(x) + 4*a) + 1/(-8*a*cosh(x) + 8*a) + 3*atanh(cosh(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_163():
    f = csch(x)**4/(a*cosh(x) + a)
    F = csch(x)**3/(5*a*cosh(x) + 5*a) - 4*coth(x)**3/(15*a) + 4*coth(x)/(5*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_164():
    f = csch(x)**5/(a*cosh(x) + a)
    F = a**2/(24*(a*cosh(x) + a)**3) + 3*a/(32*(a*cosh(x) + a)**2) - a/(32*(-a*cosh(x) + a)**2) + 3/(16*a*cosh(x) + 16*a) - 1/(-8*a*cosh(x) + 8*a) - 5*atanh(cosh(x))/(16*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_165():
    f = sinh(x)**7/(a + b*cosh(x))
    F = -a*cosh(x)**5/(5*b**2) - a*(a**2 - 3*b**2)*cosh(x)**3/(3*b**4) - a*(a**4 - 3*a**2*b**2 + 3*b**4)*cosh(x)/b**6 + cosh(x)**6/(6*b) + (a**2 - 3*b**2)*cosh(x)**4/(4*b**3) + (a**4 - 3*a**2*b**2 + 3*b**4)*cosh(x)**2/(2*b**5) + (a**2 - b**2)**3*log(a + b*cosh(x))/b**7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_166():
    f = sinh(x)**6/(a + b*cosh(x))
    F = -a*x*(8*a**4 - 20*a**2*b**2 + 15*b**4)/(8*b**6) + sinh(x)**5/(5*b) + (4*a**2 - 3*a*b*cosh(x) - 4*b**2)*sinh(x)**3/(12*b**3) + (-a*b*(4*a**2 - 7*b**2)*cosh(x) + 8*(a**2 - b**2)**2)*sinh(x)/(8*b**5) + 2*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/b**6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_167():
    f = sinh(x)**5/(a + b*cosh(x))
    F = -a*cosh(x)**3/(3*b**2) - a*(a**2 - 2*b**2)*cosh(x)/b**4 + cosh(x)**4/(4*b) + (a**2 - 2*b**2)*cosh(x)**2/(2*b**3) + (a**2 - b**2)**2*log(a + b*cosh(x))/b**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_168():
    f = sinh(x)**4/(a + b*cosh(x))
    F = -a*x*(2*a**2 - 3*b**2)/(2*b**4) + sinh(x)**3/(3*b) + (2*a**2 - a*b*cosh(x) - 2*b**2)*sinh(x)/(2*b**3) + 2*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/b**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_169():
    f = sinh(x)**3/(a + b*cosh(x))
    F = -a*cosh(x)/b**2 + cosh(x)**2/(2*b) + (a**2 - b**2)*log(a + b*cosh(x))/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_170():
    f = sinh(x)**2/(a + b*cosh(x))
    F = -a*x/b**2 + sinh(x)/b + 2*sqrt(a - b)*sqrt(a + b)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_171():
    f = sinh(x)/(a + b*cosh(x))
    F = log(a + b*cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_172():
    f = csch(x)/(a + b*cosh(x))
    F = b*log(a + b*cosh(x))/(a**2 - b**2) + log(1 - cosh(x))/(2*a + 2*b) - log(cosh(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_173():
    f = csch(x)**2/(a + b*cosh(x))
    F = 2*b**2*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + (-a*cosh(x) + b)*csch(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_174():
    f = csch(x)**3/(a + b*cosh(x))
    F = b**3*log(a + b*cosh(x))/(a**2 - b**2)**2 + (a - 2*b)*log(cosh(x) + 1)/(4*(a - b)**2) + (-a*cosh(x) + b)*csch(x)**2/(2*a**2 - 2*b**2) - (a + 2*b)*log(1 - cosh(x))/(4*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_175():
    f = csch(x)**4/(a + b*cosh(x))
    F = 2*b**4*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + (-a*cosh(x) + b)*csch(x)**3/(3*a**2 - 3*b**2) + (a*(2*a**2 - 5*b**2)*cosh(x) + 3*b**3)*csch(x)/(3*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_176():
    f = csch(x)**5/(a + b*cosh(x))
    F = b**5*log(a + b*cosh(x))/(a**2 - b**2)**3 + (-a*cosh(x) + b)*csch(x)**4/(4*a**2 - 4*b**2) + (a*(3*a**2 - 7*b**2)*cosh(x) + 4*b**3)*csch(x)**2/(8*(a**2 - b**2)**2) + (3*a**2 + 9*a*b + 8*b**2)*log(1 - cosh(x))/(16*(a + b)**3) - (3*a**2 - 9*a*b + 8*b**2)*log(cosh(x) + 1)/(16*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_177():
    f = csch(x)**6/(a + b*cosh(x))
    F = 2*b**6*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (-a*cosh(x) + b)*csch(x)**5/(5*a**2 - 5*b**2) + (a*(4*a**2 - 9*b**2)*cosh(x) + 5*b**3)*csch(x)**3/(15*(a**2 - b**2)**2) + (-a*(8*a**4 - 26*a**2*b**2 + 33*b**4)*cosh(x) + 15*b**5)*csch(x)/(15*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_178():
    f = sinh(x)**2/(a + b*cosh(x))**2
    F = -2*a*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b**2*sqrt(a - b)*sqrt(a + b)) - sinh(x)/(b*(a + b*cosh(x))) + x/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_179():
    f = tanh(x)**4/(a + b*cosh(x))
    F = tanh(x)*sech(x)**2/(3*a) - b*tanh(x)*sech(x)/(2*a**2) - (4*a**2 - 3*b**2)*tanh(x)/(3*a**3) + b*(3*a**2 - 2*b**2)*atan(sinh(x))/(2*a**4) + 2*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_180():
    f = tanh(x)**3/(a + b*cosh(x))
    F = sech(x)**2/(2*a) - b*sech(x)/a**2 - (a**2 - b**2)*log(a + b*cosh(x))/a**3 + (a**2 - b**2)*log(cosh(x))/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_181():
    f = tanh(x)**2/(a + b*cosh(x))
    F = -tanh(x)/a + b*atan(sinh(x))/a**2 + 2*sqrt(a - b)*sqrt(a + b)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_182():
    f = tanh(x)/(a + b*cosh(x))
    F = -log(a + b*cosh(x))/a + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_183():
    f = coth(x)/(a + b*cosh(x))
    F = -a*log(a + b*cosh(x))/(a**2 - b**2) + log(1 - cosh(x))/(2*a + 2*b) + log(cosh(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_184():
    f = coth(x)**2/(a + b*cosh(x))
    F = 2*a**2*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - a*coth(x)/(a**2 - b**2) + b*csch(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_185():
    f = coth(x)**3/(a + b*cosh(x))
    F = -a**3*log(a + b*cosh(x))/(a**2 - b**2)**2 - (a - b*cosh(x))*csch(x)**2/(2*a**2 - 2*b**2) + (2*a + b)*log(1 - cosh(x))/(4*(a + b)**2) + (2*a - b)*log(cosh(x) + 1)/(4*(a - b)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_186():
    f = coth(x)**4/(a + b*cosh(x))
    F = 2*a**4*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - a**3*coth(x)/(a**2 - b**2)**2 + a**2*b*csch(x)/(a**2 - b**2)**2 - a*coth(x)**3/(3*a**2 - 3*b**2) + b*csch(x)**3/(3*a**2 - 3*b**2) + b*csch(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_187():
    f = tanh(x)**6/(a*cosh(x) + a)
    F = -tanh(x)**5/(5*a) - tanh(x)**3*sech(x)/(4*a) - 3*tanh(x)*sech(x)/(8*a) + 3*atan(sinh(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_188():
    f = tanh(x)**5/(a*cosh(x) + a)
    F = -tanh(x)**4/(4*a) + sech(x)**3/(3*a) - sech(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_189():
    f = tanh(x)**4/(a*cosh(x) + a)
    F = -tanh(x)**3/(3*a) - tanh(x)*sech(x)/(2*a) + atan(sinh(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_190():
    f = tanh(x)**3/(a*cosh(x) + a)
    F = sech(x)**2/(2*a) - sech(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_191():
    f = tanh(x)**2/(a*cosh(x) + a)
    F = -tanh(x)/a + atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_192():
    f = tanh(x)/(a*cosh(x) + a)
    F = -log(cosh(x) + 1)/a + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_193():
    f = coth(x)/(a*cosh(x) + a)
    F = -coth(x)*csch(x)/(2*a) - atanh(cosh(x))/(2*a) + csch(x)**2/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_194():
    f = coth(x)**2/(a*cosh(x) + a)
    F = coth(x)**3/(3*a) - csch(x)**3/(3*a) - csch(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_195():
    f = coth(x)**3/(a*cosh(x) + a)
    F = coth(x)**4/(4*a) - coth(x)**3*csch(x)/(4*a) - 3*coth(x)*csch(x)/(8*a) - 3*atanh(cosh(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_196():
    f = coth(x)**4/(a*cosh(x) + a)
    F = coth(x)**5/(5*a) - csch(x)**5/(5*a) - 2*csch(x)**3/(3*a) - csch(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_197():
    f = sqrt(a + b*cosh(x))*tanh(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*cosh(x))/sqrt(a)) + 2*sqrt(a + b*cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_198():
    f = tanh(x)/sqrt(a + b*cosh(x))
    F = -2*atanh(sqrt(a + b*cosh(x))/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_199():
    f = (A + B*sinh(x))/(a + b*cosh(x))
    F = 2*A*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(sqrt(a - b)*sqrt(a + b)) + B*log(a + b*cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_200():
    f = (A + B*sinh(x))/(cosh(x) + 1)
    F = A*sinh(x)/(cosh(x) + 1) + B*log(cosh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_201():
    f = (A + B*sinh(x))/(1 - cosh(x))
    F = -A*sinh(x)/(1 - cosh(x)) - B*log(1 - cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_202():
    f = (A + B*tanh(x))/(a + b*cosh(x))
    F = 2*A*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(sqrt(a - b)*sqrt(a + b)) - B*log(a + b*cosh(x))/a + B*log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_203():
    f = (A + B*coth(x))/(a + b*cosh(x))
    F = 2*A*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(sqrt(a - b)*sqrt(a + b)) - B*a*log(a + b*cosh(x))/(a**2 - b**2) + B*log(1 - cosh(x))/(2*a + 2*b) + B*log(cosh(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_204():
    f = (A + B*sech(x))/(a + b*cosh(x))
    F = B*atan(sinh(x))/a + (2*A*a - 2*B*b)*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_205():
    f = (A + B*csch(x))/(a + b*cosh(x))
    F = 2*A*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(sqrt(a - b)*sqrt(a + b)) + B*b*log(a + b*cosh(x))/(a**2 - b**2) + B*log(1 - cosh(x))/(2*a + 2*b) - B*log(cosh(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_206():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + b*cosh(d + e*x))
    F = B*x/b + C*log(a + b*cosh(d + e*x))/(b*e) + (2*A*b - 2*B*a)*atanh(sqrt(a - b)*tanh(d/2 + e*x/2)/sqrt(a + b))/(b*e*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_207():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + b*cosh(d + e*x))**2
    F = -C/(b*e*(a + b*cosh(d + e*x))) - (A*b - B*a)*sinh(d + e*x)/(e*(a + b*cosh(d + e*x))*(a**2 - b**2)) + (2*A*a - 2*B*b)*atanh(sqrt(a - b)*tanh(d/2 + e*x/2)/sqrt(a + b))/(e*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_208():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + b*cosh(d + e*x))**3
    F = -C/(2*b*e*(a + b*cosh(d + e*x))**2) - (3*A*a*b - B*a**2 - 2*B*b**2)*sinh(d + e*x)/(2*e*(a + b*cosh(d + e*x))*(a**2 - b**2)**2) - (A*b - B*a)*sinh(d + e*x)/(e*(a + b*cosh(d + e*x))**2*(2*a**2 - 2*b**2)) + (2*A*a**2 + A*b**2 - 3*B*a*b)*atanh(sqrt(a - b)*tanh(d/2 + e*x/2)/sqrt(a + b))/(e*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_209():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + b*cosh(d + e*x))**4
    F = -C/(3*b*e*(a + b*cosh(d + e*x))**3) - (11*A*a**2*b + 4*A*b**3 - 2*B*a**3 - 13*B*a*b**2)*sinh(d + e*x)/(6*e*(a + b*cosh(d + e*x))*(a**2 - b**2)**3) - (5*A*a*b - 2*B*a**2 - 3*B*b**2)*sinh(d + e*x)/(6*e*(a + b*cosh(d + e*x))**2*(a**2 - b**2)**2) - (A*b - B*a)*sinh(d + e*x)/(e*(a + b*cosh(d + e*x))**3*(3*a**2 - 3*b**2)) + (2*A*a**3 + 3*A*a*b**2 - 4*B*a**2*b - B*b**3)*atanh(sqrt(a - b)*tanh(d/2 + e*x/2)/sqrt(a + b))/(e*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_210():
    f = x/(a + b*cosh(x)**2)
    F = ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_211():
    f = x**2/(a + b*cosh(x)**2)
    F = (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_212():
    f = x**3/(a + b*cosh(x)**2)
    F = (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))))**(Integer(-1)))))) * ((Integer(8) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + Symbol('b') + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b'))))))**(Integer(-1)))))) * ((Integer(8) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + Symbol('b')))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_213():
    f = cosh(sqrt(-a*x + 1)/sqrt(a*x + 1))**3/(-a**2*x**2 + 1)
    F = (Integer(-1) * ((Integer(3) * sympy.Function('CoshIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(4) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(4) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_214():
    f = cosh(sqrt(-a*x + 1)/sqrt(a*x + 1))**2/(-a**2*x**2 + 1)
    F = (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (Integer(-1) * (sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_215():
    f = cosh(sqrt(-a*x + 1)/sqrt(a*x + 1))/(-a**2*x**2 + 1)
    F = Integer(-1) * (sympy.Function('CoshIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_216():
    f = 1/((-a**2*x**2 + 1)*cosh(sqrt(-a*x + 1)/sqrt(a*x + 1)))
    F = sympy.Function('Unintegrable')((sympy.sech((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_217():
    f = 1/((-a**2*x**2 + 1)*cosh(sqrt(-a*x + 1)/sqrt(a*x + 1))**2)
    F = sympy.Function('Unintegrable')(((sympy.sech((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))))**(Integer(2)) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_218():
    f = x*sinh(x)/(a + b*cosh(x))**2
    F = -x/(b*(a + b*cosh(x))) + 2*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_219():
    f = x*sinh(x)/(a + b*cosh(x))**3
    F = a*atanh(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - sinh(x)/((a + b*cosh(x))*(2*a**2 - 2*b**2)) - x/(2*b*(a + b*cosh(x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_220():
    f = (cosh(a + b*x)**2 + 2)*sinh(a + b*x)/x
    F = ((Integer(9) * (Integer(4))**(Integer(-1))) * sympy.Function('CoshIntegral')((Symbol('b') * x)) * sympy.sinh(Symbol('a'))) + ((Integer(4))**(Integer(-1)) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x)) * sympy.sinh((Integer(3) * Symbol('a')))) + ((Integer(9) * (Integer(4))**(Integer(-1))) * sympy.cosh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x))) + ((Integer(4))**(Integer(-1)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_221():
    f = x**m*sinh(c + d*x)/(a + b*cosh(c + d*x))
    F = sympy.Function('Unintegrable')((((x)**(Symbol('m')) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_222():
    f = x**3*sinh(c + d*x)/(a + b*cosh(c + d*x))
    F = ((Integer(-1) * (x)**(Integer(4))) * ((Integer(4) * Symbol('b')))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + ((Integer(6) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(4))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_223():
    f = x**2*sinh(c + d*x)/(a + b*cosh(c + d*x))
    F = ((Integer(-1) * (x)**(Integer(3))) * ((Integer(3) * Symbol('b')))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_224():
    f = x*sinh(c + d*x)/(a + b*cosh(c + d*x))
    F = ((Integer(-1) * (x)**(Integer(2))) * ((Integer(2) * Symbol('b')))**(Integer(-1))) + ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_225():
    f = sinh(c + d*x)/(a + b*cosh(c + d*x))
    F = log(a + b*cosh(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_226():
    f = sinh(c + d*x)/(x*(a + b*cosh(c + d*x)))
    F = sympy.Function('Unintegrable')((sympy.sinh((Symbol('c') + (Symbol('d') * x))) * ((x * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_227():
    f = x**m*sinh(c + d*x)**2/(a + b*cosh(c + d*x))
    F = sympy.Function('Unintegrable')((((x)**(Symbol('m')) * (sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_228():
    f = x**3*sinh(c + d*x)**2/(a + b*cosh(c + d*x))
    F = ((Integer(-1) * (Symbol('a') * (x)**(Integer(4)))) * ((Integer(4) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(6) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + ((Integer(6) * x * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(3)) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_229():
    f = x**2*sinh(c + d*x)**2/(a + b*cosh(c + d*x))
    F = ((Integer(-1) * (Symbol('a') * (x)**(Integer(3)))) * ((Integer(3) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * x * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((Integer(2) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_230():
    f = x*sinh(c + d*x)**2/(a + b*cosh(c + d*x))
    F = ((Integer(-1) * (Symbol('a') * (x)**(Integer(2)))) * ((Integer(2) * (Symbol('b'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * (sympy.cosh((Symbol('c') + (Symbol('d') * x))) * ((Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((x * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_231():
    f = sinh(c + d*x)**2/(a + b*cosh(c + d*x))
    F = -a*x/b**2 + sinh(c + d*x)/(b*d) + 2*sqrt(a - b)*sqrt(a + b)*atanh(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_232():
    f = sinh(c + d*x)**2/(x*(a + b*cosh(c + d*x)))
    F = sympy.Function('Unintegrable')(((sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * ((x * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_233():
    f = x**m*sinh(c + d*x)**3/(a + b*cosh(c + d*x))
    F = sympy.Function('Unintegrable')((((x)**(Symbol('m')) * (sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(3))) * ((Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_234():
    f = x**3*sinh(c + d*x)**3/(a + b*cosh(c + d*x))
    F = ((Integer(3) * x) * ((Integer(8) * Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + ((x)**(Integer(3)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(4))) * ((Integer(4) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * Symbol('a') * x * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (x)**(Integer(3)) * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(6) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(6) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + ((Integer(6) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + ((Integer(6) * Symbol('a') * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(4))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * (x)**(Integer(2)) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.cosh((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(8) * Symbol('b') * (Symbol('d'))**(Integer(4))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.cosh((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * x * (sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((Integer(4) * Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(3)) * (sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_235():
    f = x**2*sinh(c + d*x)**3/(a + b*cosh(c + d*x))
    F = ((x)**(Integer(2)) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(3))) * ((Integer(3) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (x)**(Integer(2)) * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(3))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * x * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.cosh((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * ((Integer(4) * Symbol('b') * (Symbol('d'))**(Integer(3))))**(Integer(-1))) + (((x)**(Integer(2)) * (sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_236():
    f = x*sinh(c + d*x)**3/(a + b*cosh(c + d*x))
    F = (x * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (x)**(Integer(2))) * ((Integer(2) * (Symbol('b'))**(Integer(3))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * x * sympy.cosh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('c') + (Symbol('d') * x))) * sympy.sinh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((x * (sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(2))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_237():
    f = sinh(c + d*x)**3/(a + b*cosh(c + d*x))
    F = -a*cosh(c + d*x)/(b**2*d) + cosh(c + d*x)**2/(2*b*d) + (a**2 - b**2)*log(a + b*cosh(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_238():
    f = sinh(c + d*x)**3/(x*(a + b*cosh(c + d*x)))
    F = sympy.Function('Unintegrable')(((sympy.sinh((Symbol('c') + (Symbol('d') * x))))**(Integer(3)) * ((x * (Symbol('a') + (Symbol('b') * sympy.cosh((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_239():
    f = cosh(a + b*log(c*x**n))
    F = -b*n*x*sinh(a + b*log(c*x**n))/(-b**2*n**2 + 1) + x*cosh(a + b*log(c*x**n))/(-b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_240():
    f = cosh(a + b*log(c*x**n))**2
    F = -2*b**2*n**2*x/(-4*b**2*n**2 + 1) - 2*b*n*x*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(-4*b**2*n**2 + 1) + x*cosh(a + b*log(c*x**n))**2/(-4*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_241():
    f = cosh(a + b*log(c*x**n))**3
    F = 6*b**3*n**3*x*sinh(a + b*log(c*x**n))/(9*b**4*n**4 - 10*b**2*n**2 + 1) - 6*b**2*n**2*x*cosh(a + b*log(c*x**n))/(9*b**4*n**4 - 10*b**2*n**2 + 1) - 3*b*n*x*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**2/(-9*b**2*n**2 + 1) + x*cosh(a + b*log(c*x**n))**3/(-9*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_242():
    f = cosh(a + b*log(c*x**n))**4
    F = 24*b**4*n**4*x/(64*b**4*n**4 - 20*b**2*n**2 + 1) + 24*b**3*n**3*x*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(64*b**4*n**4 - 20*b**2*n**2 + 1) - 12*b**2*n**2*x*cosh(a + b*log(c*x**n))**2/(64*b**4*n**4 - 20*b**2*n**2 + 1) - 4*b*n*x*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**3/(-16*b**2*n**2 + 1) + x*cosh(a + b*log(c*x**n))**4/(-16*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_243():
    f = x**m*cosh(a + b*log(c*x**n))
    F = -b*n*x**(m + 1)*sinh(a + b*log(c*x**n))/(-b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cosh(a + b*log(c*x**n))/(-b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_244():
    f = x**m*cosh(a + b*log(c*x**n))**2
    F = -2*b**2*n**2*x**(m + 1)/((m + 1)*(-4*b**2*n**2 + (m + 1)**2)) - 2*b*n*x**(m + 1)*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(-4*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cosh(a + b*log(c*x**n))**2/(-4*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_245():
    f = x**m*cosh(a + b*log(c*x**n))**3
    F = 6*b**3*n**3*x**(m + 1)*sinh(a + b*log(c*x**n))/((-9*b**2*n**2 + (m + 1)**2)*(-b**2*n**2 + (m + 1)**2)) - 6*b**2*n**2*x**(m + 1)*(m + 1)*cosh(a + b*log(c*x**n))/((-9*b**2*n**2 + (m + 1)**2)*(-b**2*n**2 + (m + 1)**2)) - 3*b*n*x**(m + 1)*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**2/(-9*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cosh(a + b*log(c*x**n))**3/(-9*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_246():
    f = x**m*cosh(a + b*log(c*x**n))**4
    F = 24*b**4*n**4*x**(m + 1)/((m + 1)*(-16*b**2*n**2 + (m + 1)**2)*(-4*b**2*n**2 + (m + 1)**2)) + 24*b**3*n**3*x**(m + 1)*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/((-16*b**2*n**2 + (m + 1)**2)*(-4*b**2*n**2 + (m + 1)**2)) - 12*b**2*n**2*x**(m + 1)*(m + 1)*cosh(a + b*log(c*x**n))**2/((-16*b**2*n**2 + (m + 1)**2)*(-4*b**2*n**2 + (m + 1)**2)) - 4*b*n*x**(m + 1)*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**3/(-16*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*cosh(a + b*log(c*x**n))**4/(-16*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_247():
    f = cosh(a + b*log(c*x**n))/x
    F = sinh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_248():
    f = cosh(a + b*log(c*x**n))**2/x
    F = log(x)/2 + sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_249():
    f = cosh(a + b*log(c*x**n))**3/x
    F = sinh(a + b*log(c*x**n))**3/(3*b*n) + sinh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_250():
    f = cosh(a + b*log(c*x**n))**4/x
    F = 3*log(x)/8 + sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**3/(4*b*n) + 3*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(8*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_251():
    f = cosh(a + b*log(c*x**n))**5/x
    F = sinh(a + b*log(c*x**n))**5/(5*b*n) + 2*sinh(a + b*log(c*x**n))**3/(3*b*n) + sinh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_252():
    f = cosh(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = 2*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))**(sympy.S(3)/2)/(5*b*n) - 6*I*elliptic_e(I*(a + b*log(c*x**n))/2, 2)/(5*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_253():
    f = cosh(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = 2*sinh(a + b*log(c*x**n))*sqrt(cosh(a + b*log(c*x**n)))/(3*b*n) - 2*I*elliptic_f(I*(a + b*log(c*x**n))/2, 2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_254():
    f = sqrt(cosh(a + b*log(c*x**n)))/x
    F = -2*I*elliptic_e(I*(a + b*log(c*x**n))/2, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_255():
    f = 1/(x*sqrt(cosh(a + b*log(c*x**n))))
    F = -2*I*elliptic_f(I*(a + b*log(c*x**n))/2, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_256():
    f = 1/(x*cosh(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = 2*sinh(a + b*log(c*x**n))/(b*n*sqrt(cosh(a + b*log(c*x**n)))) + 2*I*elliptic_e(I*(a + b*log(c*x**n))/2, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_257():
    f = 1/(x*cosh(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 2*sinh(a + b*log(c*x**n))/(3*b*n*cosh(a + b*log(c*x**n))**(sympy.S(3)/2)) - 2*I*elliptic_f(I*(a + b*log(c*x**n))/2, 2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_258():
    f = cosh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)
    F = -x*cosh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)/4 + 5*x*cosh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)/(12 + 12*exp(-2*a)/(c*x**n)**(4/n)) + 5*x*exp(-2*a)*cosh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)/(4*(c*x**n)**(4/n)*(1 + exp(-2*a)/(c*x**n)**(4/n))**2) - 5*x*exp(-3*a)*cosh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)*acsch((c*x**n)**(2/n)*exp(a))/(4*(c*x**n)**(6/n)*(1 + exp(-2*a)/(c*x**n)**(4/n))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_259():
    f = sqrt(cosh(a + 2*log(c*x**n)/n))
    F = x*sqrt(cosh(a + 2*log(c*x**n)/n))/2 - x*exp(-a)*sqrt(cosh(a + 2*log(c*x**n)/n))*acsch((c*x**n)**(2/n)*exp(a))/(2*(c*x**n)**(2/n)*sqrt(1 + exp(-2*a)/(c*x**n)**(4/n)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_260():
    f = cosh(a + 2*log(c*x**n)/n)**(sympy.S(-3)/2)
    F = -x*(1 + exp(-2*a)/(c*x**n)**(4/n))/(2*cosh(a + 2*log(c*x**n)/n)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_261():
    f = cosh(a + 2*log(c*x**n)/n)**(sympy.S(-7)/2)
    F = -x*(1 + exp(-2*a)/(c*x**n)**(4/n))/(6*cosh(a + 2*log(c*x**n)/n)**(sympy.S(7)/2)) - x*(1 + exp(-2*a)/(c*x**n)**(4/n))*exp(-2*a)/(15*(c*x**n)**(4/n)*cosh(a + 2*log(c*x**n)/n)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_262():
    f = cosh((a + b*x)/(c + d*x))
    F = (((Symbol('c') + (Symbol('d') * x)) * sympy.cosh(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sinh((Symbol('b') * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cosh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_263():
    f = cosh((a + b*x)/(c + d*x))**2
    F = (((Symbol('c') + (Symbol('d') * x)) * (sympy.cosh(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('CoshIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sinh(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cosh(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_264():
    f = exp(a + b*x)*cosh(a + b*x)**4
    F = -exp(-3*a - 3*b*x)/(48*b) - exp(-a - b*x)/(4*b) + 3*exp(a + b*x)/(8*b) + exp(3*a + 3*b*x)/(12*b) + exp(5*a + 5*b*x)/(80*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_265():
    f = exp(a + b*x)*cosh(a + b*x)**3
    F = 3*x/8 - exp(-2*a - 2*b*x)/(16*b) + 3*exp(2*a + 2*b*x)/(16*b) + exp(4*a + 4*b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_266():
    f = exp(a + b*x)*cosh(a + b*x)**2
    F = -exp(-a - b*x)/(4*b) + exp(a + b*x)/(2*b) + exp(3*a + 3*b*x)/(12*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_267():
    f = exp(a + b*x)*cosh(a + b*x)
    F = x/2 + exp(2*a + 2*b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_268():
    f = exp(a + b*x)*sech(a + b*x)
    F = log(exp(2*a + 2*b*x) + 1)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_269():
    f = exp(a + b*x)*sech(a + b*x)**2
    F = 2*atan(exp(a + b*x))/b - 2*exp(a + b*x)/(b*(exp(2*a + 2*b*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_270():
    f = exp(a + b*x)*sech(a + b*x)**3
    F = 2*exp(4*a + 4*b*x)/(b*(exp(2*a + 2*b*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_271():
    f = exp(a + b*x)*sech(a + b*x)**4
    F = atan(exp(a + b*x))/b + exp(a + b*x)/(b*(exp(2*a + 2*b*x) + 1)) - 2*exp(a + b*x)/(b*(exp(2*a + 2*b*x) + 1)**2) - 8*exp(3*a + 3*b*x)/(3*b*(exp(2*a + 2*b*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_272():
    f = exp(a + b*x)*sech(a + b*x)**5
    F = -8/(b*(exp(2*a + 2*b*x) + 1)**2) + 32/(3*b*(exp(2*a + 2*b*x) + 1)**3) - 4/(b*(exp(2*a + 2*b*x) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_273():
    f = exp(x)*cosh(2*x)**2
    F = exp(5*x)/20 + exp(x)/2 - exp(-3*x)/12
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_274():
    f = exp(x)*cosh(2*x)
    F = exp(3*x)/6 - exp(-x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_275():
    f = exp(x)*sech(2*x)
    F = sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/4 - sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/4 + sqrt(2)*atan(sqrt(2)*exp(x) - 1)/2 + sqrt(2)*atan(sqrt(2)*exp(x) + 1)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_276():
    f = exp(x)*sech(2*x)**2
    F = -sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/8 + sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/8 + sqrt(2)*atan(sqrt(2)*exp(x) - 1)/4 + sqrt(2)*atan(sqrt(2)*exp(x) + 1)/4 - exp(x)/(exp(4*x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_277():
    f = exp(x)*cosh(3*x)**2
    F = exp(7*x)/28 + exp(x)/2 - exp(-5*x)/20
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_278():
    f = exp(x)*cosh(3*x)
    F = exp(4*x)/8 - exp(-2*x)/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_279():
    f = exp(x)*sech(3*x)
    F = -log(exp(2*x) + 1)/3 + log(exp(4*x) - exp(2*x) + 1)/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*exp(2*x))/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_280():
    f = exp(x)*sech(3*x)**2
    F = -sqrt(3)*log(exp(2*x) - sqrt(3)*exp(x) + 1)/18 + sqrt(3)*log(exp(2*x) + sqrt(3)*exp(x) + 1)/18 + atan(2*exp(x) - sqrt(3))/9 + atan(2*exp(x) + sqrt(3))/9 + 2*atan(exp(x))/9 - 2*exp(x)/(3*exp(6*x) + 3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_281():
    f = exp(x)*cosh(4*x)**2
    F = exp(9*x)/36 + exp(x)/2 - exp(-7*x)/28
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_282():
    f = exp(x)*cosh(4*x)
    F = exp(5*x)/10 - exp(-3*x)/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_283():
    f = exp(x)*sech(4*x)
    F = -log(exp(2*x) - sqrt(2 - sqrt(2))*exp(x) + 1)/(4*sqrt(4 - 2*sqrt(2))) + log(exp(2*x) + sqrt(2 - sqrt(2))*exp(x) + 1)/(4*sqrt(4 - 2*sqrt(2))) + log(exp(2*x) - sqrt(sqrt(2) + 2)*exp(x) + 1)/(4*sqrt(2*sqrt(2) + 4)) - log(exp(2*x) + sqrt(sqrt(2) + 2)*exp(x) + 1)/(4*sqrt(2*sqrt(2) + 4)) - atan((-2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(2*sqrt(4 - 2*sqrt(2))) + atan((2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(2*sqrt(4 - 2*sqrt(2))) + atan((-2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(2*sqrt(2*sqrt(2) + 4)) - atan((2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(2*sqrt(2*sqrt(2) + 4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_284():
    f = exp(x)*sech(4*x)**2
    F = -sqrt(2 - sqrt(2))*log(exp(2*x) - sqrt(2 - sqrt(2))*exp(x) + 1)/32 + sqrt(2 - sqrt(2))*log(exp(2*x) + sqrt(2 - sqrt(2))*exp(x) + 1)/32 - sqrt(sqrt(2) + 2)*log(exp(2*x) - sqrt(sqrt(2) + 2)*exp(x) + 1)/32 + sqrt(sqrt(2) + 2)*log(exp(2*x) + sqrt(sqrt(2) + 2)*exp(x) + 1)/32 - atan((-2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(8*sqrt(2*sqrt(2) + 4)) + atan((2*exp(x) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(8*sqrt(2*sqrt(2) + 4)) - atan((-2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(8*sqrt(4 - 2*sqrt(2))) + atan((2*exp(x) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(8*sqrt(4 - 2*sqrt(2))) - exp(x)/(2*exp(8*x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_285():
    f = F**(c*(a + b*x))*cosh(d + e*x)**3
    F = -6*F**(c*(a + b*x))*b*c*e**2*log(F)*cosh(d + e*x)/(b**4*c**4*log(F)**4 - 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) - F**(c*(a + b*x))*b*c*log(F)*cosh(d + e*x)**3/(-b**2*c**2*log(F)**2 + 9*e**2) + 6*F**(c*(a + b*x))*e**3*sinh(d + e*x)/(b**4*c**4*log(F)**4 - 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) + 3*F**(c*(a + b*x))*e*sinh(d + e*x)*cosh(d + e*x)**2/(-b**2*c**2*log(F)**2 + 9*e**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_286():
    f = F**(c*(a + b*x))*cosh(d + e*x)**2
    F = -F**(c*(a + b*x))*b*c*log(F)*cosh(d + e*x)**2/(-b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(c*(a + b*x))*e*sinh(d + e*x)*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(c*(a + b*x))*e**2/(b*c*(-b**2*c**2*log(F)**2 + 4*e**2)*log(F))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_287():
    f = F**(c*(a + b*x))*cosh(d + e*x)
    F = -F**(c*(a + b*x))*b*c*log(F)*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + F**(c*(a + b*x))*e*sinh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_288():
    f = F**(c*(a + b*x))*sech(d + e*x)
    F = 2*F**(c*(a + b*x))*exp(d + e*x)*hyper((1, (b*c*log(F) + e)/(2*e)), (b*c*log(F)/(2*e) + sympy.S(3)/2,), -exp(2*d + 2*e*x))/(b*c*log(F) + e)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_289():
    f = F**(c*(a + b*x))*sech(d + e*x)**2
    F = 4*F**(c*(a + b*x))*exp(2*d + 2*e*x)*hyper((2, b*c*log(F)/(2*e) + 1), (b*c*log(F)/(2*e) + 2,), -exp(2*d + 2*e*x))/(b*c*log(F) + 2*e)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_290():
    f = F**(c*(a + b*x))*sech(d + e*x)**3
    F = F**(c*(a + b*x))*b*c*log(F)*sech(d + e*x)/(2*e**2) + F**(c*(a + b*x))*tanh(d + e*x)*sech(d + e*x)/(2*e) + F**(c*(a + b*x))*(-b*c*log(F) + e)*exp(d + e*x)*hyper((1, (b*c*log(F) + e)/(2*e)), (b*c*log(F)/(2*e) + sympy.S(3)/2,), -exp(2*d + 2*e*x))/e**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_291():
    f = F**(c*(a + b*x))*sech(d + e*x)**4
    F = F**(c*(a + b*x))*b*c*log(F)*sech(d + e*x)**2/(6*e**2) + F**(c*(a + b*x))*tanh(d + e*x)*sech(d + e*x)**2/(3*e) + 2*F**(c*(a + b*x))*(-b*c*log(F) + 2*e)*exp(2*d + 2*e*x)*hyper((2, b*c*log(F)/(2*e) + 1), (b*c*log(F)/(2*e) + 2,), -exp(2*d + 2*e*x))/(3*e**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_292():
    f = (cosh(a*c + b*c*x)**2)**(sympy.S(5)/2)*exp(c*(a + b*x))
    F = 5*x*sqrt(cosh(a*c + b*c*x)**2)*sech(a*c + b*c*x)/16 + sqrt(cosh(a*c + b*c*x)**2)*exp(6*c*(a + b*x))*sech(a*c + b*c*x)/(192*b*c) + 5*sqrt(cosh(a*c + b*c*x)**2)*exp(4*c*(a + b*x))*sech(a*c + b*c*x)/(128*b*c) + 5*sqrt(cosh(a*c + b*c*x)**2)*exp(2*c*(a + b*x))*sech(a*c + b*c*x)/(32*b*c) - 5*sqrt(cosh(a*c + b*c*x)**2)*exp(-2*c*(a + b*x))*sech(a*c + b*c*x)/(64*b*c) - sqrt(cosh(a*c + b*c*x)**2)*exp(-4*c*(a + b*x))*sech(a*c + b*c*x)/(128*b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_293():
    f = (cosh(a*c + b*c*x)**2)**(sympy.S(3)/2)*exp(c*(a + b*x))
    F = 3*x*sqrt(cosh(a*c + b*c*x)**2)*sech(a*c + b*c*x)/8 + sqrt(cosh(a*c + b*c*x)**2)*exp(4*c*(a + b*x))*sech(a*c + b*c*x)/(32*b*c) + 3*sqrt(cosh(a*c + b*c*x)**2)*exp(2*c*(a + b*x))*sech(a*c + b*c*x)/(16*b*c) - sqrt(cosh(a*c + b*c*x)**2)*exp(-2*c*(a + b*x))*sech(a*c + b*c*x)/(16*b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_294():
    f = sqrt(cosh(a*c + b*c*x)**2)*exp(c*(a + b*x))
    F = x*sqrt(cosh(a*c + b*c*x)**2)*sech(a*c + b*c*x)/2 + sqrt(cosh(a*c + b*c*x)**2)*exp(2*c*(a + b*x))*sech(a*c + b*c*x)/(4*b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_295():
    f = exp(c*(a + b*x))/sqrt(cosh(a*c + b*c*x)**2)
    F = log(exp(2*c*(a + b*x)) + 1)*cosh(a*c + b*c*x)/(b*c*sqrt(cosh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_296():
    f = exp(c*(a + b*x))/(cosh(a*c + b*c*x)**2)**(sympy.S(3)/2)
    F = 2*exp(4*c*(a + b*x))*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**2*sqrt(cosh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_297():
    f = exp(c*(a + b*x))/(cosh(a*c + b*c*x)**2)**(sympy.S(5)/2)
    F = -8*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**2*sqrt(cosh(a*c + b*c*x)**2)) + 32*cosh(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**3*sqrt(cosh(a*c + b*c*x)**2)) - 4*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**4*sqrt(cosh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_298():
    f = exp(c*(a + b*x))/(cosh(a*c + b*c*x)**2)**(sympy.S(7)/2)
    F = -64*cosh(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**3*sqrt(cosh(a*c + b*c*x)**2)) + 48*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**4*sqrt(cosh(a*c + b*c*x)**2)) - 192*cosh(a*c + b*c*x)/(5*b*c*(exp(2*c*(a + b*x)) + 1)**5*sqrt(cosh(a*c + b*c*x)**2)) + 32*cosh(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**6*sqrt(cosh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_299():
    f = exp(x)*cosh(a + b*x)
    F = -b*exp(x)*sinh(a + b*x)/(1 - b**2) + exp(x)*cosh(a + b*x)/(1 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_300():
    f = exp(x)*cosh(a + c*x**2)
    F = (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Integer(4) * Symbol('c')))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_301():
    f = exp(x)*cosh(a + b*x + c*x**2)
    F = (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (((Integer(1) + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(1) + (Integer(-1) * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * (((Integer(1) + Symbol('b')))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(1) + Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_302():
    f = exp(x**2)*cosh(a + b*x)
    F = ((Integer(4))**(Integer(-1)) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * ((Integer(-1) * Symbol('b')) + (Integer(2) * x))))) + ((Integer(4))**(Integer(-1)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (Symbol('b') + (Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_303():
    f = exp(x**2)*cosh(a + c*x**2)
    F = ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c')))) * x))) * (((sympy.E)**(Symbol('a')) * (Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c')))))))**(Integer(-1))) + (((sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + Symbol('c'))) * x))) * ((Integer(4) * sympy.sqrt((Integer(1) + Symbol('c')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_304():
    f = exp(x**2)*cosh(a + b*x + c*x**2)
    F = (Integer(-1) * (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * (Integer(1) + (Integer(-1) * Symbol('c')))))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * (Integer(1) + (Integer(-1) * Symbol('c'))) * x))) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c'))))))**(Integer(-1)))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * (Integer(1) + Symbol('c'))))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * (Integer(1) + Symbol('c')) * x)) * ((Integer(2) * sympy.sqrt((Integer(1) + Symbol('c')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(1) + Symbol('c')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_305():
    f = f**(a + b*x)*cosh(d + f*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_306():
    f = f**(a + b*x)*cosh(d + f*x**2)**2
    F = (((sympy.E)**(((Integer(-2) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(8) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(4) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(8) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(4) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_307():
    f = f**(a + b*x)*cosh(d + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(12) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(6) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(12) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(6) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_308():
    f = f**(a + b*x)*cosh(d + e*x + f*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_309():
    f = f**(a + b*x)*cosh(d + e*x + f*x**2)**2
    F = (((sympy.E)**(((Integer(-2) * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(8) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(2) * Symbol('e')) + (Integer(4) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(8) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Integer(4) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_310():
    f = f**(a + b*x)*cosh(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('d')) + ((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(12) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(3) * Symbol('e')) + (Integer(6) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(12) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(6) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_311():
    f = f**(a + c*x**2)*cosh(d + e*x)
    F = ((Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_312():
    f = f**(a + c*x**2)*cosh(d + e*x)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Symbol('c') * x * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('c') * x * sympy.log(Symbol('f')))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_313():
    f = f**(a + c*x**2)*cosh(d + e*x)**3
    F = ((Integer(-3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_314():
    f = f**(a + c*x**2)*cosh(d + f*x**2)
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(4) * (sympy.E)**(Symbol('d')) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(Symbol('d')) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_315():
    f = f**(a + c*x**2)*cosh(d + f*x**2)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(8) * (sympy.E)**((Integer(2) * Symbol('d'))) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Integer(2) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_316():
    f = f**(a + c*x**2)*cosh(d + f*x**2)**3
    F = ((Integer(3) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(16) * (sympy.E)**(Symbol('d')) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(16) * (sympy.E)**((Integer(3) * Symbol('d'))) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**(Symbol('d')) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**((Integer(3) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_317():
    f = f**(a + c*x**2)*cosh(d + e*x + f*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_318():
    f = f**(a + c*x**2)*cosh(d + e*x + f*x**2)**2
    F = (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f')))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (x * ((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * (sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (x * ((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * (sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_319():
    f = f**(a + c*x**2)*cosh(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * (((Integer(12) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * Symbol('e')) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_320():
    f = f**(a + b*x + c*x**2)*cosh(d + e*x)
    F = ((Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_321():
    f = f**(a + b*x + c*x**2)*cosh(d + e*x)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_322():
    f = f**(a + b*x + c*x**2)*cosh(d + e*x)**3
    F = ((Integer(-3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_323():
    f = f**(a + b*x + c*x**2)*cosh(d + f*x**2)
    F = ((Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_324():
    f = f**(a + b*x + c*x**2)*cosh(d + f*x**2)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_325():
    f = f**(a + b*x + c*x**2)*cosh(d + f*x**2)**3
    F = ((Integer(-3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(12) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_326():
    f = f**(a + b*x + c*x**2)*cosh(d + e*x + f*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_327():
    f = f**(a + b*x + c*x**2)*cosh(d + e*x + f*x**2)**2
    F = (((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(8) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(8) * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_328():
    f = f**(a + b*x + c*x**2)*cosh(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('d')) + ((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(12) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_329():
    f = x*sqrt(cosh(x)) + x/cosh(x)**(sympy.S(3)/2)
    F = 2*x*sinh(x)/sqrt(cosh(x)) - 4*sqrt(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_330():
    f = -x/(3*sqrt(cosh(x))) + x/cosh(x)**(sympy.S(5)/2)
    F = 2*x*sinh(x)/(3*cosh(x)**(sympy.S(3)/2)) + 4/(3*sqrt(cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_331():
    f = 3*x*sqrt(cosh(x))/5 + x/cosh(x)**(sympy.S(7)/2)
    F = 6*x*sinh(x)/(5*sqrt(cosh(x))) + 2*x*sinh(x)/(5*cosh(x)**(sympy.S(5)/2)) - 12*sqrt(cosh(x))/5 + 4/(15*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_332():
    f = x**2*sqrt(cosh(x)) + x**2/cosh(x)**(sympy.S(3)/2)
    F = 2*x**2*sinh(x)/sqrt(cosh(x)) - 8*x*sqrt(cosh(x)) - 16*I*elliptic_e(I*x/2, 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_333():
    f = (x + cosh(x))**2
    F = x**3/3 + 2*x*sinh(x) + x/2 + sinh(x)*cosh(x)/2 - 2*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_334():
    f = (x + cosh(x))**3
    F = x**4/4 + 3*x**2*sinh(x) + 3*x**2/4 + 3*x*sinh(x)*cosh(x)/2 - 6*x*cosh(x) + sinh(x)**3/3 + 7*sinh(x) - 3*cosh(x)**2/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_335():
    f = cosh(a + b*x)/(c + d*x**2)
    F = ((sympy.cosh((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('CoshIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('SinhIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_2_Hyperbolic_cosine_6_2_5_Hyperbolic_cosine_functions_336():
    f = cosh(a + b*x)/(c + d*x + e*x**2)
    F = ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('CoshIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1)))) + ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))))
    assert integrate(f, x) == F

