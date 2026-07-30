"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.1 Hyperbolic sine/6.1.5 Hyperbolic sine functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

A, B, C, F, a, b, c, d, e, f, m, n = symbols('A B C F a b c d e f m n')

def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_1():
    f = sinh(a + b*x)
    F = cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_2():
    f = sinh(a + b*x)**2
    F = -x/2 + sinh(a + b*x)*cosh(a + b*x)/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_3():
    f = sinh(a + b*x)**3
    F = cosh(a + b*x)**3/(3*b) - cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_4():
    f = sinh(a + b*x)**4
    F = 3*x/8 + sinh(a + b*x)**3*cosh(a + b*x)/(4*b) - 3*sinh(a + b*x)*cosh(a + b*x)/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_5():
    f = sinh(a + b*x)**5
    F = cosh(a + b*x)**5/(5*b) - 2*cosh(a + b*x)**3/(3*b) + cosh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_6():
    f = sinh(a + b*x)**6
    F = -5*x/16 + sinh(a + b*x)**5*cosh(a + b*x)/(6*b) - 5*sinh(a + b*x)**3*cosh(a + b*x)/(24*b) + 5*sinh(a + b*x)*cosh(a + b*x)/(16*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_7():
    f = sinh(a + b*x)**(sympy.S(7)/2)
    F = -10*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(21*b*sqrt(sinh(a + b*x))) + 2*sinh(a + b*x)**(sympy.S(5)/2)*cosh(a + b*x)/(7*b) - 10*sqrt(sinh(a + b*x))*cosh(a + b*x)/(21*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_8():
    f = sinh(a + b*x)**(sympy.S(5)/2)
    F = 2*sinh(a + b*x)**(sympy.S(3)/2)*cosh(a + b*x)/(5*b) + 6*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(5*b*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_9():
    f = sinh(a + b*x)**(sympy.S(3)/2)
    F = 2*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(3*b*sqrt(sinh(a + b*x))) + 2*sqrt(sinh(a + b*x))*cosh(a + b*x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_10():
    f = sqrt(sinh(a + b*x))
    F = -2*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(b*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_11():
    f = 1/sqrt(sinh(a + b*x))
    F = -2*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(b*sqrt(sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_12():
    f = sinh(a + b*x)**(sympy.S(-3)/2)
    F = -2*cosh(a + b*x)/(b*sqrt(sinh(a + b*x))) - 2*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(b*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_13():
    f = sinh(a + b*x)**(sympy.S(-5)/2)
    F = 2*I*sqrt(I*sinh(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(3*b*sqrt(sinh(a + b*x))) - 2*cosh(a + b*x)/(3*b*sinh(a + b*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_14():
    f = sinh(a + b*x)**(sympy.S(-7)/2)
    F = 6*cosh(a + b*x)/(5*b*sqrt(sinh(a + b*x))) - 2*cosh(a + b*x)/(5*b*sinh(a + b*x)**(sympy.S(5)/2)) + 6*I*sqrt(sinh(a + b*x))*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(5*b*sqrt(I*sinh(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_15():
    f = (b*sinh(c + d*x))**(sympy.S(7)/2)
    F = -10*I*b**4*sqrt(I*sinh(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(21*d*sqrt(b*sinh(c + d*x))) - 10*b**3*sqrt(b*sinh(c + d*x))*cosh(c + d*x)/(21*d) + 2*b*(b*sinh(c + d*x))**(sympy.S(5)/2)*cosh(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_16():
    f = (b*sinh(c + d*x))**(sympy.S(5)/2)
    F = 6*I*b**2*sqrt(b*sinh(c + d*x))*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(5*d*sqrt(I*sinh(c + d*x))) + 2*b*(b*sinh(c + d*x))**(sympy.S(3)/2)*cosh(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_17():
    f = (b*sinh(c + d*x))**(sympy.S(3)/2)
    F = 2*I*b**2*sqrt(I*sinh(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(3*d*sqrt(b*sinh(c + d*x))) + 2*b*sqrt(b*sinh(c + d*x))*cosh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_18():
    f = sqrt(b*sinh(c + d*x))
    F = -2*I*sqrt(b*sinh(c + d*x))*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(d*sqrt(I*sinh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_19():
    f = 1/sqrt(b*sinh(c + d*x))
    F = -2*I*sqrt(I*sinh(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(d*sqrt(b*sinh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_20():
    f = (b*sinh(c + d*x))**(sympy.S(-3)/2)
    F = -2*cosh(c + d*x)/(b*d*sqrt(b*sinh(c + d*x))) - 2*I*sqrt(b*sinh(c + d*x))*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(b**2*d*sqrt(I*sinh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_21():
    f = (b*sinh(c + d*x))**(sympy.S(-5)/2)
    F = -2*cosh(c + d*x)/(3*b*d*(b*sinh(c + d*x))**(sympy.S(3)/2)) + 2*I*sqrt(I*sinh(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(3*b**2*d*sqrt(b*sinh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_22():
    f = (b*sinh(c + d*x))**(sympy.S(-7)/2)
    F = -2*cosh(c + d*x)/(5*b*d*(b*sinh(c + d*x))**(sympy.S(5)/2)) + 6*cosh(c + d*x)/(5*b**3*d*sqrt(b*sinh(c + d*x))) + 6*I*sqrt(b*sinh(c + d*x))*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(5*b**4*d*sqrt(I*sinh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_23():
    f = (I*sinh(c + d*x))**(sympy.S(7)/2)
    F = 2*I*(I*sinh(c + d*x))**(sympy.S(5)/2)*cosh(c + d*x)/(7*d) + 10*I*sqrt(I*sinh(c + d*x))*cosh(c + d*x)/(21*d) - 10*I*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_24():
    f = (I*sinh(c + d*x))**(sympy.S(5)/2)
    F = 2*I*(I*sinh(c + d*x))**(sympy.S(3)/2)*cosh(c + d*x)/(5*d) - 6*I*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_25():
    f = (I*sinh(c + d*x))**(sympy.S(3)/2)
    F = 2*I*sqrt(I*sinh(c + d*x))*cosh(c + d*x)/(3*d) - 2*I*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_26():
    f = sqrt(I*sinh(c + d*x))
    F = -2*I*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_27():
    f = 1/sqrt(I*sinh(c + d*x))
    F = -2*I*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_28():
    f = (I*sinh(c + d*x))**(sympy.S(-3)/2)
    F = 2*I*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/d + 2*I*cosh(c + d*x)/(d*sqrt(I*sinh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_29():
    f = (I*sinh(c + d*x))**(sympy.S(-5)/2)
    F = -2*I*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(3*d) + 2*I*cosh(c + d*x)/(3*d*(I*sinh(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_30():
    f = (I*sinh(c + d*x))**(sympy.S(-7)/2)
    F = 6*I*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(5*d) + 6*I*cosh(c + d*x)/(5*d*sqrt(I*sinh(c + d*x))) + 2*I*cosh(c + d*x)/(5*d*(I*sinh(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_31():
    f = (b*sinh(c + d*x))**(sympy.S(4)/3)
    F = 3*(b*sinh(c + d*x))**(sympy.S(7)/3)*cosh(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), -sinh(c + d*x)**2)/(7*b*d*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_32():
    f = (b*sinh(c + d*x))**(sympy.S(2)/3)
    F = 3*(b*sinh(c + d*x))**(sympy.S(5)/3)*cosh(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), -sinh(c + d*x)**2)/(5*b*d*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_33():
    f = (b*sinh(c + d*x))**(sympy.S(1)/3)
    F = 3*(b*sinh(c + d*x))**(sympy.S(4)/3)*cosh(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), -sinh(c + d*x)**2)/(4*b*d*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_34():
    f = (b*sinh(c + d*x))**(sympy.S(-1)/3)
    F = 3*(b*sinh(c + d*x))**(sympy.S(2)/3)*cosh(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), -sinh(c + d*x)**2)/(2*b*d*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_35():
    f = (b*sinh(c + d*x))**(sympy.S(-2)/3)
    F = 3*(b*sinh(c + d*x))**(sympy.S(1)/3)*cosh(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), -sinh(c + d*x)**2)/(b*d*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_36():
    f = (b*sinh(c + d*x))**(sympy.S(-4)/3)
    F = -3*cosh(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), -sinh(c + d*x)**2)/(b*d*(b*sinh(c + d*x))**(sympy.S(1)/3)*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_37():
    f = (b*sinh(c + d*x))**n
    F = (b*sinh(c + d*x))**(n + 1)*cosh(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -sinh(c + d*x)**2)/(b*d*(n + 1)*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_38():
    f = (I*sinh(c + d*x))**n
    F = -I*(I*sinh(c + d*x))**(n + 1)*cosh(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -sinh(c + d*x)**2)/(d*(n + 1)*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_39():
    f = (-I*sinh(c + d*x))**n
    F = I*(-I*sinh(c + d*x))**(n + 1)*cosh(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -sinh(c + d*x)**2)/(d*(n + 1)*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_40():
    f = sinh(x)**4/(sinh(x) + I)
    F = 3*I*x/2 - 3*I*sinh(x)*cosh(x)/2 + 4*cosh(x)**3/3 - 4*cosh(x) - sinh(x)**3*cosh(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_41():
    f = sinh(x)**3/(sinh(x) + I)
    F = -3*x/2 + 3*sinh(x)*cosh(x)/2 - 2*I*cosh(x) - sinh(x)**2*cosh(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_42():
    f = sinh(x)**2/(sinh(x) + I)
    F = -I*x + cosh(x) + I*cosh(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_43():
    f = sinh(x)/(sinh(x) + I)
    F = x - cosh(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_44():
    f = csch(x)/(sinh(x) + I)
    F = I*atanh(cosh(x)) + cosh(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_45():
    f = csch(x)**2/(sinh(x) + I)
    F = 2*I*coth(x) - atanh(cosh(x)) + coth(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_46():
    f = csch(x)**3/(sinh(x) + I)
    F = 3*I*coth(x)*csch(x)/2 - 2*coth(x) - 3*I*atanh(cosh(x))/2 + coth(x)*csch(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_47():
    f = csch(x)**4/(sinh(x) + I)
    F = 4*I*coth(x)**3/3 - 3*coth(x)*csch(x)/2 - 4*I*coth(x) + 3*atanh(cosh(x))/2 + coth(x)*csch(x)**2/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_48():
    f = sinh(x)**4/(sinh(x) + I)**2
    F = -7*x/2 + 7*sinh(x)*cosh(x)/2 - 16*I*cosh(x)/3 - 8*sinh(x)**2*cosh(x)/(3*sinh(x) + 3*I) - sinh(x)**3*cosh(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_49():
    f = sinh(x)**3/(sinh(x) + I)**2
    F = -2*I*x + 4*cosh(x)/3 + 2*I*cosh(x)/(sinh(x) + I) - sinh(x)**2*cosh(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_50():
    f = sinh(x)**2/(sinh(x) + I)**2
    F = x - 5*cosh(x)/(3*sinh(x) + 3*I) + I*cosh(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_51():
    f = sinh(x)/(sinh(x) + I)**2
    F = -2*I*cosh(x)/(3*sinh(x) + 3*I) - cosh(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_52():
    f = csch(x)/(sinh(x) + I)**2
    F = atanh(cosh(x)) - 4*I*cosh(x)/(3*sinh(x) + 3*I) + cosh(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_53():
    f = csch(x)**2/(sinh(x) + I)**2
    F = 10*coth(x)/3 + 2*I*atanh(cosh(x)) - 2*I*coth(x)/(sinh(x) + I) + coth(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_54():
    f = csch(x)**3/(sinh(x) + I)**2
    F = 7*coth(x)*csch(x)/2 + 16*I*coth(x)/3 - 7*atanh(cosh(x))/2 - 8*I*coth(x)*csch(x)/(3*sinh(x) + 3*I) + coth(x)*csch(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_55():
    f = csch(x)**4/(sinh(x) + I)**2
    F = 4*coth(x)**3 + 5*I*coth(x)*csch(x) - 12*coth(x) - 5*I*atanh(cosh(x)) - 10*I*coth(x)*csch(x)**2/(3*sinh(x) + 3*I) + coth(x)*csch(x)**2/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_56():
    f = 1/(I*sinh(c + d*x) + 1)
    F = I*cosh(c + d*x)/(d*(I*sinh(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_57():
    f = (I*sinh(c + d*x) + 1)**(-2)
    F = I*cosh(c + d*x)/(3*d*(I*sinh(c + d*x) + 1)) + I*cosh(c + d*x)/(3*d*(I*sinh(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_58():
    f = (I*sinh(c + d*x) + 1)**(-3)
    F = 2*I*cosh(c + d*x)/(15*d*(I*sinh(c + d*x) + 1)) + 2*I*cosh(c + d*x)/(15*d*(I*sinh(c + d*x) + 1)**2) + I*cosh(c + d*x)/(5*d*(I*sinh(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_59():
    f = (I*sinh(c + d*x) + 1)**(-4)
    F = 2*I*cosh(c + d*x)/(35*d*(I*sinh(c + d*x) + 1)) + 2*I*cosh(c + d*x)/(35*d*(I*sinh(c + d*x) + 1)**2) + 3*I*cosh(c + d*x)/(35*d*(I*sinh(c + d*x) + 1)**3) + I*cosh(c + d*x)/(7*d*(I*sinh(c + d*x) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_60():
    f = 1/(-I*sinh(c + d*x) + 1)
    F = -I*cosh(c + d*x)/(d*(-I*sinh(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_61():
    f = (-I*sinh(c + d*x) + 1)**(-2)
    F = -I*cosh(c + d*x)/(3*d*(-I*sinh(c + d*x) + 1)) - I*cosh(c + d*x)/(3*d*(-I*sinh(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_62():
    f = (-I*sinh(c + d*x) + 1)**(-3)
    F = -2*I*cosh(c + d*x)/(15*d*(-I*sinh(c + d*x) + 1)) - 2*I*cosh(c + d*x)/(15*d*(-I*sinh(c + d*x) + 1)**2) - I*cosh(c + d*x)/(5*d*(-I*sinh(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_63():
    f = (-I*sinh(c + d*x) + 1)**(-4)
    F = -2*I*cosh(c + d*x)/(35*d*(-I*sinh(c + d*x) + 1)) - 2*I*cosh(c + d*x)/(35*d*(-I*sinh(c + d*x) + 1)**2) - 3*I*cosh(c + d*x)/(35*d*(-I*sinh(c + d*x) + 1)**3) - I*cosh(c + d*x)/(7*d*(-I*sinh(c + d*x) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_64():
    f = sinh(x)/sqrt(I*a*sinh(x) + a)
    F = 2*cosh(x)/sqrt(I*a*sinh(x) + a) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cosh(x)/(2*sqrt(I*a*sinh(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_65():
    f = sinh(x)/sqrt(-I*a*sinh(x) + a)
    F = 2*cosh(x)/sqrt(-I*a*sinh(x) + a) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cosh(x)/(2*sqrt(-I*a*sinh(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_66():
    f = (I*a*sinh(c + d*x) + a)**(sympy.S(5)/2)
    F = 64*I*a**3*cosh(c + d*x)/(15*d*sqrt(I*a*sinh(c + d*x) + a)) + 16*I*a**2*sqrt(I*a*sinh(c + d*x) + a)*cosh(c + d*x)/(15*d) + 2*I*a*(I*a*sinh(c + d*x) + a)**(sympy.S(3)/2)*cosh(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_67():
    f = (I*a*sinh(c + d*x) + a)**(sympy.S(3)/2)
    F = 8*I*a**2*cosh(c + d*x)/(3*d*sqrt(I*a*sinh(c + d*x) + a)) + 2*I*a*sqrt(I*a*sinh(c + d*x) + a)*cosh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_68():
    f = sqrt(I*a*sinh(c + d*x) + a)
    F = 2*I*a*cosh(c + d*x)/(d*sqrt(I*a*sinh(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_69():
    f = 1/sqrt(I*a*sinh(c + d*x) + a)
    F = sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*cosh(c + d*x)/(2*sqrt(I*a*sinh(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_70():
    f = (I*a*sinh(c + d*x) + a)**(sympy.S(-3)/2)
    F = I*cosh(c + d*x)/(2*d*(I*a*sinh(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*cosh(c + d*x)/(2*sqrt(I*a*sinh(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_71():
    f = (I*a*sinh(c + d*x) + a)**(sympy.S(-5)/2)
    F = I*cosh(c + d*x)/(4*d*(I*a*sinh(c + d*x) + a)**(sympy.S(5)/2)) + 3*I*cosh(c + d*x)/(16*a*d*(I*a*sinh(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*cosh(c + d*x)/(2*sqrt(I*a*sinh(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_72():
    f = sinh(x)**4/(a + b*sinh(x))
    F = -2*a**4*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**4*sqrt(a**2 + b**2)) - a*sinh(x)*cosh(x)/(2*b**2) - a*x*(2*a**2 - b**2)/(2*b**4) - (-3*a**2/b**2 + 2)*cosh(x)/(3*b) + sinh(x)**2*cosh(x)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_73():
    f = sinh(x)**3/(a + b*sinh(x))
    F = 2*a**3*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**3*sqrt(a**2 + b**2)) - a*cosh(x)/b**2 + sinh(x)*cosh(x)/(2*b) + x*(2*a**2 - b**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_74():
    f = sinh(x)**2/(a + b*sinh(x))
    F = -2*a**2*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**2*sqrt(a**2 + b**2)) - a*x/b**2 + cosh(x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_75():
    f = sinh(x)/(a + b*sinh(x))
    F = 2*a*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b*sqrt(a**2 + b**2)) + x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_76():
    f = csch(x)/(a + b*sinh(x))
    F = 2*b*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a*sqrt(a**2 + b**2)) - atanh(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_77():
    f = csch(x)**2/(a + b*sinh(x))
    F = -coth(x)/a - 2*b**2*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2*sqrt(a**2 + b**2)) + b*atanh(cosh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_78():
    f = csch(x)**3/(a + b*sinh(x))
    F = -coth(x)*csch(x)/(2*a) + b*coth(x)/a**2 + 2*b**3*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**3*sqrt(a**2 + b**2)) + (a**2 - 2*b**2)*atanh(cosh(x))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_79():
    f = csch(x)**4/(a + b*sinh(x))
    F = -coth(x)*csch(x)**2/(3*a) + b*coth(x)*csch(x)/(2*a**2) + (2*a**2 - 3*b**2)*coth(x)/(3*a**3) - 2*b**4*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**4*sqrt(a**2 + b**2)) - b*(a**2 - 2*b**2)*atanh(cosh(x))/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_80():
    f = sinh(x)**4/(a + b*sinh(x))**2
    F = 2*a**3*(3*a**2 + 4*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**4*(a**2 + b**2)**(sympy.S(3)/2)) - a**2*sinh(x)**2*cosh(x)/(b*(a + b*sinh(x))*(a**2 + b**2)) - a*(3*a**2 + 2*b**2)*cosh(x)/(b**3*(a**2 + b**2)) + (3*a**2 + b**2)*sinh(x)*cosh(x)/(2*b**2*(a**2 + b**2)) + x*(6*a**2 - b**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_81():
    f = sinh(x)**3/(a + b*sinh(x))**2
    F = -a**2*sinh(x)*cosh(x)/(b*(a + b*sinh(x))*(a**2 + b**2)) - 2*a**2*(2*a**2 + 3*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**3*(a**2 + b**2)**(sympy.S(3)/2)) - 2*a*x/b**3 + (2*a**2 + b**2)*cosh(x)/(b**2*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_82():
    f = sinh(x)**2/(a + b*sinh(x))**2
    F = -a**2*cosh(x)/(b*(a + b*sinh(x))*(a**2 + b**2)) + 2*a*(a**2 + 2*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**2*(a**2 + b**2)**(sympy.S(3)/2)) + x/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_83():
    f = sinh(x)/(a + b*sinh(x))**2
    F = a*cosh(x)/((a + b*sinh(x))*(a**2 + b**2)) - 2*b*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_84():
    f = csch(x)/(a + b*sinh(x))**2
    F = b**2*cosh(x)/(a*(a + b*sinh(x))*(a**2 + b**2)) + 2*b*(2*a**2 + b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2*(a**2 + b**2)**(sympy.S(3)/2)) - atanh(cosh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_85():
    f = csch(x)**2/(a + b*sinh(x))**2
    F = b**2*coth(x)/(a*(a + b*sinh(x))*(a**2 + b**2)) - (a**2 + 2*b**2)*coth(x)/(a**2*(a**2 + b**2)) - 2*b**2*(3*a**2 + 2*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**3*(a**2 + b**2)**(sympy.S(3)/2)) + 2*b*atanh(cosh(x))/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_86():
    f = csch(x)**3/(a + b*sinh(x))**2
    F = b**2*coth(x)*csch(x)/(a*(a + b*sinh(x))*(a**2 + b**2)) - (a**2 + 3*b**2)*coth(x)*csch(x)/(2*a**2*(a**2 + b**2)) + b*(2*a**2 + 3*b**2)*coth(x)/(a**3*(a**2 + b**2)) + 2*b**3*(4*a**2 + 3*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**4*(a**2 + b**2)**(sympy.S(3)/2)) + (a**2 - 6*b**2)*atanh(cosh(x))/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_87():
    f = csch(x)**4/(a + b*sinh(x))**2
    F = b**2*coth(x)*csch(x)**2/(a*(a + b*sinh(x))*(a**2 + b**2)) - (a**2 + 4*b**2)*coth(x)*csch(x)**2/(3*a**2*(a**2 + b**2)) + b*(a**2 + 2*b**2)*coth(x)*csch(x)/(a**3*(a**2 + b**2)) + (2*a**4 - 7*a**2*b**2 - 12*b**4)*coth(x)/(3*a**4*(a**2 + b**2)) - 2*b**4*(5*a**2 + 4*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**5*(a**2 + b**2)**(sympy.S(3)/2)) - b*(a**2 - 4*b**2)*atanh(cosh(x))/a**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_88():
    f = 1/(5*I*sinh(c + d*x) + 3)
    F = I*log(I*sinh(c/2 + d*x/2) + 3*cosh(c/2 + d*x/2))/(4*d) - I*log(3*I*sinh(c/2 + d*x/2) + cosh(c/2 + d*x/2))/(4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_89():
    f = (5*I*sinh(c + d*x) + 3)**(-2)
    F = -3*I*log(I*sinh(c/2 + d*x/2) + 3*cosh(c/2 + d*x/2))/(64*d) + 3*I*log(3*I*sinh(c/2 + d*x/2) + cosh(c/2 + d*x/2))/(64*d) + 5*I*cosh(c + d*x)/(16*d*(5*I*sinh(c + d*x) + 3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_90():
    f = (5*I*sinh(c + d*x) + 3)**(-3)
    F = 43*I*log(I*sinh(c/2 + d*x/2) + 3*cosh(c/2 + d*x/2))/(2048*d) - 43*I*log(3*I*sinh(c/2 + d*x/2) + cosh(c/2 + d*x/2))/(2048*d) - 45*I*cosh(c + d*x)/(512*d*(5*I*sinh(c + d*x) + 3)) + 5*I*cosh(c + d*x)/(32*d*(5*I*sinh(c + d*x) + 3)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_91():
    f = (5*I*sinh(c + d*x) + 3)**(-4)
    F = -279*I*log(I*sinh(c/2 + d*x/2) + 3*cosh(c/2 + d*x/2))/(32768*d) + 279*I*log(3*I*sinh(c/2 + d*x/2) + cosh(c/2 + d*x/2))/(32768*d) + 995*I*cosh(c + d*x)/(24576*d*(5*I*sinh(c + d*x) + 3)) - 25*I*cosh(c + d*x)/(512*d*(5*I*sinh(c + d*x) + 3)**2) + 5*I*cosh(c + d*x)/(48*d*(5*I*sinh(c + d*x) + 3)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_92():
    f = 1/(3*I*sinh(c + d*x) + 5)
    F = x/4 - I*atan(cosh(c + d*x)/(I*sinh(c + d*x) + 3))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_93():
    f = (3*I*sinh(c + d*x) + 5)**(-2)
    F = 5*x/64 - 5*I*atan(cosh(c + d*x)/(I*sinh(c + d*x) + 3))/(32*d) - 3*I*cosh(c + d*x)/(16*d*(3*I*sinh(c + d*x) + 5))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_94():
    f = (3*I*sinh(c + d*x) + 5)**(-3)
    F = 59*x/2048 - 59*I*atan(cosh(c + d*x)/(I*sinh(c + d*x) + 3))/(1024*d) - 45*I*cosh(c + d*x)/(512*d*(3*I*sinh(c + d*x) + 5)) - 3*I*cosh(c + d*x)/(32*d*(3*I*sinh(c + d*x) + 5)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_95():
    f = (3*I*sinh(c + d*x) + 5)**(-4)
    F = 385*x/32768 - 385*I*atan(cosh(c + d*x)/(I*sinh(c + d*x) + 3))/(16384*d) - 311*I*cosh(c + d*x)/(8192*d*(3*I*sinh(c + d*x) + 5)) - 25*I*cosh(c + d*x)/(512*d*(3*I*sinh(c + d*x) + 5)**2) - I*cosh(c + d*x)/(16*d*(3*I*sinh(c + d*x) + 5)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_96():
    f = (a + b*sinh(c + d*x))**5
    F = 7*a*b**2*(22*a**2 - 23*b**2)*sinh(c + d*x)*cosh(c + d*x)/(120*d) + 9*a*b*(a + b*sinh(c + d*x))**3*cosh(c + d*x)/(20*d) + a*x*(8*a**4 - 40*a**2*b**2 + 15*b**4)/8 + b*(a + b*sinh(c + d*x))**4*cosh(c + d*x)/(5*d) + b*(a + b*sinh(c + d*x))**2*(47*a**2 - 16*b**2)*cosh(c + d*x)/(60*d) + b*(107*a**4 - 192*a**2*b**2 + 16*b**4)*cosh(c + d*x)/(30*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_97():
    f = (a + b*sinh(c + d*x))**4
    F = 7*a*b*(a + b*sinh(c + d*x))**2*cosh(c + d*x)/(12*d) + a*b*(19*a**2 - 16*b**2)*cosh(c + d*x)/(6*d) + b**2*(26*a**2 - 9*b**2)*sinh(c + d*x)*cosh(c + d*x)/(24*d) + b*(a + b*sinh(c + d*x))**3*cosh(c + d*x)/(4*d) + x*(a**4 - 3*a**2*b**2 + 3*b**4/8)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_98():
    f = (a + b*sinh(c + d*x))**3
    F = 5*a*b**2*sinh(c + d*x)*cosh(c + d*x)/(6*d) + a*x*(2*a**2 - 3*b**2)/2 + b*(a + b*sinh(c + d*x))**2*cosh(c + d*x)/(3*d) + 2*b*(4*a**2 - b**2)*cosh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_99():
    f = (a + b*sinh(c + d*x))**2
    F = 2*a*b*cosh(c + d*x)/d + b**2*sinh(c + d*x)*cosh(c + d*x)/(2*d) + x*(a**2 - b**2/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_100():
    f = a + b*sinh(c + d*x)
    F = a*x + b*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_101():
    f = 1/(a + b*sinh(c + d*x))
    F = -2*atanh((-a*tanh(c/2 + d*x/2) + b)/sqrt(a**2 + b**2))/(d*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_102():
    f = (a + b*sinh(c + d*x))**(-2)
    F = -2*a*atanh((-a*tanh(c/2 + d*x/2) + b)/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(3)/2)) - b*cosh(c + d*x)/(d*(a + b*sinh(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_103():
    f = (a + b*sinh(c + d*x))**(-3)
    F = -3*a*b*cosh(c + d*x)/(2*d*(a + b*sinh(c + d*x))*(a**2 + b**2)**2) - b*cosh(c + d*x)/(d*(a + b*sinh(c + d*x))**2*(2*a**2 + 2*b**2)) - (2*a**2 - b**2)*atanh((-a*tanh(c/2 + d*x/2) + b)/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_104():
    f = (a + b*sinh(c + d*x))**(-4)
    F = -5*a*b*cosh(c + d*x)/(6*d*(a + b*sinh(c + d*x))**2*(a**2 + b**2)**2) - a*(2*a**2 - 3*b**2)*atanh((-a*tanh(c/2 + d*x/2) + b)/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(7)/2)) - b*(11*a**2 - 4*b**2)*cosh(c + d*x)/(6*d*(a + b*sinh(c + d*x))*(a**2 + b**2)**3) - b*cosh(c + d*x)/(d*(a + b*sinh(c + d*x))**3*(3*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_105():
    f = (a + b*sinh(x))**(sympy.S(5)/2)
    F = 16*a*b*sqrt(a + b*sinh(x))*cosh(x)/15 + 16*I*a*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2)*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(15*sqrt(a + b*sinh(x))) + 2*b*(a + b*sinh(x))**(sympy.S(3)/2)*cosh(x)/5 - 2*I*sqrt(a + b*sinh(x))*(23*a**2 - 9*b**2)*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(15*sqrt((a + b*sinh(x))/(a - I*b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_106():
    f = (a + b*sinh(x))**(sympy.S(3)/2)
    F = -8*I*a*sqrt(a + b*sinh(x))*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(3*sqrt((a + b*sinh(x))/(a - I*b))) + 2*b*sqrt(a + b*sinh(x))*cosh(x)/3 + 2*I*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2)*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(3*sqrt(a + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_107():
    f = sqrt(a + b*sinh(x))
    F = -2*I*sqrt(a + b*sinh(x))*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/sqrt((a + b*sinh(x))/(a - I*b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_108():
    f = 1/sqrt(a + b*sinh(x))
    F = -2*I*sqrt((a + b*sinh(x))/(a - I*b))*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/sqrt(a + b*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_109():
    f = (a + b*sinh(x))**(sympy.S(-3)/2)
    F = -2*b*cosh(x)/(sqrt(a + b*sinh(x))*(a**2 + b**2)) - 2*I*sqrt(a + b*sinh(x))*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_110():
    f = (a + b*sinh(x))**(sympy.S(-5)/2)
    F = -8*a*b*cosh(x)/(3*sqrt(a + b*sinh(x))*(a**2 + b**2)**2) - 8*I*a*sqrt(a + b*sinh(x))*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(3*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2)**2) - 2*b*cosh(x)/((a + b*sinh(x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + 2*I*sqrt((a + b*sinh(x))/(a - I*b))*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(sqrt(a + b*sinh(x))*(3*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_111():
    f = sinh(x)/sqrt(a + b*sinh(x))
    F = 2*I*a*sqrt((a + b*sinh(x))/(a - I*b))*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(b*sqrt(a + b*sinh(x))) - 2*I*sqrt(a + b*sinh(x))*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(b*sqrt((a + b*sinh(x))/(a - I*b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_112():
    f = (A + B*sinh(x))*(I*a*sinh(x) + a)**(sympy.S(5)/2)
    F = 2*B*(I*a*sinh(x) + a)**(sympy.S(5)/2)*cosh(x)/7 + 64*a**3*(7*I*A + 5*B)*cosh(x)/(105*sqrt(I*a*sinh(x) + a)) + 16*a**2*(7*I*A + 5*B)*sqrt(I*a*sinh(x) + a)*cosh(x)/105 + 2*a*(7*I*A + 5*B)*(I*a*sinh(x) + a)**(sympy.S(3)/2)*cosh(x)/35
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_113():
    f = (A + B*sinh(x))*(I*a*sinh(x) + a)**(sympy.S(3)/2)
    F = 2*B*(I*a*sinh(x) + a)**(sympy.S(3)/2)*cosh(x)/5 + 8*a**2*(5*I*A + 3*B)*cosh(x)/(15*sqrt(I*a*sinh(x) + a)) + 2*a*(5*I*A + 3*B)*sqrt(I*a*sinh(x) + a)*cosh(x)/15
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_114():
    f = (A + B*sinh(x))*sqrt(I*a*sinh(x) + a)
    F = 2*B*sqrt(I*a*sinh(x) + a)*cosh(x)/3 + 2*a*(3*I*A + B)*cosh(x)/(3*sqrt(I*a*sinh(x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_115():
    f = (A + B*sinh(x))/(sinh(x) + I)
    F = B*x - (I*A + B)*cosh(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_116():
    f = (A + B*sinh(x))/(sinh(x) + I)**2
    F = -(A + 2*I*B)*cosh(x)/(3*sinh(x) + 3*I) - (I*A + B)*cosh(x)/(3*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_117():
    f = (A + B*sinh(x))/(sinh(x) + I)**3
    F = -(2*A + 3*I*B)*cosh(x)/(15*(sinh(x) + I)**2) - (I*A + B)*cosh(x)/(5*(sinh(x) + I)**3) + (2*I*A - 3*B)*cosh(x)/(15*sinh(x) + 15*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_118():
    f = (A + B*sinh(x))/(sinh(x) + I)**4
    F = -(3*A + 4*I*B)*cosh(x)/(35*(sinh(x) + I)**3) + (6*A + 8*I*B)*cosh(x)/(105*sinh(x) + 105*I) - (I*A + B)*cosh(x)/(7*(sinh(x) + I)**4) + (6*I*A - 8*B)*cosh(x)/(105*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_119():
    f = (A + B*sinh(x))/(-sinh(x) + I)
    F = -B*x + (I*A - B)*cosh(x)/(-sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_120():
    f = (A + B*sinh(x))/(-sinh(x) + I)**2
    F = (A - 2*I*B)*cosh(x)/(-3*sinh(x) + 3*I) + (I*A - B)*cosh(x)/(3*(-sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_121():
    f = (A + B*sinh(x))/(-sinh(x) + I)**3
    F = (2*A - 3*I*B)*cosh(x)/(15*(-sinh(x) + I)**2) + (I*A - B)*cosh(x)/(5*(-sinh(x) + I)**3) - (2*I*A + 3*B)*cosh(x)/(-15*sinh(x) + 15*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_122():
    f = (A + B*sinh(x))/(-sinh(x) + I)**4
    F = (3*A - 4*I*B)*cosh(x)/(35*(-sinh(x) + I)**3) - (6*A - 8*I*B)*cosh(x)/(-105*sinh(x) + 105*I) + (I*A - B)*cosh(x)/(7*(-sinh(x) + I)**4) - (6*I*A + 8*B)*cosh(x)/(105*(-sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_123():
    f = (A + B*sinh(x))/sqrt(I*a*sinh(x) + a)
    F = 2*B*cosh(x)/sqrt(I*a*sinh(x) + a) + sqrt(2)*(I*A - B)*atanh(sqrt(2)*sqrt(a)*cosh(x)/(2*sqrt(I*a*sinh(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_124():
    f = (A + B*sinh(x))/(I*a*sinh(x) + a)**(sympy.S(3)/2)
    F = (I*A - B)*cosh(x)/(2*(I*a*sinh(x) + a)**(sympy.S(3)/2)) + sqrt(2)*(I*A + 3*B)*atanh(sqrt(2)*sqrt(a)*cosh(x)/(2*sqrt(I*a*sinh(x) + a)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_125():
    f = (A + B*sinh(x))/(I*a*sinh(x) + a)**(sympy.S(5)/2)
    F = (I*A - B)*cosh(x)/(4*(I*a*sinh(x) + a)**(sympy.S(5)/2)) + (3*I*A + 5*B)*cosh(x)/(16*a*(I*a*sinh(x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*I*A + 5*B)*atanh(sqrt(2)*sqrt(a)*cosh(x)/(2*sqrt(I*a*sinh(x) + a)))/(32*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_126():
    f = (A + B*sinh(x))*(a + b*sinh(x))**(sympy.S(5)/2)
    F = 2*B*(a + b*sinh(x))**(sympy.S(5)/2)*cosh(x)/7 + (a + b*sinh(x))**(sympy.S(3)/2)*(2*A*b/5 + 2*B*a/7)*cosh(x) + sqrt(a + b*sinh(x))*(16*A*a*b/15 + 2*B*a**2/7 - 10*B*b**2/21)*cosh(x) + 2*I*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2)*(56*A*a*b + 15*B*a**2 - 25*B*b**2)*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(105*b*sqrt(a + b*sinh(x))) - 2*I*sqrt(a + b*sinh(x))*(161*A*a**2*b - 63*A*b**3 + 15*B*a**3 - 145*B*a*b**2)*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(105*b*sqrt((a + b*sinh(x))/(a - I*b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_127():
    f = (A + B*sinh(x))*(a + b*sinh(x))**(sympy.S(3)/2)
    F = 2*B*(a + b*sinh(x))**(sympy.S(3)/2)*cosh(x)/5 + sqrt(a + b*sinh(x))*(2*A*b/3 + 2*B*a/5)*cosh(x) + 2*I*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2)*(5*A*b + 3*B*a)*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(15*b*sqrt(a + b*sinh(x))) - 2*I*sqrt(a + b*sinh(x))*(20*A*a*b + 3*B*a**2 - 9*B*b**2)*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(15*b*sqrt((a + b*sinh(x))/(a - I*b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_128():
    f = (A + B*sinh(x))*sqrt(a + b*sinh(x))
    F = 2*B*sqrt(a + b*sinh(x))*cosh(x)/3 + 2*I*B*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2)*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(3*b*sqrt(a + b*sinh(x))) - 2*I*sqrt(a + b*sinh(x))*(3*A*b + B*a)*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(3*b*sqrt((a + b*sinh(x))/(a - I*b)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_129():
    f = (A + B*sinh(x))/(a + b*sinh(x))
    F = B*x/b - (2*A*b - 2*B*a)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_130():
    f = (A + B*sinh(x))/(a + b*sinh(x))**2
    F = -(2*A*a + 2*B*b)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) - (A*b - B*a)*cosh(x)/((a + b*sinh(x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_131():
    f = (A + B*sinh(x))/(a + b*sinh(x))**3
    F = -(2*A*a**2 - A*b**2 + 3*B*a*b)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - (3*A*a*b - B*a**2 + 2*B*b**2)*cosh(x)/(2*(a + b*sinh(x))*(a**2 + b**2)**2) - (A*b - B*a)*cosh(x)/((a + b*sinh(x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_132():
    f = (A + B*sinh(x))/(a + b*sinh(x))**4
    F = -(2*A*a**3 - 3*A*a*b**2 + 4*B*a**2*b - B*b**3)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(7)/2) - (11*A*a**2*b - 4*A*b**3 - 2*B*a**3 + 13*B*a*b**2)*cosh(x)/(6*(a + b*sinh(x))*(a**2 + b**2)**3) - (5*A*a*b - 2*B*a**2 + 3*B*b**2)*cosh(x)/(6*(a + b*sinh(x))**2*(a**2 + b**2)**2) - (A*b - B*a)*cosh(x)/((a + b*sinh(x))**3*(3*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_133():
    f = (B*sinh(x) + B*b/a)/(a + b*sinh(x))
    F = B*x/b + B*(2*a**2 - 2*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a*b*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_134():
    f = (B*a/b + B*sinh(x))/(a + b*sinh(x))
    F = B*x/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_135():
    f = (a - b*sinh(x))/(a*sinh(x) + b)**2
    F = -cosh(x)/(a*sinh(x) + b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_136():
    f = (2 - sinh(x))/(sinh(x) + 2)
    F = -x + 4*sqrt(5)*x/5 - 8*sqrt(5)*atanh(cosh(x)/(sinh(x) + 2 + sqrt(5)))/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_137():
    f = (A + B*sinh(x))/sqrt(a + b*sinh(x))
    F = -2*I*B*sqrt(a + b*sinh(x))*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(b*sqrt((a + b*sinh(x))/(a - I*b))) - 2*I*sqrt((a + b*sinh(x))/(a - I*b))*(A*b - B*a)*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(b*sqrt(a + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_138():
    f = (A + B*sinh(x))/(a + b*sinh(x))**(sympy.S(3)/2)
    F = -2*I*B*sqrt((a + b*sinh(x))/(a - I*b))*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(b*sqrt(a + b*sinh(x))) - (2*A*b - 2*B*a)*cosh(x)/(sqrt(a + b*sinh(x))*(a**2 + b**2)) - 2*I*sqrt(a + b*sinh(x))*(A*b - B*a)*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(b*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_139():
    f = (A + B*sinh(x))/(a + b*sinh(x))**(sympy.S(5)/2)
    F = -(8*A*a*b - 2*B*a**2 + 6*B*b**2)*cosh(x)/(3*sqrt(a + b*sinh(x))*(a**2 + b**2)**2) - (2*A*b - 2*B*a)*cosh(x)/((a + b*sinh(x))**(sympy.S(3)/2)*(3*a**2 + 3*b**2)) + 2*I*sqrt((a + b*sinh(x))/(a - I*b))*(A*b - B*a)*elliptic_f(I*x/2 - pi/4, 2*b/(I*a + b))/(3*b*sqrt(a + b*sinh(x))*(a**2 + b**2)) - 2*I*sqrt(a + b*sinh(x))*(4*A*a*b - B*a**2 + 3*B*b**2)*elliptic_e(I*x/2 - pi/4, 2*b/(I*a + b))/(3*b*sqrt((a + b*sinh(x))/(a - I*b))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_140():
    f = (a*sinh(x)**2)**(sympy.S(5)/2)
    F = 8*a**2*sqrt(a*sinh(x)**2)*coth(x)/15 - 4*a*(a*sinh(x)**2)**(sympy.S(3)/2)*coth(x)/15 + (a*sinh(x)**2)**(sympy.S(5)/2)*coth(x)/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_141():
    f = (a*sinh(x)**2)**(sympy.S(3)/2)
    F = -2*a*sqrt(a*sinh(x)**2)*coth(x)/3 + (a*sinh(x)**2)**(sympy.S(3)/2)*coth(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_142():
    f = sqrt(a*sinh(x)**2)
    F = sqrt(a*sinh(x)**2)*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_143():
    f = 1/sqrt(a*sinh(x)**2)
    F = -sinh(x)*atanh(cosh(x))/sqrt(a*sinh(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_144():
    f = (a*sinh(x)**2)**(sympy.S(-3)/2)
    F = sinh(x)*atanh(cosh(x))/(2*a*sqrt(a*sinh(x)**2)) - coth(x)/(2*a*sqrt(a*sinh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_145():
    f = (a*sinh(x)**2)**(sympy.S(-5)/2)
    F = -coth(x)/(4*a*(a*sinh(x)**2)**(sympy.S(3)/2)) - 3*sinh(x)*atanh(cosh(x))/(8*a**2*sqrt(a*sinh(x)**2)) + 3*coth(x)/(8*a**2*sqrt(a*sinh(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_146():
    f = (a*sinh(x)**3)**(sympy.S(5)/2)
    F = -26*I*a**2*sqrt(I*sinh(x))*sqrt(a*sinh(x)**3)*csch(x)**2*elliptic_f(I*x/2 - pi/4, 2)/77 + 2*a**2*sqrt(a*sinh(x)**3)*sinh(x)**5*cosh(x)/15 - 26*a**2*sqrt(a*sinh(x)**3)*sinh(x)**3*cosh(x)/165 + 78*a**2*sqrt(a*sinh(x)**3)*sinh(x)*cosh(x)/385 - 26*a**2*sqrt(a*sinh(x)**3)*coth(x)/77
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_147():
    f = (a*sinh(x)**3)**(sympy.S(3)/2)
    F = 2*a*sqrt(a*sinh(x)**3)*sinh(x)**2*cosh(x)/9 - 14*a*sqrt(a*sinh(x)**3)*cosh(x)/45 - 14*I*a*sqrt(a*sinh(x)**3)*csch(x)*elliptic_e(I*x/2 - pi/4, 2)/(15*sqrt(I*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_148():
    f = sqrt(a*sinh(x)**3)
    F = 2*I*sqrt(I*sinh(x))*sqrt(a*sinh(x)**3)*csch(x)**2*elliptic_f(I*x/2 - pi/4, 2)/3 + 2*sqrt(a*sinh(x)**3)*coth(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_149():
    f = 1/sqrt(a*sinh(x)**3)
    F = -2*sinh(x)*cosh(x)/sqrt(a*sinh(x)**3) - 2*I*sinh(x)**2*elliptic_e(I*x/2 - pi/4, 2)/(sqrt(I*sinh(x))*sqrt(a*sinh(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_150():
    f = (a*sinh(x)**3)**(sympy.S(-3)/2)
    F = -10*I*sqrt(I*sinh(x))*sinh(x)*elliptic_f(I*x/2 - pi/4, 2)/(21*a*sqrt(a*sinh(x)**3)) + 10*cosh(x)/(21*a*sqrt(a*sinh(x)**3)) - 2*coth(x)*csch(x)/(7*a*sqrt(a*sinh(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_151():
    f = (a*sinh(x)**3)**(sympy.S(-5)/2)
    F = 154*sinh(x)*cosh(x)/(195*a**2*sqrt(a*sinh(x)**3)) - 2*coth(x)*csch(x)**4/(13*a**2*sqrt(a*sinh(x)**3)) + 22*coth(x)*csch(x)**2/(117*a**2*sqrt(a*sinh(x)**3)) - 154*coth(x)/(585*a**2*sqrt(a*sinh(x)**3)) + 154*I*sinh(x)**2*elliptic_e(I*x/2 - pi/4, 2)/(195*a**2*sqrt(I*sinh(x))*sqrt(a*sinh(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_152():
    f = (a*sinh(x)**4)**(sympy.S(5)/2)
    F = -63*a**2*x*sqrt(a*sinh(x)**4)*csch(x)**2/256 + a**2*sqrt(a*sinh(x)**4)*sinh(x)**7*cosh(x)/10 - 9*a**2*sqrt(a*sinh(x)**4)*sinh(x)**5*cosh(x)/80 + 21*a**2*sqrt(a*sinh(x)**4)*sinh(x)**3*cosh(x)/160 - 21*a**2*sqrt(a*sinh(x)**4)*sinh(x)*cosh(x)/128 + 63*a**2*sqrt(a*sinh(x)**4)*coth(x)/256
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_153():
    f = (a*sinh(x)**4)**(sympy.S(3)/2)
    F = -5*a*x*sqrt(a*sinh(x)**4)*csch(x)**2/16 + a*sqrt(a*sinh(x)**4)*sinh(x)**3*cosh(x)/6 - 5*a*sqrt(a*sinh(x)**4)*sinh(x)*cosh(x)/24 + 5*a*sqrt(a*sinh(x)**4)*coth(x)/16
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_154():
    f = sqrt(a*sinh(x)**4)
    F = -x*sqrt(a*sinh(x)**4)*csch(x)**2/2 + sqrt(a*sinh(x)**4)*coth(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_155():
    f = 1/sqrt(a*sinh(x)**4)
    F = -sinh(x)*cosh(x)/sqrt(a*sinh(x)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_156():
    f = (a*sinh(x)**4)**(sympy.S(-3)/2)
    F = -sinh(x)*cosh(x)/(a*sqrt(a*sinh(x)**4)) - cosh(x)**2*coth(x)**3/(5*a*sqrt(a*sinh(x)**4)) + 2*cosh(x)**2*coth(x)/(3*a*sqrt(a*sinh(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_157():
    f = (a*sinh(x)**4)**(sympy.S(-5)/2)
    F = -sinh(x)*cosh(x)/(a**2*sqrt(a*sinh(x)**4)) - cosh(x)**2*coth(x)**7/(9*a**2*sqrt(a*sinh(x)**4)) + 4*cosh(x)**2*coth(x)**5/(7*a**2*sqrt(a*sinh(x)**4)) - 6*cosh(x)**2*coth(x)**3/(5*a**2*sqrt(a*sinh(x)**4)) + 4*cosh(x)**2*coth(x)/(3*a**2*sqrt(a*sinh(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_158():
    f = cosh(x)**8/(sinh(x) + I)
    F = -5*I*x/16 - I*sinh(x)*cosh(x)**5/6 - 5*I*sinh(x)*cosh(x)**3/24 - 5*I*sinh(x)*cosh(x)/16 + cosh(x)**7/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_159():
    f = cosh(x)**7/(sinh(x) + I)
    F = (-sinh(x) + I)**6/6 - 4*I*(-sinh(x) + I)**5/5 - (-sinh(x) + I)**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_160():
    f = cosh(x)**6/(sinh(x) + I)
    F = -3*I*x/8 - I*sinh(x)*cosh(x)**3/4 - 3*I*sinh(x)*cosh(x)/8 + cosh(x)**5/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_161():
    f = cosh(x)**5/(sinh(x) + I)
    F = sinh(x)**4/4 - I*sinh(x)**3/3 + sinh(x)**2/2 - I*sinh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_162():
    f = cosh(x)**4/(sinh(x) + I)
    F = -I*x/2 - I*sinh(x)*cosh(x)/2 + cosh(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_163():
    f = cosh(x)**3/(sinh(x) + I)
    F = sinh(x)**2/2 - I*sinh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_164():
    f = cosh(x)**2/(sinh(x) + I)
    F = -I*x + cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_165():
    f = cosh(x)/(sinh(x) + I)
    F = log(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_166():
    f = sech(x)/(sinh(x) + I)
    F = -I*atan(sinh(x))/2 - I/(2*sinh(x) + 2*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_167():
    f = sech(x)**2/(sinh(x) + I)
    F = -2*I*tanh(x)/3 - I*sech(x)/(3*sinh(x) + 3*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_168():
    f = sech(x)**3/(sinh(x) + I)
    F = -3*I*atan(sinh(x))/8 - I/(4*sinh(x) + 4*I) + 1/(8*(sinh(x) + I)**2) + I/(-8*sinh(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_169():
    f = sech(x)**4/(sinh(x) + I)
    F = 4*I*tanh(x)**3/15 - 4*I*tanh(x)/5 - I*sech(x)**3/(5*sinh(x) + 5*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_170():
    f = sech(x)**5/(sinh(x) + I)
    F = -5*I*atan(sinh(x))/16 - 3*I/(16*sinh(x) + 16*I) + 3/(32*(sinh(x) + I)**2) + I/(24*(sinh(x) + I)**3) - 1/(32*(-sinh(x) + I)**2) + I/(-8*sinh(x) + 8*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_171():
    f = cosh(x)**6/(sinh(x) + I)**2
    F = -5*x/8 - 5*sinh(x)*cosh(x)/8 - 5*I*cosh(x)**3/12 + cosh(x)**5/(4*sinh(x) + 4*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_172():
    f = cosh(x)**5/(sinh(x) + I)**2
    F = -(-sinh(x) + I)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_173():
    f = cosh(x)**4/(sinh(x) + I)**2
    F = -3*x/2 - 3*I*cosh(x)/2 + cosh(x)**3/(2*sinh(x) + 2*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_174():
    f = cosh(x)**3/(sinh(x) + I)**2
    F = -2*I*log(sinh(x) + I) + sinh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_175():
    f = cosh(x)**2/(sinh(x) + I)**2
    F = x - 2*cosh(x)/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_176():
    f = cosh(x)/(sinh(x) + I)**2
    F = -1/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_177():
    f = sech(x)/(sinh(x) + I)**2
    F = -atan(sinh(x))/4 - 1/(4*sinh(x) + 4*I) - I/(4*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_178():
    f = sech(x)**2/(sinh(x) + I)**2
    F = -2*tanh(x)/5 - sech(x)/(5*sinh(x) + 5*I) - I*sech(x)/(5*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_179():
    f = sech(x)**3/(sinh(x) + I)**2
    F = -atan(sinh(x))/4 - 3/(16*sinh(x) + 16*I) - I/(8*(sinh(x) + I)**2) + 1/(12*(sinh(x) + I)**3) + 1/(-16*sinh(x) + 16*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_180():
    f = sech(x)**4/(sinh(x) + I)**2
    F = 4*tanh(x)**3/21 - 4*tanh(x)/7 - sech(x)**3/(7*sinh(x) + 7*I) - I*sech(x)**3/(7*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_181():
    f = cosh(x)**3/(I*sinh(x) + 1)**3
    F = I*log(-sinh(x) + I) + 2*I/(I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_182():
    f = cosh(x)**2/(I*sinh(x) + 1)**3
    F = I*cosh(x)**3/(3*(I*sinh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_183():
    f = cosh(x)/(I*sinh(x) + 1)**3
    F = I/(2*(I*sinh(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_184():
    f = cosh(x)**3/(-I*sinh(x) + 1)**3
    F = -I*log(sinh(x) + I) - 2*I/(-I*sinh(x) + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_185():
    f = cosh(x)**2/(-I*sinh(x) + 1)**3
    F = -I*cosh(x)**3/(3*(-I*sinh(x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_186():
    f = cosh(x)/(-I*sinh(x) + 1)**3
    F = -I/(2*(-I*sinh(x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_187():
    f = cosh(x)**7/(a + b*sinh(x))
    F = -a*sinh(x)**5/(5*b**2) - a*(a**2 + 3*b**2)*sinh(x)**3/(3*b**4) - a*(a**4 + 3*a**2*b**2 + 3*b**4)*sinh(x)/b**6 + sinh(x)**6/(6*b) + (a**2 + 3*b**2)*sinh(x)**4/(4*b**3) + (a**4 + 3*a**2*b**2 + 3*b**4)*sinh(x)**2/(2*b**5) + (a**2 + b**2)**3*log(a + b*sinh(x))/b**7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_188():
    f = cosh(x)**6/(a + b*sinh(x))
    F = -a*x*(8*a**4 + 20*a**2*b**2 + 15*b**4)/(8*b**6) + cosh(x)**5/(5*b) + (4*a**2 - 3*a*b*sinh(x) + 4*b**2)*cosh(x)**3/(12*b**3) + (-a*b*(4*a**2 + 7*b**2)*sinh(x) + 8*(a**2 + b**2)**2)*cosh(x)/(8*b**5) - 2*(a**2 + b**2)**(sympy.S(5)/2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/b**6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_189():
    f = cosh(x)**5/(a + b*sinh(x))
    F = -a*sinh(x)**3/(3*b**2) - a*(a**2 + 2*b**2)*sinh(x)/b**4 + sinh(x)**4/(4*b) + (a**2 + 2*b**2)*sinh(x)**2/(2*b**3) + (a**2 + b**2)**2*log(a + b*sinh(x))/b**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_190():
    f = cosh(x)**4/(a + b*sinh(x))
    F = -a*x*(2*a**2 + 3*b**2)/(2*b**4) + cosh(x)**3/(3*b) + (2*a**2 - a*b*sinh(x) + 2*b**2)*cosh(x)/(2*b**3) - 2*(a**2 + b**2)**(sympy.S(3)/2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/b**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_191():
    f = cosh(x)**3/(a + b*sinh(x))
    F = -a*sinh(x)/b**2 + sinh(x)**2/(2*b) + (a**2 + b**2)*log(a + b*sinh(x))/b**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_192():
    f = cosh(x)**2/(a + b*sinh(x))
    F = -a*x/b**2 + cosh(x)/b - 2*sqrt(a**2 + b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_193():
    f = cosh(x)/(a + b*sinh(x))
    F = log(a + b*sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_194():
    f = sech(x)/(a + b*sinh(x))
    F = a*atan(sinh(x))/(a**2 + b**2) + b*log(a + b*sinh(x))/(a**2 + b**2) - b*log(cosh(x))/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_195():
    f = sech(x)**2/(a + b*sinh(x))
    F = -2*b**2*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) + (a*sinh(x) + b)*sech(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_196():
    f = sech(x)**3/(a + b*sinh(x))
    F = a*(a**2 + 3*b**2)*atan(sinh(x))/(2*(a**2 + b**2)**2) + b**3*log(a + b*sinh(x))/(a**2 + b**2)**2 - b**3*log(cosh(x))/(a**2 + b**2)**2 + (a*sinh(x) + b)*sech(x)**2/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_197():
    f = sech(x)**4/(a + b*sinh(x))
    F = -2*b**4*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) + (a*sinh(x) + b)*sech(x)**3/(3*a**2 + 3*b**2) + (a*(2*a**2 + 5*b**2)*sinh(x) + 3*b**3)*sech(x)/(3*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_198():
    f = sech(x)**5/(a + b*sinh(x))
    F = a*(3*a**4 + 10*a**2*b**2 + 15*b**4)*atan(sinh(x))/(8*(a**2 + b**2)**3) + b**5*log(a + b*sinh(x))/(a**2 + b**2)**3 - b**5*log(cosh(x))/(a**2 + b**2)**3 + (a*sinh(x) + b)*sech(x)**4/(4*a**2 + 4*b**2) + (a*(3*a**2 + 7*b**2)*sinh(x) + 4*b**3)*sech(x)**2/(8*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_199():
    f = sech(x)**6/(a + b*sinh(x))
    F = -2*b**6*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(7)/2) + (a*sinh(x) + b)*sech(x)**5/(5*a**2 + 5*b**2) + (a*(4*a**2 + 9*b**2)*sinh(x) + 5*b**3)*sech(x)**3/(15*(a**2 + b**2)**2) + (a*(8*a**4 + 26*a**2*b**2 + 33*b**4)*sinh(x) + 15*b**5)*sech(x)/(15*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_200():
    f = cosh(x)**4/(a + b*sinh(x))**2
    F = 6*a*sqrt(a**2 + b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/b**4 - cosh(x)**3/(b*(a + b*sinh(x))) - 3*(2*a - b*sinh(x))*cosh(x)/(2*b**3) + x*(6*a**2 + 3*b**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_201():
    f = cosh(x)**3/(a + b*sinh(x))**2
    F = -2*a*log(a + b*sinh(x))/b**3 + sinh(x)/b**2 - (a**2 + b**2)/(b**3*(a + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_202():
    f = cosh(x)**2/(a + b*sinh(x))**2
    F = 2*a*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(b**2*sqrt(a**2 + b**2)) - cosh(x)/(b*(a + b*sinh(x))) + x/b**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_203():
    f = cosh(x)/(a + b*sinh(x))**2
    F = -1/(b*(a + b*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_204():
    f = sech(x)/(a + b*sinh(x))**2
    F = 2*a*b*log(a + b*sinh(x))/(a**2 + b**2)**2 - 2*a*b*log(cosh(x))/(a**2 + b**2)**2 - b/((a + b*sinh(x))*(a**2 + b**2)) + (a**2 - b**2)*atan(sinh(x))/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_205():
    f = sech(x)**2/(a + b*sinh(x))**2
    F = -6*a*b**2*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - b*sech(x)/((a + b*sinh(x))*(a**2 + b**2)) + (3*a*b + (a**2 - 2*b**2)*sinh(x))*sech(x)/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_206():
    f = sech(x)**3/(a + b*sinh(x))**2
    F = 4*a*b**3*log(a + b*sinh(x))/(a**2 + b**2)**3 - 4*a*b**3*log(cosh(x))/(a**2 + b**2)**3 + b*(a**2 - 3*b**2)/(2*(a + b*sinh(x))*(a**2 + b**2)**2) + (a**4 + 6*a**2*b**2 - 3*b**4)*atan(sinh(x))/(2*(a**2 + b**2)**3) + (a*sinh(x) + b)*sech(x)**2/((a + b*sinh(x))*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_207():
    f = sech(x)**4/(a + b*sinh(x))**2
    F = -10*a*b**4*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(7)/2) - b*sech(x)**3/((a + b*sinh(x))*(a**2 + b**2)) + (5*a*b + (a**2 - 4*b**2)*sinh(x))*sech(x)**3/(3*(a**2 + b**2)**2) + (15*a*b**3 + (2*a**4 + 9*a**2*b**2 - 8*b**4)*sinh(x))*sech(x)/(3*(a**2 + b**2)**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_208():
    f = tanh(x)**4/(sinh(x) + I)
    F = -I*tanh(x)**5/5 - sech(x)**5/5 + 2*sech(x)**3/3 - sech(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_209():
    f = tanh(x)**3/(sinh(x) + I)
    F = -I*tanh(x)**4/4 - tanh(x)**3*sech(x)/4 - 3*tanh(x)*sech(x)/8 + 3*atan(sinh(x))/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_210():
    f = tanh(x)**2/(sinh(x) + I)
    F = -I*tanh(x)**3/3 + sech(x)**3/3 - sech(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_211():
    f = tanh(x)/(sinh(x) + I)
    F = -tanh(x)*sech(x)/2 + atan(sinh(x))/2 + I*sech(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_212():
    f = coth(x)/(sinh(x) + I)
    F = I*log(sinh(x) + I) - I*log(sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_213():
    f = coth(x)**2/(sinh(x) + I)
    F = I*coth(x) - atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_214():
    f = coth(x)**3/(sinh(x) + I)
    F = I*csch(x)**2/2 - csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_215():
    f = coth(x)**4/(sinh(x) + I)
    F = I*coth(x)**3/3 - coth(x)*csch(x)/2 - atanh(cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_216():
    f = coth(x)**5/(sinh(x) + I)
    F = I*coth(x)**4/4 - csch(x)**3/3 - csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_217():
    f = coth(x)**6/(sinh(x) + I)
    F = I*coth(x)**5/5 - coth(x)**3*csch(x)/4 - 3*coth(x)*csch(x)/8 - 3*atanh(cosh(x))/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_218():
    f = tanh(x)**4/(sinh(x) + I)**2
    F = 2*tanh(x)**7/7 - tanh(x)**5/5 + 2*I*sech(x)**7/7 - 4*I*sech(x)**5/5 + 2*I*sech(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_219():
    f = tanh(x)**3/(sinh(x) + I)**2
    F = -I*atan(sinh(x))/8 - 3*I/(16*sinh(x) + 16*I) - 1/(4*(sinh(x) + I)**2) + I/(12*(sinh(x) + I)**3) - I/(-16*sinh(x) + 16*I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_220():
    f = tanh(x)**2/(sinh(x) + I)**2
    F = 2*tanh(x)**5/5 - tanh(x)**3/3 - 2*I*sech(x)**5/5 + 2*I*sech(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_221():
    f = tanh(x)/(sinh(x) + I)**2
    F = -I*atan(sinh(x))/4 - I/(4*sinh(x) + 4*I) - 1/(4*(sinh(x) + I)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_222():
    f = coth(x)/(sinh(x) + I)**2
    F = log(sinh(x) + I) - log(sinh(x)) - I/(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_223():
    f = coth(x)**3/(sinh(x) + I)**2
    F = -2*log(sinh(x) + I) + 2*log(sinh(x)) + csch(x)**2/2 + 2*I*csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_224():
    f = coth(x)**4/(sinh(x) + I)**2
    F = coth(x)**3/3 + I*coth(x)*csch(x) - 2*coth(x) - I*atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_225():
    f = coth(x)**5/(sinh(x) + I)**2
    F = csch(x)**4/4 + 2*I*csch(x)**3/3 - csch(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_226():
    f = coth(x)**6/(sinh(x) + I)**2
    F = coth(x)**5/5 - 2*coth(x)**3/3 + I*coth(x)*csch(x)**3/2 + I*coth(x)*csch(x)/4 - I*atanh(cosh(x))/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_227():
    f = tanh(x)**4/(a + b*sinh(x))
    F = -2*a**4*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - a**3*tanh(x)/(a**2 + b**2)**2 - a**2*b*sech(x)/(a**2 + b**2)**2 - a*tanh(x)**3/(3*a**2 + 3*b**2) + b*sech(x)**3/(3*a**2 + 3*b**2) - b*sech(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_228():
    f = tanh(x)**3/(a + b*sinh(x))
    F = -a**3*log(a + b*sinh(x))/(a**2 + b**2)**2 + a**3*log(cosh(x))/(a**2 + b**2)**2 + b*(3*a**2 + b**2)*atan(sinh(x))/(2*(a**2 + b**2)**2) + (a - b*sinh(x))*sech(x)**2/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_229():
    f = tanh(x)**2/(a + b*sinh(x))
    F = -2*a**2*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) - a*tanh(x)/(a**2 + b**2) - b*sech(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_230():
    f = tanh(x)/(a + b*sinh(x))
    F = -a*log(a + b*sinh(x))/(a**2 + b**2) + a*log(cosh(x))/(a**2 + b**2) + b*atan(sinh(x))/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_231():
    f = coth(x)/(a + b*sinh(x))
    F = -log(a + b*sinh(x))/a + log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_232():
    f = coth(x)**2/(a + b*sinh(x))
    F = -coth(x)/a + b*atanh(cosh(x))/a**2 - 2*sqrt(a**2 + b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_233():
    f = coth(x)**3/(a + b*sinh(x))
    F = -csch(x)**2/(2*a) + b*csch(x)/a**2 - (a**2 + b**2)*log(a + b*sinh(x))/a**3 + (a**2 + b**2)*log(sinh(x))/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_234():
    f = coth(x)**4/(a + b*sinh(x))
    F = -coth(x)*csch(x)**2/(3*a) + b*coth(x)*csch(x)/(2*a**2) - (4*a**2 + 3*b**2)*coth(x)/(3*a**3) + b*(3*a**2 + 2*b**2)*atanh(cosh(x))/(2*a**4) - 2*(a**2 + b**2)**(sympy.S(3)/2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_235():
    f = tanh(x)**4/(a + b*sinh(x))**2
    F = -2*a**5*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(7)/2) - a**4*b*cosh(x)/((a + b*sinh(x))*(a**2 + b**2)**3) + 8*a**3*b**2*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(7)/2) - 4*a**3*b*sech(x)/(a**2 + b**2)**3 + 2*a*b*sech(x)**3/(3*(a**2 + b**2)**2) - (a**2 - b**2)*tanh(x)**3/(3*(a**2 + b**2)**2) + (a**2 - b**2)*tanh(x)/(a**2 + b**2)**2 - (2*a**4 - 3*a**2*b**2 - b**4)*tanh(x)/(a**2 + b**2)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_236():
    f = tanh(x)**3/(a + b*sinh(x))**2
    F = a**3/((a + b*sinh(x))*(a**2 + b**2)**2) - a**2*(a**2 - 3*b**2)*log(a + b*sinh(x))/(a**2 + b**2)**3 + a**2*(a**2 - 3*b**2)*log(cosh(x))/(a**2 + b**2)**3 + a*b*(3*a**2 - b**2)*atan(sinh(x))/(a**2 + b**2)**3 + (a**2 - 2*a*b*sinh(x) - b**2)*sech(x)**2/(2*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_237():
    f = tanh(x)**2/(a + b*sinh(x))**2
    F = -2*a**3*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - a**2*b*cosh(x)/((a + b*sinh(x))*(a**2 + b**2)**2) + 4*a*b**2*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - 2*a*b*sech(x)/(a**2 + b**2)**2 - (a**2 - b**2)*tanh(x)/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_238():
    f = tanh(x)/(a + b*sinh(x))**2
    F = 2*a*b*atan(sinh(x))/(a**2 + b**2)**2 + a/((a + b*sinh(x))*(a**2 + b**2)) - (a**2 - b**2)*log(a + b*sinh(x))/(a**2 + b**2)**2 + (a**2 - b**2)*log(cosh(x))/(a**2 + b**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_239():
    f = coth(x)/(a + b*sinh(x))**2
    F = 1/(a*(a + b*sinh(x))) - log(a + b*sinh(x))/a**2 + log(sinh(x))/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_240():
    f = coth(x)**2/(a + b*sinh(x))**2
    F = coth(x)/(a*(a + b*sinh(x))) - 2*coth(x)/a**2 + 2*b*atanh(cosh(x))/a**3 - (2*a**2 + 4*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a**3*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_241():
    f = coth(x)**3/(a + b*sinh(x))**2
    F = -csch(x)**2/(2*a**2) + 2*b*csch(x)/a**3 + (a**2 + b**2)/(a**3*(a + b*sinh(x))) - (a**2 + 3*b**2)*log(a + b*sinh(x))/a**4 + (a**2 + 3*b**2)*log(sinh(x))/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_242():
    f = coth(x)**4/(a + b*sinh(x))**2
    F = -(3 + 4*b**2/a**2)*coth(x)*csch(x)/(3*b*(a + b*sinh(x))) - coth(x)*csch(x)**2/(3*a*(a + b*sinh(x))) + (a**2 + 2*b**2)*coth(x)*csch(x)/(a**3*b) - (7*a**2 + 12*b**2)*coth(x)/(3*a**4) + b*(3*a**2 + 4*b**2)*atanh(cosh(x))/a**5 - 2*sqrt(a**2 + b**2)*(a**2 + 4*b**2)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/a**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_243():
    f = sqrt(a + b*sinh(x))*coth(x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*sinh(x))/sqrt(a)) + 2*sqrt(a + b*sinh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_244():
    f = coth(x)/sqrt(a + b*sinh(x))
    F = -2*atanh(sqrt(a + b*sinh(x))/sqrt(a))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_245():
    f = (A + B*cosh(x))/(a + b*sinh(x))
    F = -2*A*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2) + B*log(a + b*sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_246():
    f = (A + B*cosh(x))/(sinh(x) + I)
    F = -A*cosh(x)/(-I*sinh(x) + 1) + B*log(sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_247():
    f = (A + B*cosh(x))/(-sinh(x) + I)
    F = A*cosh(x)/(I*sinh(x) + 1) - B*log(-sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_248():
    f = (A + B*tanh(x))/(a + b*sinh(x))
    F = -2*A*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2) - B*a*log(a + b*sinh(x))/(a**2 + b**2) + B*a*log(cosh(x))/(a**2 + b**2) + B*b*atan(sinh(x))/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_249():
    f = (A + B*coth(x))/(a + b*sinh(x))
    F = -2*A*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2) - B*log(a + b*sinh(x))/a + B*log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_250():
    f = (A + B*sech(x))/(a + b*sinh(x))
    F = -2*A*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/sqrt(a**2 + b**2) + B*a*atan(sinh(x))/(a**2 + b**2) + B*b*log(a + b*sinh(x))/(a**2 + b**2) - B*b*log(cosh(x))/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_251():
    f = (A + B*csch(x))/(a + b*sinh(x))
    F = -B*atanh(cosh(x))/a - (2*A*a - 2*B*b)*atanh((-a*tanh(x/2) + b)/sqrt(a**2 + b**2))/(a*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_252():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + c*sinh(d + e*x))
    F = B*log(a + c*sinh(d + e*x))/(c*e) + C*x/c - (2*A*c - 2*C*a)*atanh((-a*tanh(d/2 + e*x/2) + c)/sqrt(a**2 + c**2))/(c*e*sqrt(a**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_253():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + c*sinh(d + e*x))**2
    F = -B/(c*e*(a + c*sinh(d + e*x))) - (2*A*a + 2*C*c)*atanh((-a*tanh(d/2 + e*x/2) + c)/sqrt(a**2 + c**2))/(e*(a**2 + c**2)**(sympy.S(3)/2)) - (A*c - C*a)*cosh(d + e*x)/(e*(a + c*sinh(d + e*x))*(a**2 + c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_254():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + c*sinh(d + e*x))**3
    F = -B/(2*c*e*(a + c*sinh(d + e*x))**2) - (2*A*a**2 - A*c**2 + 3*C*a*c)*atanh((-a*tanh(d/2 + e*x/2) + c)/sqrt(a**2 + c**2))/(e*(a**2 + c**2)**(sympy.S(5)/2)) - (3*A*a*c - C*a**2 + 2*C*c**2)*cosh(d + e*x)/(2*e*(a + c*sinh(d + e*x))*(a**2 + c**2)**2) - (A*c - C*a)*cosh(d + e*x)/(e*(a + c*sinh(d + e*x))**2*(2*a**2 + 2*c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_255():
    f = (A + B*cosh(d + e*x) + C*sinh(d + e*x))/(a + c*sinh(d + e*x))**4
    F = -B/(3*c*e*(a + c*sinh(d + e*x))**3) - (2*A*a**3 - 3*A*a*c**2 + 4*C*a**2*c - C*c**3)*atanh((-a*tanh(d/2 + e*x/2) + c)/sqrt(a**2 + c**2))/(e*(a**2 + c**2)**(sympy.S(7)/2)) - (11*A*a**2*c - 4*A*c**3 - 2*C*a**3 + 13*C*a*c**2)*cosh(d + e*x)/(6*e*(a + c*sinh(d + e*x))*(a**2 + c**2)**3) - (5*A*a*c - 2*C*a**2 + 3*C*c**2)*cosh(d + e*x)/(6*e*(a + c*sinh(d + e*x))**2*(a**2 + c**2)**2) - (A*c - C*a)*cosh(d + e*x)/(e*(a + c*sinh(d + e*x))**3*(3*a**2 + 3*c**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_256():
    f = x**3/(a + b*sinh(x)**2)
    F = (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(3)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))) + ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (x)**(Integer(2)) * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))) + ((Integer(3) * x * sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.Function('PolyLog')(Integer(4), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_257():
    f = x**2/(a + b*sinh(x)**2)
    F = (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + (Integer(-1) * (((x)**(Integer(2)) * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))) + ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(3), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_258():
    f = x/(a + b*sinh(x)**2)
    F = ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + (Integer(-1) * ((x * sympy.log((Integer(1) + ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1)))) + (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))) + (Integer(-1) * (sympy.Function('PolyLog')(Integer(2), (Integer(-1) * ((Symbol('b') * (sympy.E)**((Integer(2) * x))) * (((Integer(2) * Symbol('a')) + (Integer(2) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))) + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b'))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_259():
    f = (sinh(a + b*x)**2 - 2)*cosh(a + b*x)/x
    F = ((Integer(-1) * (Integer(9) * (Integer(4))**(Integer(-1)))) * sympy.cosh(Symbol('a')) * sympy.Function('CoshIntegral')((Symbol('b') * x))) + ((Integer(4))**(Integer(-1)) * sympy.cosh((Integer(3) * Symbol('a'))) * sympy.Function('CoshIntegral')((Integer(3) * Symbol('b') * x))) + (Integer(-1) * ((Integer(9) * (Integer(4))**(Integer(-1))) * sympy.sinh(Symbol('a')) * sympy.Function('SinhIntegral')((Symbol('b') * x)))) + ((Integer(4))**(Integer(-1)) * sympy.sinh((Integer(3) * Symbol('a'))) * sympy.Function('SinhIntegral')((Integer(3) * Symbol('b') * x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_260():
    f = sinh(sqrt(-a*x + 1)/sqrt(a*x + 1))**3/(-a**2*x**2 + 1)
    F = ((Integer(3) * sympy.Function('SinhIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1))))) * ((Integer(4) * Symbol('a')))**(Integer(-1))) + (Integer(-1) * (sympy.Function('SinhIntegral')(((Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(4) * Symbol('a')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_261():
    f = sinh(sqrt(-a*x + 1)/sqrt(a*x + 1))**2/(-a**2*x**2 + 1)
    F = (Integer(-1) * (sympy.Function('CoshIntegral')(((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x))))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))) + (sympy.log((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * ((Integer(2) * Symbol('a')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_262():
    f = sinh(sqrt(-a*x + 1)/sqrt(a*x + 1))/(-a**2*x**2 + 1)
    F = Integer(-1) * (sympy.Function('SinhIntegral')((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (Symbol('a'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_263():
    f = 1/((-a**2*x**2 + 1)*sinh(sqrt(-a*x + 1)/sqrt(a*x + 1)))
    F = sympy.Function('Unintegrable')((sympy.csch((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_264():
    f = 1/((-a**2*x**2 + 1)*sinh(sqrt(-a*x + 1)/sqrt(a*x + 1))**2)
    F = sympy.Function('Unintegrable')(((sympy.csch((sympy.sqrt((Integer(1) + (Integer(-1) * (Symbol('a') * x)))) * (sympy.sqrt((Integer(1) + (Symbol('a') * x))))**(Integer(-1)))))**(Integer(2)) * (((Integer(1) + (Integer(-1) * (Symbol('a') * x))) * (Integer(1) + (Symbol('a') * x))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_265():
    f = sinh(a + b*log(c*x**n))
    F = -b*n*x*cosh(a + b*log(c*x**n))/(-b**2*n**2 + 1) + x*sinh(a + b*log(c*x**n))/(-b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_266():
    f = sinh(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x/(-4*b**2*n**2 + 1) - 2*b*n*x*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(-4*b**2*n**2 + 1) + x*sinh(a + b*log(c*x**n))**2/(-4*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_267():
    f = sinh(a + b*log(c*x**n))**3
    F = -6*b**3*n**3*x*cosh(a + b*log(c*x**n))/(9*b**4*n**4 - 10*b**2*n**2 + 1) + 6*b**2*n**2*x*sinh(a + b*log(c*x**n))/(9*b**4*n**4 - 10*b**2*n**2 + 1) - 3*b*n*x*sinh(a + b*log(c*x**n))**2*cosh(a + b*log(c*x**n))/(-9*b**2*n**2 + 1) + x*sinh(a + b*log(c*x**n))**3/(-9*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_268():
    f = sinh(a + b*log(c*x**n))**4
    F = 24*b**4*n**4*x/(64*b**4*n**4 - 20*b**2*n**2 + 1) - 24*b**3*n**3*x*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(64*b**4*n**4 - 20*b**2*n**2 + 1) + 12*b**2*n**2*x*sinh(a + b*log(c*x**n))**2/(64*b**4*n**4 - 20*b**2*n**2 + 1) - 4*b*n*x*sinh(a + b*log(c*x**n))**3*cosh(a + b*log(c*x**n))/(-16*b**2*n**2 + 1) + x*sinh(a + b*log(c*x**n))**4/(-16*b**2*n**2 + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_269():
    f = x**m*sinh(a + b*log(c*x**n))
    F = -b*n*x**(m + 1)*cosh(a + b*log(c*x**n))/(-b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*sinh(a + b*log(c*x**n))/(-b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_270():
    f = x**m*sinh(a + b*log(c*x**n))**2
    F = 2*b**2*n**2*x**(m + 1)/((m + 1)*(-4*b**2*n**2 + (m + 1)**2)) - 2*b*n*x**(m + 1)*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(-4*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*sinh(a + b*log(c*x**n))**2/(-4*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_271():
    f = x**m*sinh(a + b*log(c*x**n))**3
    F = -6*b**3*n**3*x**(m + 1)*cosh(a + b*log(c*x**n))/((-9*b**2*n**2 + (m + 1)**2)*(-b**2*n**2 + (m + 1)**2)) + 6*b**2*n**2*x**(m + 1)*(m + 1)*sinh(a + b*log(c*x**n))/((-9*b**2*n**2 + (m + 1)**2)*(-b**2*n**2 + (m + 1)**2)) - 3*b*n*x**(m + 1)*sinh(a + b*log(c*x**n))**2*cosh(a + b*log(c*x**n))/(-9*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*sinh(a + b*log(c*x**n))**3/(-9*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_272():
    f = x**m*sinh(a + b*log(c*x**n))**4
    F = 24*b**4*n**4*x**(m + 1)/((m + 1)*(-16*b**2*n**2 + (m + 1)**2)*(-4*b**2*n**2 + (m + 1)**2)) - 24*b**3*n**3*x**(m + 1)*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/((-16*b**2*n**2 + (m + 1)**2)*(-4*b**2*n**2 + (m + 1)**2)) + 12*b**2*n**2*x**(m + 1)*(m + 1)*sinh(a + b*log(c*x**n))**2/((-16*b**2*n**2 + (m + 1)**2)*(-4*b**2*n**2 + (m + 1)**2)) - 4*b*n*x**(m + 1)*sinh(a + b*log(c*x**n))**3*cosh(a + b*log(c*x**n))/(-16*b**2*n**2 + (m + 1)**2) + x**(m + 1)*(m + 1)*sinh(a + b*log(c*x**n))**4/(-16*b**2*n**2 + (m + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_273():
    f = sinh(a + b*log(c*x**n))/x
    F = cosh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_274():
    f = sinh(a + b*log(c*x**n))**2/x
    F = -log(x)/2 + sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_275():
    f = sinh(a + b*log(c*x**n))**3/x
    F = cosh(a + b*log(c*x**n))**3/(3*b*n) - cosh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_276():
    f = sinh(a + b*log(c*x**n))**4/x
    F = 3*log(x)/8 + sinh(a + b*log(c*x**n))**3*cosh(a + b*log(c*x**n))/(4*b*n) - 3*sinh(a + b*log(c*x**n))*cosh(a + b*log(c*x**n))/(8*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_277():
    f = sinh(a + b*log(c*x**n))**5/x
    F = cosh(a + b*log(c*x**n))**5/(5*b*n) - 2*cosh(a + b*log(c*x**n))**3/(3*b*n) + cosh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_278():
    f = sinh(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = 2*sinh(a + b*log(c*x**n))**(sympy.S(3)/2)*cosh(a + b*log(c*x**n))/(5*b*n) + 6*I*sqrt(sinh(a + b*log(c*x**n)))*elliptic_e(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(5*b*n*sqrt(I*sinh(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_279():
    f = sinh(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = 2*I*sqrt(I*sinh(a + b*log(c*x**n)))*elliptic_f(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(3*b*n*sqrt(sinh(a + b*log(c*x**n)))) + 2*sqrt(sinh(a + b*log(c*x**n)))*cosh(a + b*log(c*x**n))/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_280():
    f = sqrt(sinh(a + b*log(c*x**n)))/x
    F = -2*I*sqrt(sinh(a + b*log(c*x**n)))*elliptic_e(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(b*n*sqrt(I*sinh(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_281():
    f = 1/(x*sqrt(sinh(a + b*log(c*x**n))))
    F = -2*I*sqrt(I*sinh(a + b*log(c*x**n)))*elliptic_f(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(b*n*sqrt(sinh(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_282():
    f = 1/(x*sinh(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = -2*cosh(a + b*log(c*x**n))/(b*n*sqrt(sinh(a + b*log(c*x**n)))) - 2*I*sqrt(sinh(a + b*log(c*x**n)))*elliptic_e(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(b*n*sqrt(I*sinh(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_283():
    f = 1/(x*sinh(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 2*I*sqrt(I*sinh(a + b*log(c*x**n)))*elliptic_f(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(3*b*n*sqrt(sinh(a + b*log(c*x**n)))) - 2*cosh(a + b*log(c*x**n))/(3*b*n*sinh(a + b*log(c*x**n))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_284():
    f = sinh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)
    F = -x*sinh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)/4 + 5*x*sinh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)/(12 - 12*exp(-2*a)/(c*x**n)**(4/n)) - 5*x*exp(-2*a)*sinh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)/(4*(c*x**n)**(4/n)*(1 - exp(-2*a)/(c*x**n)**(4/n))**2) - 5*x*exp(-3*a)*sinh(a + 2*log(c*x**n)/n)**(sympy.S(5)/2)*acsc((c*x**n)**(2/n)*exp(a))/(4*(c*x**n)**(6/n)*(1 - exp(-2*a)/(c*x**n)**(4/n))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_285():
    f = sqrt(sinh(a + 2*log(c*x**n)/n))
    F = x*sqrt(sinh(a + 2*log(c*x**n)/n))/2 + x*exp(-a)*sqrt(sinh(a + 2*log(c*x**n)/n))*acsc((c*x**n)**(2/n)*exp(a))/(2*(c*x**n)**(2/n)*sqrt(1 - exp(-2*a)/(c*x**n)**(4/n)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_286():
    f = sinh(a + 2*log(c*x**n)/n)**(sympy.S(-3)/2)
    F = -x*(1 - exp(-2*a)/(c*x**n)**(4/n))/(2*sinh(a + 2*log(c*x**n)/n)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_287():
    f = sinh(a + 2*log(c*x**n)/n)**(sympy.S(-7)/2)
    F = -x*(1 - exp(-2*a)/(c*x**n)**(4/n))/(6*sinh(a + 2*log(c*x**n)/n)**(sympy.S(7)/2)) + x*(1 - exp(-2*a)/(c*x**n)**(4/n))*exp(-2*a)/(15*(c*x**n)**(4/n)*sinh(a + 2*log(c*x**n)/n)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_288():
    f = sinh(a/(c + d*x))
    F = (Integer(-1) * ((Symbol('a') * sympy.Function('CoshIntegral')((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sinh((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_289():
    f = sinh(a/(c + d*x))**2
    F = (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_290():
    f = sinh(a/(c + d*x))**3
    F = ((Integer(3) * Symbol('a') * sympy.Function('CoshIntegral')((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.Function('CoshIntegral')(((Integer(3) * Symbol('a')) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh((Symbol('a') * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(3))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_291():
    f = sinh(b*x/(c + d*x))
    F = ((Symbol('b') * Symbol('c') * sympy.cosh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CoshIntegral')(((Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sinh(((Symbol('b') * x) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * sympy.sinh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_292():
    f = sinh(b*x/(c + d*x))**2
    F = ((Symbol('b') * Symbol('c') * sympy.Function('CoshIntegral')(((Integer(2) * Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sinh(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh(((Symbol('b') * x) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * Symbol('c') * sympy.cosh(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Integer(2) * Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_293():
    f = sinh(b*x/(c + d*x))**3
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c') * sympy.cosh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CoshIntegral')(((Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('b') * Symbol('c') * sympy.cosh(((Integer(3) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CoshIntegral')(((Integer(3) * Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh(((Symbol('b') * x) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + ((Integer(3) * Symbol('b') * Symbol('c') * sympy.sinh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * Symbol('c') * sympy.sinh(((Integer(3) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Integer(3) * Symbol('b') * Symbol('c')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_294():
    f = sinh((a + b*x)/(c + d*x))
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cosh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sinh(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sinh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_295():
    f = sinh((a + b*x)/(c + d*x))**2
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('CoshIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sinh(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cosh(((Integer(2) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_296():
    f = sinh((a + b*x)/(c + d*x))**3
    F = (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cosh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CoshIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cosh(((Integer(3) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('CoshIntegral')(((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh(((Symbol('a') + (Symbol('b') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sinh((Symbol('b') * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sinh(((Integer(3) * Symbol('b')) * (Symbol('d'))**(Integer(-1)))) * sympy.Function('SinhIntegral')(((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d'))))) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_297():
    f = sinh(e + f*(a + b*x)/(c + d*x))
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.cosh((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('CoshIntegral')(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * sympy.sinh((((Symbol('c') * Symbol('e')) + (Symbol('a') * Symbol('f')) + (Symbol('d') * Symbol('e') * x) + (Symbol('b') * Symbol('f') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sinh((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('SinhIntegral')(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_298():
    f = sinh(e + f*(a + b*x)/(c + d*x))**2
    F = ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.Function('CoshIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1)))) * sympy.sinh((Integer(2) * (Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh((((Symbol('c') * Symbol('e')) + (Symbol('a') * Symbol('f')) + (Symbol('d') * Symbol('e') * x) + (Symbol('b') * Symbol('f') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(2))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.cosh((Integer(2) * (Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')(((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Symbol('d'))**(Integer(2)))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_299():
    f = sinh(e + f*(a + b*x)/(c + d*x))**3
    F = (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.cosh((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('CoshIntegral')(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.cosh((Integer(3) * (Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('CoshIntegral')(((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (((Symbol('c') + (Symbol('d') * x)) * (sympy.sinh((((Symbol('c') * Symbol('e')) + (Symbol('a') * Symbol('f')) + (Symbol('d') * Symbol('e') * x) + (Symbol('b') * Symbol('f') * x)) * ((Symbol('c') + (Symbol('d') * x)))**(Integer(-1)))))**(Integer(3))) * (Symbol('d'))**(Integer(-1))) + ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sinh((Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1))))) * sympy.Function('SinhIntegral')(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sinh((Integer(3) * (Symbol('e') + ((Symbol('b') * Symbol('f')) * (Symbol('d'))**(Integer(-1)))))) * sympy.Function('SinhIntegral')(((Integer(3) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')) * ((Symbol('d') * (Symbol('c') + (Symbol('d') * x))))**(Integer(-1))))) * ((Integer(4) * (Symbol('d'))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_300():
    f = exp(a + b*x)*sinh(a + b*x)**4
    F = -exp(-3*a - 3*b*x)/(48*b) + exp(-a - b*x)/(4*b) + 3*exp(a + b*x)/(8*b) - exp(3*a + 3*b*x)/(12*b) + exp(5*a + 5*b*x)/(80*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_301():
    f = exp(a + b*x)*sinh(a + b*x)**3
    F = 3*x/8 + exp(-2*a - 2*b*x)/(16*b) - 3*exp(2*a + 2*b*x)/(16*b) + exp(4*a + 4*b*x)/(32*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_302():
    f = exp(a + b*x)*sinh(a + b*x)**2
    F = -exp(-a - b*x)/(4*b) - exp(a + b*x)/(2*b) + exp(3*a + 3*b*x)/(12*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_303():
    f = exp(a + b*x)*sinh(a + b*x)
    F = -x/2 + exp(2*a + 2*b*x)/(4*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_304():
    f = exp(a + b*x)*csch(a + b*x)
    F = log(1 - exp(2*a + 2*b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_305():
    f = exp(a + b*x)*csch(a + b*x)**2
    F = -2*atanh(exp(a + b*x))/b + 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_306():
    f = exp(a + b*x)*csch(a + b*x)**3
    F = -2*exp(4*a + 4*b*x)/(b*(1 - exp(2*a + 2*b*x))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_307():
    f = exp(a + b*x)*csch(a + b*x)**4
    F = atanh(exp(a + b*x))/b + exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))) - 2*exp(a + b*x)/(b*(1 - exp(2*a + 2*b*x))**2) + 8*exp(3*a + 3*b*x)/(3*b*(1 - exp(2*a + 2*b*x))**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_308():
    f = exp(a + b*x)*csch(a + b*x)**5
    F = -8/(b*(1 - exp(2*a + 2*b*x))**2) + 32/(3*b*(1 - exp(2*a + 2*b*x))**3) - 4/(b*(1 - exp(2*a + 2*b*x))**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_309():
    f = exp(x)*sinh(2*x)**2
    F = exp(5*x)/20 - exp(x)/2 - exp(-3*x)/12
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_310():
    f = exp(x)*sinh(2*x)
    F = exp(3*x)/6 + exp(-x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_311():
    f = exp(x)*csch(2*x)
    F = atan(exp(x)) - atanh(exp(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_312():
    f = exp(x)*csch(2*x)**2
    F = -atan(exp(x))/2 - atanh(exp(x))/2 + exp(x)/(1 - exp(4*x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_313():
    f = exp(x)*sinh(3*x)**2
    F = exp(7*x)/28 - exp(x)/2 - exp(-5*x)/20
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_314():
    f = exp(x)*sinh(3*x)
    F = exp(4*x)/8 + exp(-2*x)/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_315():
    f = exp(x)*csch(3*x)
    F = log(1 - exp(2*x))/3 - log(exp(4*x) + exp(2*x) + 1)/6 + sqrt(3)*atan(sqrt(3)*(2*exp(2*x) + 1)/3)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_316():
    f = exp(x)*csch(3*x)**2
    F = log(exp(2*x) - exp(x) + 1)/18 - log(exp(2*x) + exp(x) + 1)/18 + sqrt(3)*atan(sqrt(3)*(1 - 2*exp(x))/3)/9 - sqrt(3)*atan(sqrt(3)*(2*exp(x) + 1)/3)/9 - 2*atanh(exp(x))/9 + 2*exp(x)/(3 - 3*exp(6*x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_317():
    f = exp(x)*sinh(4*x)**2
    F = exp(9*x)/36 - exp(x)/2 - exp(-7*x)/28
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_318():
    f = exp(x)*sinh(4*x)
    F = exp(5*x)/10 + exp(-3*x)/6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_319():
    f = exp(x)*csch(4*x)
    F = -sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/8 + sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/8 + sqrt(2)*atan(sqrt(2)*exp(x) - 1)/4 + sqrt(2)*atan(sqrt(2)*exp(x) + 1)/4 - atan(exp(x))/2 - atanh(exp(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_320():
    f = exp(x)*csch(4*x)**2
    F = sqrt(2)*log(exp(2*x) - sqrt(2)*exp(x) + 1)/32 - sqrt(2)*log(exp(2*x) + sqrt(2)*exp(x) + 1)/32 - sqrt(2)*atan(sqrt(2)*exp(x) - 1)/16 - sqrt(2)*atan(sqrt(2)*exp(x) + 1)/16 - atan(exp(x))/8 - atanh(exp(x))/8 + exp(x)/(2 - 2*exp(8*x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_321():
    f = F**(c*(a + b*x))*sinh(d + e*x)**3
    F = 6*F**(c*(a + b*x))*b*c*e**2*log(F)*sinh(d + e*x)/(b**4*c**4*log(F)**4 - 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) - F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**3/(-b**2*c**2*log(F)**2 + 9*e**2) - 6*F**(c*(a + b*x))*e**3*cosh(d + e*x)/(b**4*c**4*log(F)**4 - 10*b**2*c**2*e**2*log(F)**2 + 9*e**4) + 3*F**(c*(a + b*x))*e*sinh(d + e*x)**2*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + 9*e**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_322():
    f = F**(c*(a + b*x))*sinh(d + e*x)**2
    F = -F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)**2/(-b**2*c**2*log(F)**2 + 4*e**2) + 2*F**(c*(a + b*x))*e*sinh(d + e*x)*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + 4*e**2) - 2*F**(c*(a + b*x))*e**2/(b*c*(-b**2*c**2*log(F)**2 + 4*e**2)*log(F))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_323():
    f = F**(c*(a + b*x))*sinh(d + e*x)
    F = -F**(c*(a + b*x))*b*c*log(F)*sinh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2) + F**(c*(a + b*x))*e*cosh(d + e*x)/(-b**2*c**2*log(F)**2 + e**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_324():
    f = F**(c*(a + b*x))*csch(d + e*x)
    F = -2*F**(c*(a + b*x))*exp(d + e*x)*hyper((1, (b*c*log(F) + e)/(2*e)), (b*c*log(F)/(2*e) + sympy.S(3)/2,), exp(2*d + 2*e*x))/(b*c*log(F) + e)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_325():
    f = F**(c*(a + b*x))*csch(d + e*x)**2
    F = 4*F**(c*(a + b*x))*exp(2*d + 2*e*x)*hyper((2, b*c*log(F)/(2*e) + 1), (b*c*log(F)/(2*e) + 2,), exp(2*d + 2*e*x))/(b*c*log(F) + 2*e)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_326():
    f = F**(c*(a + b*x))*csch(d + e*x)**3
    F = -F**(c*(a + b*x))*b*c*log(F)*csch(d + e*x)/(2*e**2) - F**(c*(a + b*x))*coth(d + e*x)*csch(d + e*x)/(2*e) + F**(c*(a + b*x))*(-b*c*log(F) + e)*exp(d + e*x)*hyper((1, (b*c*log(F) + e)/(2*e)), (b*c*log(F)/(2*e) + sympy.S(3)/2,), exp(2*d + 2*e*x))/e**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_327():
    f = F**(c*(a + b*x))*csch(d + e*x)**4
    F = -F**(c*(a + b*x))*b*c*log(F)*csch(d + e*x)**2/(6*e**2) - F**(c*(a + b*x))*coth(d + e*x)*csch(d + e*x)**2/(3*e) - 2*F**(c*(a + b*x))*(-b*c*log(F) + 2*e)*exp(2*d + 2*e*x)*hyper((2, b*c*log(F)/(2*e) + 1), (b*c*log(F)/(2*e) + 2,), exp(2*d + 2*e*x))/(3*e**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_328():
    f = (sinh(a*c + b*c*x)**2)**(sympy.S(5)/2)*exp(c*(a + b*x))
    F = -5*x*sqrt(sinh(a*c + b*c*x)**2)*csch(a*c + b*c*x)/16 + sqrt(sinh(a*c + b*c*x)**2)*exp(6*c*(a + b*x))*csch(a*c + b*c*x)/(192*b*c) - 5*sqrt(sinh(a*c + b*c*x)**2)*exp(4*c*(a + b*x))*csch(a*c + b*c*x)/(128*b*c) + 5*sqrt(sinh(a*c + b*c*x)**2)*exp(2*c*(a + b*x))*csch(a*c + b*c*x)/(32*b*c) - 5*sqrt(sinh(a*c + b*c*x)**2)*exp(-2*c*(a + b*x))*csch(a*c + b*c*x)/(64*b*c) + sqrt(sinh(a*c + b*c*x)**2)*exp(-4*c*(a + b*x))*csch(a*c + b*c*x)/(128*b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_329():
    f = (sinh(a*c + b*c*x)**2)**(sympy.S(3)/2)*exp(c*(a + b*x))
    F = 3*x*sqrt(sinh(a*c + b*c*x)**2)*csch(a*c + b*c*x)/8 + sqrt(sinh(a*c + b*c*x)**2)*exp(4*c*(a + b*x))*csch(a*c + b*c*x)/(32*b*c) - 3*sqrt(sinh(a*c + b*c*x)**2)*exp(2*c*(a + b*x))*csch(a*c + b*c*x)/(16*b*c) + sqrt(sinh(a*c + b*c*x)**2)*exp(-2*c*(a + b*x))*csch(a*c + b*c*x)/(16*b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_330():
    f = sqrt(sinh(a*c + b*c*x)**2)*exp(c*(a + b*x))
    F = -x*sqrt(sinh(a*c + b*c*x)**2)*csch(a*c + b*c*x)/2 + sqrt(sinh(a*c + b*c*x)**2)*exp(2*c*(a + b*x))*csch(a*c + b*c*x)/(4*b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_331():
    f = exp(c*(a + b*x))/sqrt(sinh(a*c + b*c*x)**2)
    F = log(1 - exp(2*c*(a + b*x)))*sinh(a*c + b*c*x)/(b*c*sqrt(sinh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_332():
    f = exp(c*(a + b*x))/(sinh(a*c + b*c*x)**2)**(sympy.S(3)/2)
    F = -2*exp(4*c*(a + b*x))*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**2*sqrt(sinh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_333():
    f = exp(c*(a + b*x))/(sinh(a*c + b*c*x)**2)**(sympy.S(5)/2)
    F = -8*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**2*sqrt(sinh(a*c + b*c*x)**2)) + 32*sinh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**3*sqrt(sinh(a*c + b*c*x)**2)) - 4*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**4*sqrt(sinh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_334():
    f = exp(c*(a + b*x))/(sinh(a*c + b*c*x)**2)**(sympy.S(7)/2)
    F = 64*sinh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**3*sqrt(sinh(a*c + b*c*x)**2)) - 48*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**4*sqrt(sinh(a*c + b*c*x)**2)) + 192*sinh(a*c + b*c*x)/(5*b*c*(1 - exp(2*c*(a + b*x)))**5*sqrt(sinh(a*c + b*c*x)**2)) - 32*sinh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**6*sqrt(sinh(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_335():
    f = exp(x)*sinh(a + b*x)
    F = -b*exp(x)*cosh(a + b*x)/(1 - b**2) + exp(x)*sinh(a + b*x)/(1 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_336():
    f = exp(x)*sinh(a + c*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('a')) + ((Integer(4) * Symbol('c')))**(Integer(-1)))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(1) + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(1) + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_337():
    f = exp(x)*sinh(a + b*x + c*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('a')) + (((Integer(1) + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Integer(1) + (Integer(-1) * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('c') * x))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * (((Integer(1) + Symbol('b')))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(1) + Symbol('b') + (Integer(2) * Symbol('c') * x)) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_338():
    f = exp(x**2)*sinh(a + b*x)
    F = ((Integer(-1) * (Integer(4))**(Integer(-1))) * (sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * ((Integer(-1) * Symbol('b')) + (Integer(2) * x))))) + ((Integer(4))**(Integer(-1)) * (sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * (Integer(4))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Integer(2))**(Integer(-1)) * (Symbol('b') + (Integer(2) * x)))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_339():
    f = exp(x**2)*sinh(a + c*x**2)
    F = (Integer(-1) * ((sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c')))) * x))) * (((sympy.E)**(Symbol('a')) * (Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c')))))))**(Integer(-1)))) + (((sympy.E)**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt((Integer(1) + Symbol('c'))) * x))) * ((Integer(4) * sympy.sqrt((Integer(1) + Symbol('c')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_340():
    f = exp(x**2)*sinh(a + b*x + c*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('a')) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * (Integer(1) + (Integer(-1) * Symbol('c')))))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(-1) * (Integer(2) * (Integer(1) + (Integer(-1) * Symbol('c'))) * x))) * ((Integer(2) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c'))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * Symbol('c'))))))**(Integer(-1))) + (((sympy.E)**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * (Integer(1) + Symbol('c'))))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('b') + (Integer(2) * (Integer(1) + Symbol('c')) * x)) * ((Integer(2) * sympy.sqrt((Integer(1) + Symbol('c')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Integer(1) + Symbol('c')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_341():
    f = f**(a + b*x)*sinh(d + f*x**2)
    F = ((Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))))) * (Integer(4))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_342():
    f = f**(a + b*x)*sinh(d + f*x**2)**2
    F = (((sympy.E)**(((Integer(-2) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(8) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(4) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(8) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(4) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_343():
    f = f**(a + b*x)*sinh(d + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(12) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(6) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(12) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(6) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_344():
    f = f**(a + b*x)*sinh(d + e*x + f*x**2)
    F = ((Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1)))))) * (Integer(4))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(4))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_345():
    f = f**(a + b*x)*sinh(d + e*x + f*x**2)**2
    F = (((sympy.E)**(((Integer(-2) * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(8) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(2) * Symbol('e')) + (Integer(4) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(8) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Integer(4) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(8))**(Integer(-1))) + (Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Symbol('b') * x))) * ((Integer(2) * Symbol('b') * sympy.log(Symbol('f'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_346():
    f = f**(a + b*x)*sinh(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * Symbol('d')) + ((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(12) * Symbol('f')))**(Integer(-1))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erf')((((Integer(3) * Symbol('e')) + (Integer(6) * Symbol('f') * x) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(12) * Symbol('f')))**(Integer(-1)))))) * (Symbol('f'))**(((Integer(-1) * (Integer(2))**(Integer(-1))) + Symbol('a'))) * sympy.sqrt((sympy.pi * (Integer(3))**(Integer(-1)))) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(6) * Symbol('f') * x) + (Symbol('b') * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Integer(3)) * sympy.sqrt(Symbol('f'))))**(Integer(-1))))) * (Integer(16))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_347():
    f = f**(a + c*x**2)*sinh(d + e*x)
    F = (((sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_348():
    f = f**(a + c*x**2)*sinh(d + e*x)**2
    F = ((Integer(-1) * ((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f'))))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Symbol('c') * x * sympy.log(Symbol('f'))))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('c') * x * sympy.log(Symbol('f')))) * ((sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_349():
    f = f**(a + c*x**2)*sinh(d + e*x)**3
    F = ((Integer(-3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_350():
    f = f**(a + c*x**2)*sinh(d + f*x**2)
    F = ((Integer(-1) * ((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))))) * ((Integer(4) * (sympy.E)**(Symbol('d')) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(Symbol('d')) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_351():
    f = f**(a + c*x**2)*sinh(d + f*x**2)**2
    F = ((Integer(-1) * ((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f'))))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(8) * (sympy.E)**((Integer(2) * Symbol('d'))) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Integer(2) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_352():
    f = f**(a + c*x**2)*sinh(d + f*x**2)**3
    F = ((Integer(3) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(16) * (sympy.E)**(Symbol('d')) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((x * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))) * ((Integer(16) * (sympy.E)**((Integer(3) * Symbol('d'))) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**(Symbol('d')) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))) + (((sympy.E)**((Integer(3) * Symbol('d'))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((x * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_353():
    f = f**(a + c*x**2)*sinh(d + e*x + f*x**2)
    F = ((Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_354():
    f = f**(a + c*x**2)*sinh(d + e*x + f*x**2)**2
    F = ((Integer(-1) * ((Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((sympy.sqrt(Symbol('c')) * x * sympy.sqrt(sympy.log(Symbol('f'))))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (x * ((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * (sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * (((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (x * ((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * (sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_355():
    f = f**(a + c*x**2)*sinh(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + ((Symbol('e'))**(Integer(2)) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * Symbol('d')) + ((Integer(9) * (Symbol('e'))**(Integer(2))) * (((Integer(12) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * Symbol('e')) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * ((Symbol('e'))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((Integer(9) * (Symbol('e'))**(Integer(2))) * ((Integer(4) * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_356():
    f = f**(a + b*x + c*x**2)*sinh(d + e*x)
    F = (((sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_357():
    f = f**(a + b*x + c*x**2)*sinh(d + e*x)**2
    F = ((Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_358():
    f = f**(a + b*x + c*x**2)*sinh(d + e*x)**3
    F = ((Integer(-3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (Integer(-1) * (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(-1) * (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f'))))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * Symbol('c') * x * sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_359():
    f = f**(a + b*x + c*x**2)*sinh(d + f*x**2)
    F = (((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_360():
    f = f**(a + b*x + c*x**2)*sinh(d + f*x**2)**2
    F = ((Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-2) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(8) * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_361():
    f = f**(a + b*x + c*x**2)*sinh(d + f*x**2)**3
    F = ((Integer(-3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(4) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(-3) * Symbol('d')) + (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * (((Integer(12) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(-1) * (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * (sympy.log(Symbol('f')))**(Integer(2))) * ((Integer(4) * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_362():
    f = f**(a + b*x + c*x**2)*sinh(d + e*x + f*x**2)
    F = ((Integer(-1) * ((sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(4) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_363():
    f = f**(a + b*x + c*x**2)*sinh(d + e*x + f*x**2)**2
    F = ((Integer(-1) * ((Symbol('f'))**((Symbol('a') + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(4) * Symbol('c')))**(Integer(-1)))))) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Symbol('b') + (Integer(2) * Symbol('c') * x)) * sympy.sqrt(sympy.log(Symbol('f')))) * ((Integer(2) * sympy.sqrt(Symbol('c'))))**(Integer(-1)))))) * ((Integer(4) * sympy.sqrt(Symbol('c')) * sympy.sqrt(sympy.log(Symbol('f')))))**(Integer(-1))) + (((sympy.E)**(((Integer(-2) * Symbol('d')) + ((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(8) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(2) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (((sympy.E)**(((Integer(2) * Symbol('d')) + (Integer(-1) * ((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * (((Integer(8) * Symbol('f')) + (Integer(4) * Symbol('c') * sympy.log(Symbol('f')))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(2) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(((Integer(2) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_364():
    f = f**(a + b*x + c*x**2)*sinh(d + e*x + f*x**2)**3
    F = ((Integer(3) * (sympy.E)**(((Integer(-1) * Symbol('d')) + (((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')(((Symbol('e') + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * (Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))) + (Integer(-1) * (((sympy.E)**(((Integer(-3) * Symbol('d')) + ((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f'))))))**(Integer(2)) * (((Integer(12) * Symbol('f')) + (Integer(-1) * (Integer(4) * Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erf')((((Integer(3) * Symbol('e')) + (Integer(-1) * (Symbol('b') * sympy.log(Symbol('f')))) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Integer(-1) * (Symbol('c') * sympy.log(Symbol('f'))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (sympy.E)**((Symbol('d') + (Integer(-1) * (((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')(((Symbol('e') + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * (Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt((Symbol('f') + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))) + (((sympy.E)**(((Integer(3) * Symbol('d')) + (Integer(-1) * ((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f')))))**(Integer(2)) * ((Integer(4) * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f'))))))**(Integer(-1)))))) * (Symbol('f'))**(Symbol('a')) * sympy.sqrt(sympy.pi) * sympy.Function('Erfi')((((Integer(3) * Symbol('e')) + (Symbol('b') * sympy.log(Symbol('f'))) + (Integer(2) * x * ((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))) * ((Integer(2) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1))))) * ((Integer(16) * sympy.sqrt(((Integer(3) * Symbol('f')) + (Symbol('c') * sympy.log(Symbol('f')))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_365():
    f = (x + sinh(x))**2
    F = x**3/3 + 2*x*cosh(x) - x/2 + sinh(x)*cosh(x)/2 - 2*sinh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_366():
    f = (x + sinh(x))**3
    F = x**4/4 + 3*x**2*cosh(x) - 3*x**2/4 + 3*x*sinh(x)*cosh(x)/2 - 6*x*sinh(x) - 3*sinh(x)**2/4 + cosh(x)**3/3 + 5*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_367():
    f = sinh(a + b*x)/(c + d*x**2)
    F = (Integer(-1) * ((sympy.Function('CoshIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + ((sympy.Function('CoshIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x)))) * sympy.sinh((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))))) * sympy.Function('SinhIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Integer(-1) * (Symbol('b') * x))))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1)))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * sympy.sqrt((Integer(-1) * Symbol('c')))) * (sympy.sqrt(Symbol('d')))**(Integer(-1))) + (Symbol('b') * x)))) * ((Integer(2) * sympy.sqrt((Integer(-1) * Symbol('c'))) * sympy.sqrt(Symbol('d'))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_1_Hyperbolic_sine_6_1_5_Hyperbolic_sine_functions_368():
    f = sinh(a + b*x)/(c + d*x + e*x**2)
    F = ((sympy.Function('CoshIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))))))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.Function('CoshIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x))) * sympy.sinh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))))))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1)))) + ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * (Symbol('d') + (Integer(-1) * sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e')))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))) + (Integer(-1) * ((sympy.cosh((Symbol('a') + (Integer(-1) * ((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1)))))) * sympy.Function('SinhIntegral')((((Symbol('b') * (Symbol('d') + sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))) * ((Integer(2) * Symbol('e')))**(Integer(-1))) + (Symbol('b') * x)))) * (sympy.sqrt(((Symbol('d'))**(Integer(2)) + (Integer(-1) * (Integer(4) * Symbol('c') * Symbol('e'))))))**(Integer(-1))))
    assert integrate(f, x) == F

