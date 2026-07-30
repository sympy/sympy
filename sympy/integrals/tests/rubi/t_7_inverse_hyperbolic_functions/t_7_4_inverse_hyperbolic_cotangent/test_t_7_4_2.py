"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.4 Inverse hyperbolic cotangent/7.4.2 Exponentials of inverse hyperbolic cotangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, c, m, n, p = symbols('a c m n p')

def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_1():
    f = x**3*exp(acoth(a*x))
    F = x**4*sqrt(1 - 1/(a**2*x**2))/4 + x**3*sqrt(1 - 1/(a**2*x**2))/(3*a) + 3*x**2*sqrt(1 - 1/(a**2*x**2))/(8*a**2) + 2*x*sqrt(1 - 1/(a**2*x**2))/(3*a**3) + 3*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_2():
    f = x**2*exp(acoth(a*x))
    F = x**3*sqrt(1 - 1/(a**2*x**2))/3 + x**2*sqrt(1 - 1/(a**2*x**2))/(2*a) + 2*x*sqrt(1 - 1/(a**2*x**2))/(3*a**2) + atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_3():
    f = x*exp(acoth(a*x))
    F = x**2*sqrt(1 - 1/(a**2*x**2))/2 + x*sqrt(1 - 1/(a**2*x**2))/a + atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_4():
    f = exp(acoth(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2)) + atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_5():
    f = exp(acoth(a*x))/x
    F = -acsc(a*x) + atanh(sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_6():
    f = exp(acoth(a*x))/x**2
    F = a*sqrt(1 - 1/(a**2*x**2)) - a*acsc(a*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_7():
    f = exp(acoth(a*x))/x**3
    F = -a**2*acsc(a*x)/2 + a*sqrt(1 - 1/(a**2*x**2))*(2*a + 1/x)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_8():
    f = exp(acoth(a*x))/x**4
    F = -a**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/3 + a**3*sqrt(1 - 1/(a**2*x**2)) - a**3*acsc(a*x)/2 + a**2*sqrt(1 - 1/(a**2*x**2))/(2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_9():
    f = exp(acoth(a*x))/x**5
    F = -3*a**4*acsc(a*x)/8 + a**3*sqrt(1 - 1/(a**2*x**2))*(16*a + 9/x)/24 + a**2*sqrt(1 - 1/(a**2*x**2))/(3*x**2) + a*sqrt(1 - 1/(a**2*x**2))/(4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_10():
    f = x**3*exp(2*acoth(a*x))
    F = x**4/4 + 2*x**3/(3*a) + x**2/a**2 + 2*x/a**3 + 2*log(-a*x + 1)/a**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_11():
    f = x**2*exp(2*acoth(a*x))
    F = x**3/3 + x**2/a + 2*x/a**2 + 2*log(-a*x + 1)/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_12():
    f = x*exp(2*acoth(a*x))
    F = x**2/2 + 2*x/a + 2*log(-a*x + 1)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_13():
    f = exp(2*acoth(a*x))
    F = x + 2*log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_14():
    f = exp(2*acoth(a*x))/x
    F = -log(x) + 2*log(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_15():
    f = exp(2*acoth(a*x))/x**2
    F = -2*a*log(x) + 2*a*log(-a*x + 1) + 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_16():
    f = exp(2*acoth(a*x))/x**3
    F = -2*a**2*log(x) + 2*a**2*log(-a*x + 1) + 2*a/x + 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_17():
    f = exp(2*acoth(a*x))/x**4
    F = -2*a**3*log(x) + 2*a**3*log(-a*x + 1) + 2*a**2/x + a/x**2 + 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_18():
    f = x**2*exp(3*acoth(a*x))
    F = x**3*sqrt(1 - 1/(a**2*x**2))/3 + 3*x**2*sqrt(1 - 1/(a**2*x**2))/(2*a) + 14*x*sqrt(1 - 1/(a**2*x**2))/(3*a**2) - 4*sqrt(1 - 1/(a**2*x**2))/(a**2*(a - 1/x)) + 11*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_19():
    f = x*exp(3*acoth(a*x))
    F = x**2*sqrt(1 - 1/(a**2*x**2))/2 + 3*x*sqrt(1 - 1/(a**2*x**2))/a - 4*sqrt(1 - 1/(a**2*x**2))/(a*(a - 1/x)) + 9*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_20():
    f = exp(3*acoth(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2)) - 4*sqrt(1 - 1/(a**2*x**2))/(a - 1/x) + 3*atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_21():
    f = exp(3*acoth(a*x))/x
    F = -4*a*sqrt(1 - 1/(a**2*x**2))/(a - 1/x) + acsc(a*x) + atanh(sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_22():
    f = exp(3*acoth(a*x))/x**2
    F = -3*a*sqrt(1 - 1/(a**2*x**2)) + 3*a*acsc(a*x) - 2*(a + 1/x)**2/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_23():
    f = exp(3*acoth(a*x))/x**3
    F = -a**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(a - 1/x)**3 - 3*a**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*a - 2/x) - 9*a**2*sqrt(1 - 1/(a**2*x**2))/2 + 9*a**2*acsc(a*x)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_24():
    f = exp(3*acoth(a*x))/x**4
    F = 11*a**3*acsc(a*x)/2 - a**2*sqrt(1 - 1/(a**2*x**2))*(28*a + 3/x)/6 - a*sqrt(1 - 1/(a**2*x**2))*(3*a + 1/x)**2/3 - (a + 1/x)**3/sqrt(1 - 1/(a**2*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_25():
    f = x**3*exp(4*acoth(a*x))
    F = x**4/4 + 4*x**3/(3*a) + 4*x**2/a**2 + 12*x/a**3 + 16*log(-a*x + 1)/a**4 + 4/(a**4*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_26():
    f = x**2*exp(4*acoth(a*x))
    F = x**3/3 + 2*x**2/a + 8*x/a**2 + 12*log(-a*x + 1)/a**3 + 4/(a**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_27():
    f = x*exp(4*acoth(a*x))
    F = x**2/2 + 4*x/a + 8*log(-a*x + 1)/a**2 + 4/(a**2*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_28():
    f = exp(4*acoth(a*x))
    F = x + 4*log(-a*x + 1)/a + 4/(a*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_29():
    f = exp(4*acoth(a*x))/x
    F = log(x) + 4/(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_30():
    f = exp(4*acoth(a*x))/x**2
    F = 4*a*log(x) - 4*a*log(-a*x + 1) + 4*a/(-a*x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_31():
    f = exp(4*acoth(a*x))/x**3
    F = 8*a**2*log(x) - 8*a**2*log(-a*x + 1) + 4*a**2/(-a*x + 1) - 4*a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_32():
    f = exp(4*acoth(a*x))/x**4
    F = 12*a**3*log(x) - 12*a**3*log(-a*x + 1) + 4*a**3/(-a*x + 1) - 8*a**2/x - 2*a/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_33():
    f = x**3*exp(-acoth(a*x))
    F = x**4*sqrt(1 - 1/(a**2*x**2))/4 - x**3*sqrt(1 - 1/(a**2*x**2))/(3*a) + 3*x**2*sqrt(1 - 1/(a**2*x**2))/(8*a**2) - 2*x*sqrt(1 - 1/(a**2*x**2))/(3*a**3) + 3*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_34():
    f = x**2*exp(-acoth(a*x))
    F = x**3*sqrt(1 - 1/(a**2*x**2))/3 - x**2*sqrt(1 - 1/(a**2*x**2))/(2*a) + 2*x*sqrt(1 - 1/(a**2*x**2))/(3*a**2) - atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_35():
    f = x*exp(-acoth(a*x))
    F = x**2*sqrt(1 - 1/(a**2*x**2))/2 - x*sqrt(1 - 1/(a**2*x**2))/a + atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_36():
    f = exp(-acoth(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2)) - atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_37():
    f = exp(-acoth(a*x))/x
    F = acsc(a*x) + atanh(sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_38():
    f = exp(-acoth(a*x))/x**2
    F = -a*sqrt(1 - 1/(a**2*x**2)) - a*acsc(a*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_39():
    f = exp(-acoth(a*x))/x**3
    F = a**2*acsc(a*x)/2 + a*sqrt(1 - 1/(a**2*x**2))*(2*a - 1/x)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_40():
    f = exp(-acoth(a*x))/x**4
    F = a**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/3 - a**3*sqrt(1 - 1/(a**2*x**2)) - a**3*acsc(a*x)/2 + a**2*sqrt(1 - 1/(a**2*x**2))/(2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_41():
    f = exp(-acoth(a*x))/x**5
    F = 3*a**4*acsc(a*x)/8 + a**3*sqrt(1 - 1/(a**2*x**2))*(16*a - 9/x)/24 + a**2*sqrt(1 - 1/(a**2*x**2))/(3*x**2) - a*sqrt(1 - 1/(a**2*x**2))/(4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_42():
    f = x**3*exp(-2*acoth(a*x))
    F = x**4/4 - 2*x**3/(3*a) + x**2/a**2 - 2*x/a**3 + 2*log(a*x + 1)/a**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_43():
    f = x**2*exp(-2*acoth(a*x))
    F = x**3/3 - x**2/a + 2*x/a**2 - 2*log(a*x + 1)/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_44():
    f = x*exp(-2*acoth(a*x))
    F = x**2/2 - 2*x/a + 2*log(a*x + 1)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_45():
    f = exp(-2*acoth(a*x))
    F = x - 2*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_46():
    f = exp(-2*acoth(a*x))/x
    F = -log(x) + 2*log(a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_47():
    f = exp(-2*acoth(a*x))/x**2
    F = 2*a*log(x) - 2*a*log(a*x + 1) + 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_48():
    f = exp(-2*acoth(a*x))/x**3
    F = -2*a**2*log(x) + 2*a**2*log(a*x + 1) - 2*a/x + 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_49():
    f = exp(-2*acoth(a*x))/x**4
    F = 2*a**3*log(x) - 2*a**3*log(a*x + 1) + 2*a**2/x - a/x**2 + 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_50():
    f = x**3*exp(-3*acoth(a*x))
    F = x**4*sqrt(1 - 1/(a**2*x**2))/4 - x**3*sqrt(1 - 1/(a**2*x**2))/a + 19*x**2*sqrt(1 - 1/(a**2*x**2))/(8*a**2) - 6*x*sqrt(1 - 1/(a**2*x**2))/a**3 - 4*sqrt(1 - 1/(a**2*x**2))/(a**3*(a + 1/x)) + 51*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_51():
    f = x**2*exp(-3*acoth(a*x))
    F = x**3*sqrt(1 - 1/(a**2*x**2))/3 - 3*x**2*sqrt(1 - 1/(a**2*x**2))/(2*a) + 14*x*sqrt(1 - 1/(a**2*x**2))/(3*a**2) + 4*sqrt(1 - 1/(a**2*x**2))/(a**2*(a + 1/x)) - 11*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_52():
    f = x*exp(-3*acoth(a*x))
    F = x**2*sqrt(1 - 1/(a**2*x**2))/2 - 3*x*sqrt(1 - 1/(a**2*x**2))/a - 4*sqrt(1 - 1/(a**2*x**2))/(a*(a + 1/x)) + 9*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_53():
    f = exp(-3*acoth(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2)) + 4*sqrt(1 - 1/(a**2*x**2))/(a + 1/x) - 3*atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_54():
    f = exp(-3*acoth(a*x))/x
    F = -4*a*sqrt(1 - 1/(a**2*x**2))/(a + 1/x) - acsc(a*x) + atanh(sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_55():
    f = exp(-3*acoth(a*x))/x**2
    F = 3*a*sqrt(1 - 1/(a**2*x**2)) + 3*a*acsc(a*x) + 2*(a - 1/x)**2/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_56():
    f = exp(-3*acoth(a*x))/x**3
    F = -a**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(a + 1/x)**3 - 3*a**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*a + 2/x) - 9*a**2*sqrt(1 - 1/(a**2*x**2))/2 - 9*a**2*acsc(a*x)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_57():
    f = exp(-3*acoth(a*x))/x**4
    F = 11*a**3*acsc(a*x)/2 + a**2*sqrt(1 - 1/(a**2*x**2))*(28*a - 3/x)/6 + a*sqrt(1 - 1/(a**2*x**2))*(3*a - 1/x)**2/3 + (a - 1/x)**3/sqrt(1 - 1/(a**2*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_58():
    f = exp(-3*acoth(a*x))/x**5
    F = -27*a**4*sqrt(1 - 1/(a**2*x**2))/4 - 51*a**4*acsc(a*x)/8 - 9*a**3*sqrt(1 - 1/(a**2*x**2))*(2*a - 3/x)/8 - a**2*sqrt(1 - 1/(a**2*x**2))/x**2 - a*(a - 1/x)**3/sqrt(1 - 1/(a**2*x**2)) + a*sqrt(1 - 1/(a**2*x**2))/(4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_59():
    f = x**4*exp(acoth(a*x)/2)
    F = x**5*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/5 + 9*x**4*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(40*a) + 11*x**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(48*a**2) + 269*x**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(960*a**3) + 611*x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(1920*a**4) + 31*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5) + 31*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_60():
    f = x**3*exp(acoth(a*x)/2)
    F = x**4*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/4 + 7*x**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(24*a) + 29*x**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(96*a**2) + 83*x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(192*a**3) + 11*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4) + 11*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_61():
    f = x**2*exp(acoth(a*x)/2)
    F = x**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/3 + 5*x**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(12*a) + 11*x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(24*a**2) + 3*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3) + 3*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_62():
    f = x*exp(acoth(a*x)/2)
    F = x**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(5)/4)/2 + x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(4*a) + atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2) + atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_63():
    f = exp(acoth(a*x)/2)
    F = x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4) + atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a + atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_64():
    f = exp(acoth(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 - sqrt(2)*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 + 2*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1) + 2*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_65():
    f = exp(acoth(a*x)/2)/x**2
    F = a*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(2)*a*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 - sqrt(2)*a*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 + sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/2 + sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_66():
    f = exp(acoth(a*x)/2)/x**3
    F = a**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(5)/4)/2 + a**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/4 + sqrt(2)*a**2*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 - sqrt(2)*a**2*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 + sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/8 + sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_67():
    f = exp(acoth(a*x)/2)/x**4
    F = a**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(5)/4)/12 + 3*a**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/8 + 3*sqrt(2)*a**3*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 - 3*sqrt(2)*a**3*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 + 3*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/16 + 3*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/16 + a**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(5)/4)/(3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_68():
    f = x**4*exp(3*acoth(a*x)/2)
    F = x**5*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/5 + 11*x**4*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(40*a) + 5*x**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(16*a**2) + 157*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(320*a**3) + 557*x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(640*a**4) - 237*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5) + 237*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_69():
    f = x**3*exp(3*acoth(a*x)/2)
    F = x**4*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/4 + 3*x**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(8*a) + 15*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(32*a**2) + 63*x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(64*a**3) - 123*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4) + 123*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_70():
    f = x**2*exp(3*acoth(a*x)/2)
    F = x**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/3 + 7*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(12*a) + 23*x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(24*a**2) - 17*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3) + 17*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_71():
    f = x*exp(3*acoth(a*x)/2)
    F = x**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(7)/4)/2 + 3*x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(4*a) - 9*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2) + 9*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_72():
    f = exp(3*acoth(a*x)/2)
    F = x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4) - 3*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a + 3*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_73():
    f = exp(3*acoth(a*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 + sqrt(2)*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 - 2*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1) + 2*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_74():
    f = exp(3*acoth(a*x)/2)/x**2
    F = a*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4) - 3*sqrt(2)*a*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 + 3*sqrt(2)*a*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 + 3*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/2 + 3*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_75():
    f = exp(3*acoth(a*x)/2)/x**3
    F = a**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(7)/4)/2 + 3*a**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/4 - 9*sqrt(2)*a**2*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 + 9*sqrt(2)*a**2*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 + 9*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/8 + 9*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_76():
    f = exp(3*acoth(a*x)/2)/x**4
    F = a**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(7)/4)/4 + 17*a**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/24 - 17*sqrt(2)*a**3*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 + 17*sqrt(2)*a**3*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 + 17*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/16 + 17*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/16 + a**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(7)/4)/(3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_77():
    f = x**4*exp(5*acoth(a*x)/2)
    F = x**5*(1 + 1/(a*x))**(sympy.S(1)/4)/(5*(1 - 1/(a*x))**(sympy.S(1)/4)) + 21*x**4*(1 + 1/(a*x))**(sympy.S(1)/4)/(40*a*(1 - 1/(a*x))**(sympy.S(1)/4)) + 181*x**3*(1 + 1/(a*x))**(sympy.S(1)/4)/(240*a**2*(1 - 1/(a*x))**(sympy.S(1)/4)) + 1189*x**2*(1 + 1/(a*x))**(sympy.S(1)/4)/(960*a**3*(1 - 1/(a*x))**(sympy.S(1)/4)) + 5533*x*(1 + 1/(a*x))**(sympy.S(1)/4)/(1920*a**4*(1 - 1/(a*x))**(sympy.S(1)/4)) + 1003*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5) + 1003*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5) - 26111*(1 + 1/(a*x))**(sympy.S(1)/4)/(1920*a**5*(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_78():
    f = x**3*exp(5*acoth(a*x)/2)
    F = x**4*(1 + 1/(a*x))**(sympy.S(1)/4)/(4*(1 - 1/(a*x))**(sympy.S(1)/4)) + 17*x**3*(1 + 1/(a*x))**(sympy.S(1)/4)/(24*a*(1 - 1/(a*x))**(sympy.S(1)/4)) + 113*x**2*(1 + 1/(a*x))**(sympy.S(1)/4)/(96*a**2*(1 - 1/(a*x))**(sympy.S(1)/4)) + 521*x*(1 + 1/(a*x))**(sympy.S(1)/4)/(192*a**3*(1 - 1/(a*x))**(sympy.S(1)/4)) + 475*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4) + 475*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4) - 2467*(1 + 1/(a*x))**(sympy.S(1)/4)/(192*a**4*(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_79():
    f = x**2*exp(5*acoth(a*x)/2)
    F = x**3*(1 + 1/(a*x))**(sympy.S(1)/4)/(3*(1 - 1/(a*x))**(sympy.S(1)/4)) + 13*x**2*(1 + 1/(a*x))**(sympy.S(1)/4)/(12*a*(1 - 1/(a*x))**(sympy.S(1)/4)) + 61*x*(1 + 1/(a*x))**(sympy.S(1)/4)/(24*a**2*(1 - 1/(a*x))**(sympy.S(1)/4)) + 55*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3) + 55*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3) - 287*(1 + 1/(a*x))**(sympy.S(1)/4)/(24*a**3*(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_80():
    f = x*exp(5*acoth(a*x)/2)
    F = x**2*(1 + 1/(a*x))**(sympy.S(9)/4)/(2*(1 - 1/(a*x))**(sympy.S(1)/4)) + 5*x*(1 + 1/(a*x))**(sympy.S(5)/4)/(4*a*(1 - 1/(a*x))**(sympy.S(1)/4)) + 25*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2) + 25*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2) - 25*(1 + 1/(a*x))**(sympy.S(1)/4)/(2*a**2*(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_81():
    f = exp(5*acoth(a*x)/2)
    F = x*(1 + 1/(a*x))**(sympy.S(5)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) + 5*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a + 5*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a - 10*(1 + 1/(a*x))**(sympy.S(1)/4)/(a*(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_82():
    f = exp(5*acoth(a*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 + sqrt(2)*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 + 2*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1) + 2*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)) - 8*(1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_83():
    f = exp(5*acoth(a*x)/2)/x**2
    F = -5*a*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4) - 5*sqrt(2)*a*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 + 5*sqrt(2)*a*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 - 5*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/2 - 5*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2 - 4*a*(1 + 1/(a*x))**(sympy.S(5)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_84():
    f = exp(5*acoth(a*x)/2)/x**3
    F = -5*a**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(5)/4)/2 - 25*a**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/4 - 25*sqrt(2)*a**2*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 + 25*sqrt(2)*a**2*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 - 25*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/8 - 25*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8 - 2*a**2*(1 + 1/(a*x))**(sympy.S(9)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_85():
    f = exp(5*acoth(a*x)/2)/x**4
    F = -a**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(9)/4)/3 - 11*a**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(5)/4)/4 - 55*a**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/8 - 55*sqrt(2)*a**3*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 + 55*sqrt(2)*a**3*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 - 55*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/16 - 55*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/16 - 2*a**3*(1 + 1/(a*x))**(sympy.S(9)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_86():
    f = x**4*exp(-acoth(a*x)/2)
    F = x**5*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/5 - 9*x**4*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(40*a) + 11*x**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(48*a**2) - 269*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(960*a**3) + 611*x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(1920*a**4) + 31*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5) - 31*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_87():
    f = x**3*exp(-acoth(a*x)/2)
    F = x**4*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/4 - 7*x**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(24*a) + 29*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(96*a**2) - 83*x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(192*a**3) - 11*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4) + 11*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_88():
    f = x**2*exp(-acoth(a*x)/2)
    F = x**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/3 - 5*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(12*a) + 11*x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(24*a**2) + 3*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3) - 3*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_89():
    f = x*exp(-acoth(a*x)/2)
    F = x**2*(1 - 1/(a*x))**(sympy.S(5)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/2 - x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(4*a) - atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2) + atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_90():
    f = exp(-acoth(a*x)/2)
    F = x*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4) + atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a - atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_91():
    f = exp(-acoth(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 - sqrt(2)*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 - 2*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1) + 2*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_92():
    f = exp(-acoth(a*x)/2)/x**2
    F = -a*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4) - sqrt(2)*a*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 + sqrt(2)*a*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 + sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/2 + sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_93():
    f = exp(-acoth(a*x)/2)/x**3
    F = a**2*(1 - 1/(a*x))**(sympy.S(5)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/2 + a**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/4 + sqrt(2)*a**2*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 - sqrt(2)*a**2*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 - sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/8 - sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_94():
    f = exp(-acoth(a*x)/2)/x**4
    F = -a**3*(1 - 1/(a*x))**(sympy.S(5)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/12 - 3*a**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/8 - 3*sqrt(2)*a**3*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 + 3*sqrt(2)*a**3*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 + 3*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/16 + 3*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/16 + a**2*(1 - 1/(a*x))**(sympy.S(5)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/(3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_95():
    f = x**4*exp(-3*acoth(a*x)/2)
    F = x**5*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/5 - 11*x**4*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(40*a) + 5*x**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(16*a**2) - 157*x**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(320*a**3) + 557*x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(640*a**4) - 237*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5) - 237*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_96():
    f = x**3*exp(-3*acoth(a*x)/2)
    F = x**4*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/4 - 3*x**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(8*a) + 15*x**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(32*a**2) - 63*x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(64*a**3) + 123*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4) + 123*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_97():
    f = x**2*exp(-3*acoth(a*x)/2)
    F = x**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/3 - 7*x**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(12*a) + 23*x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(24*a**2) - 17*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3) - 17*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_98():
    f = x*exp(-3*acoth(a*x)/2)
    F = x**2*(1 - 1/(a*x))**(sympy.S(7)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/2 - 3*x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(4*a) + 9*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2) + 9*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_99():
    f = exp(-3*acoth(a*x)/2)
    F = x*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4) - 3*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a - 3*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_100():
    f = exp(-3*acoth(a*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 + sqrt(2)*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 + 2*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1) + 2*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_101():
    f = exp(-3*acoth(a*x)/2)/x**2
    F = -a*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4) + 3*sqrt(2)*a*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 - 3*sqrt(2)*a*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 + 3*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/2 + 3*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_102():
    f = exp(-3*acoth(a*x)/2)/x**3
    F = a**2*(1 - 1/(a*x))**(sympy.S(7)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/2 + 3*a**2*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/4 - 9*sqrt(2)*a**2*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 + 9*sqrt(2)*a**2*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 - 9*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/8 - 9*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_103():
    f = exp(-3*acoth(a*x)/2)/x**4
    F = -a**3*(1 - 1/(a*x))**(sympy.S(7)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/4 - 17*a**3*(1 - 1/(a*x))**(sympy.S(3)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/24 + 17*sqrt(2)*a**3*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 - 17*sqrt(2)*a**3*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 + 17*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/16 + 17*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/16 + a**2*(1 - 1/(a*x))**(sympy.S(7)/4)*(1 + 1/(a*x))**(sympy.S(1)/4)/(3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_104():
    f = x**4*exp(-5*acoth(a*x)/2)
    F = x**5*(1 - 1/(a*x))**(sympy.S(1)/4)/(5*(1 + 1/(a*x))**(sympy.S(1)/4)) - 21*x**4*(1 - 1/(a*x))**(sympy.S(1)/4)/(40*a*(1 + 1/(a*x))**(sympy.S(1)/4)) + 181*x**3*(1 - 1/(a*x))**(sympy.S(1)/4)/(240*a**2*(1 + 1/(a*x))**(sympy.S(1)/4)) - 1189*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)/(960*a**3*(1 + 1/(a*x))**(sympy.S(1)/4)) + 5533*x*(1 - 1/(a*x))**(sympy.S(1)/4)/(1920*a**4*(1 + 1/(a*x))**(sympy.S(1)/4)) + 26111*(1 - 1/(a*x))**(sympy.S(1)/4)/(1920*a**5*(1 + 1/(a*x))**(sympy.S(1)/4)) + 1003*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5) - 1003*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(128*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_105():
    f = x**3*exp(-5*acoth(a*x)/2)
    F = x**4*(1 - 1/(a*x))**(sympy.S(1)/4)/(4*(1 + 1/(a*x))**(sympy.S(1)/4)) - 17*x**3*(1 - 1/(a*x))**(sympy.S(1)/4)/(24*a*(1 + 1/(a*x))**(sympy.S(1)/4)) + 113*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)/(96*a**2*(1 + 1/(a*x))**(sympy.S(1)/4)) - 521*x*(1 - 1/(a*x))**(sympy.S(1)/4)/(192*a**3*(1 + 1/(a*x))**(sympy.S(1)/4)) - 2467*(1 - 1/(a*x))**(sympy.S(1)/4)/(192*a**4*(1 + 1/(a*x))**(sympy.S(1)/4)) - 475*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4) + 475*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_106():
    f = x**2*exp(-5*acoth(a*x)/2)
    F = x**3*(1 - 1/(a*x))**(sympy.S(1)/4)/(3*(1 + 1/(a*x))**(sympy.S(1)/4)) - 13*x**2*(1 - 1/(a*x))**(sympy.S(1)/4)/(12*a*(1 + 1/(a*x))**(sympy.S(1)/4)) + 61*x*(1 - 1/(a*x))**(sympy.S(1)/4)/(24*a**2*(1 + 1/(a*x))**(sympy.S(1)/4)) + 287*(1 - 1/(a*x))**(sympy.S(1)/4)/(24*a**3*(1 + 1/(a*x))**(sympy.S(1)/4)) + 55*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3) - 55*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_107():
    f = x*exp(-5*acoth(a*x)/2)
    F = x**2*(1 - 1/(a*x))**(sympy.S(9)/4)/(2*(1 + 1/(a*x))**(sympy.S(1)/4)) - 5*x*(1 - 1/(a*x))**(sympy.S(5)/4)/(4*a*(1 + 1/(a*x))**(sympy.S(1)/4)) - 25*(1 - 1/(a*x))**(sympy.S(1)/4)/(2*a**2*(1 + 1/(a*x))**(sympy.S(1)/4)) - 25*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2) + 25*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_108():
    f = exp(-5*acoth(a*x)/2)
    F = x*(1 - 1/(a*x))**(sympy.S(5)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 10*(1 - 1/(a*x))**(sympy.S(1)/4)/(a*(1 + 1/(a*x))**(sympy.S(1)/4)) + 5*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a - 5*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_109():
    f = exp(-5*acoth(a*x)/2)/x
    F = -8*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - sqrt(2)*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 + sqrt(2)*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/2 - 2*atan((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1) + 2*atanh((1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_110():
    f = exp(-5*acoth(a*x)/2)/x**2
    F = 4*a*(1 - 1/(a*x))**(sympy.S(5)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 5*a*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4) + 5*sqrt(2)*a*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 - 5*sqrt(2)*a*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/4 - 5*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/2 - 5*sqrt(2)*a*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_111():
    f = exp(-5*acoth(a*x)/2)/x**3
    F = -2*a**2*(1 - 1/(a*x))**(sympy.S(9)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 5*a**2*(1 - 1/(a*x))**(sympy.S(5)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/2 - 25*a**2*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/4 - 25*sqrt(2)*a**2*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 + 25*sqrt(2)*a**2*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/16 + 25*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/8 + 25*sqrt(2)*a**2*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_112():
    f = exp(-5*acoth(a*x)/2)/x**4
    F = a**3*(1 - 1/(a*x))**(sympy.S(9)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/3 + 2*a**3*(1 - 1/(a*x))**(sympy.S(9)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 11*a**3*(1 - 1/(a*x))**(sympy.S(5)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/4 + 55*a**3*(1 - 1/(a*x))**(sympy.S(1)/4)*(1 + 1/(a*x))**(sympy.S(3)/4)/8 + 55*sqrt(2)*a**3*log(-sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 - 55*sqrt(2)*a**3*log(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + sqrt(1 - 1/(a*x))/sqrt(1 + 1/(a*x)) + 1)/32 - 55*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) - 1)/16 - 55*sqrt(2)*a**3*atan(sqrt(2)*(1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/16
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_113():
    f = x**2*exp(acoth(x)/3)
    F = x**3*((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6)/3 + 7*x**2*((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6)/18 + 11*x*((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6)/27 - 19*log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) - (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/324 + 19*log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/324 - 19*sqrt(3)*atan(sqrt(3)*(1 - 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3)/162 + 19*sqrt(3)*atan(sqrt(3)*(1 + 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3)/162 + 19*atanh((1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/81
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_114():
    f = x*exp(acoth(x)/3)
    F = x**2*((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(7)/6)/2 + x*((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6)/6 - log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) - (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/36 + log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/36 - sqrt(3)*atan(sqrt(3)*(1 - 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3)/18 + sqrt(3)*atan(sqrt(3)*(1 + 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3)/18 + atanh((1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/9
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_115():
    f = exp(acoth(x)/3)
    F = x*((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6) - log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) - (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/6 + log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/6 - sqrt(3)*atan(sqrt(3)*(1 - 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3)/3 + sqrt(3)*atan(sqrt(3)*(1 + 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3)/3 + 2*atanh((1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_116():
    f = exp(acoth(x)/3)/x
    F = -log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) - (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/2 + log(1 + (1 + 1/x)**(sympy.S(1)/3)/((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/2 + sqrt(3)*log(-sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/2 - sqrt(3)*log(sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/2 - sqrt(3)*atan(sqrt(3)*(1 - 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3) + sqrt(3)*atan(sqrt(3)*(1 + 2*(1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))/3) + 2*atan(((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6)) + atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) - sqrt(3)) + atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + sqrt(3)) + 2*atanh((1 + 1/x)**(sympy.S(1)/6)/((x - 1)/x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_117():
    f = exp(acoth(x)/3)/x**2
    F = ((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6) + sqrt(3)*log(-sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/6 - sqrt(3)*log(sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/6 + 2*atan(((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6))/3 + atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) - sqrt(3))/3 + atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + sqrt(3))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_118():
    f = exp(acoth(x)/3)/x**3
    F = ((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(7)/6)/2 + ((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6)/6 + sqrt(3)*log(-sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/36 - sqrt(3)*log(sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/36 + atan(((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6))/9 + atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) - sqrt(3))/18 + atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + sqrt(3))/18
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_119():
    f = exp(acoth(x)/3)/x**4
    F = ((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(7)/6)/18 + 19*((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(1)/6)/54 + 19*sqrt(3)*log(-sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/324 - 19*sqrt(3)*log(sqrt(3)*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + ((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/324 + 19*atan(((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6))/81 + 19*atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) - sqrt(3))/162 + 19*atan(2*((x - 1)/x)**(sympy.S(1)/6)/(1 + 1/x)**(sympy.S(1)/6) + sqrt(3))/162 + ((x - 1)/x)**(sympy.S(5)/6)*(1 + 1/x)**(sympy.S(7)/6)/(3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_120():
    f = x**2*exp(2*acoth(x)/3)
    F = x**3*((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(1)/3)/3 + 4*x**2*((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(1)/3)/9 + 14*x*((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(1)/3)/27 - 11*log(x)/81 - 11*log(-((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/3))/27 - 22*sqrt(3)*atan(2*sqrt(3)*((x - 1)/x)**(sympy.S(1)/3)/(3*(1 + 1/x)**(sympy.S(1)/3)) + sqrt(3)/3)/81
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_121():
    f = x*exp(2*acoth(x)/3)
    F = x**2*((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(4)/3)/2 + x*((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(1)/3)/3 - log(x)/9 - log(-((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/3))/3 - 2*sqrt(3)*atan(2*sqrt(3)*((x - 1)/x)**(sympy.S(1)/3)/(3*(1 + 1/x)**(sympy.S(1)/3)) + sqrt(3)/3)/9
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_122():
    f = exp(2*acoth(x)/3)
    F = x*((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(1)/3) - log(x)/3 - log(-((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/3)) - 2*sqrt(3)*atan(2*sqrt(3)*((x - 1)/x)**(sympy.S(1)/3)/(3*(1 + 1/x)**(sympy.S(1)/3)) + sqrt(3)/3)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_123():
    f = exp(2*acoth(x)/3)/x
    F = -log(x)/2 - log(1 + 1/x)/2 - 3*log(-((x - 1)/x)**(sympy.S(1)/3) + (1 + 1/x)**(sympy.S(1)/3))/2 - 3*log(((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/2 + sqrt(3)*atan(2*sqrt(3)*((x - 1)/x)**(sympy.S(1)/3)/(3*(1 + 1/x)**(sympy.S(1)/3)) - sqrt(3)/3) - sqrt(3)*atan(2*sqrt(3)*((x - 1)/x)**(sympy.S(1)/3)/(3*(1 + 1/x)**(sympy.S(1)/3)) + sqrt(3)/3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_124():
    f = exp(2*acoth(x)/3)/x**2
    F = ((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(1)/3) - log(1 + 1/x)/3 - log(((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1) + 2*sqrt(3)*atan(2*sqrt(3)*((x - 1)/x)**(sympy.S(1)/3)/(3*(1 + 1/x)**(sympy.S(1)/3)) - sqrt(3)/3)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_125():
    f = exp(2*acoth(x)/3)/x**3
    F = ((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(4)/3)/2 + ((x - 1)/x)**(sympy.S(2)/3)*(1 + 1/x)**(sympy.S(1)/3)/3 - log(1 + 1/x)/9 - log(((x - 1)/x)**(sympy.S(1)/3)/(1 + 1/x)**(sympy.S(1)/3) + 1)/3 + 2*sqrt(3)*atan(2*sqrt(3)*((x - 1)/x)**(sympy.S(1)/3)/(3*(1 + 1/x)**(sympy.S(1)/3)) - sqrt(3)/3)/9
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_126():
    f = x**2*exp(acoth(a*x)/4)
    F = x**3*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(1)/8)/3 + 3*x**2*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(1)/8)/(8*a) + 37*x*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(1)/8)/(96*a**2) - 11*sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(256*a**3) + 11*sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(256*a**3) + 11*atan((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(64*a**3) - 11*sqrt(2)*atan(1 - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(128*a**3) + 11*sqrt(2)*atan(1 + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(128*a**3) + 11*atanh((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(64*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_127():
    f = x*exp(acoth(a*x)/4)
    F = x**2*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(9)/8)/2 + x*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(1)/8)/(8*a) - sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(64*a**2) + sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(64*a**2) + atan((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(16*a**2) - sqrt(2)*atan(1 - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(32*a**2) + sqrt(2)*atan(1 + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(32*a**2) + atanh((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(16*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_128():
    f = exp(acoth(a*x)/4)
    F = x*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(1)/8) - sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(8*a) + sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(8*a) + atan((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(2*a) - sqrt(2)*atan(1 - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(4*a) + sqrt(2)*atan(1 + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(4*a) + atanh((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_129():
    f = exp(acoth(a*x)/4)/x
    F = -sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/2 + sqrt(2)*log(1 + (1 + 1/(a*x))**(sympy.S(1)/4)/(1 - 1/(a*x))**(sympy.S(1)/4) + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))/2 + sqrt(2 - sqrt(2))*log(-(1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(2 - sqrt(2))/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2 - sqrt(2 - sqrt(2))*log((1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(2 - sqrt(2))/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2 + sqrt(sqrt(2) + 2)*log(-(1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2 - sqrt(sqrt(2) + 2)*log((1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/2 + 2*atan((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8)) - sqrt(2 - sqrt(2))*atan((-2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) + sqrt(2 - sqrt(2))*atan((2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) - sqrt(sqrt(2) + 2)*atan((-2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2)) + sqrt(sqrt(2) + 2)*atan((2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2)) - sqrt(2)*atan(1 - sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8)) + sqrt(2)*atan(1 + sqrt(2)*(1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8)) + 2*atanh((1 + 1/(a*x))**(sympy.S(1)/8)/(1 - 1/(a*x))**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_130():
    f = exp(acoth(a*x)/4)/x**2
    F = a*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(1)/8) + a*sqrt(2 - sqrt(2))*log(-(1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(2 - sqrt(2))/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8 - a*sqrt(2 - sqrt(2))*log((1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(2 - sqrt(2))/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8 + a*sqrt(sqrt(2) + 2)*log(-(1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8 - a*sqrt(sqrt(2) + 2)*log((1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/8 - a*sqrt(2 - sqrt(2))*atan((-2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/4 + a*sqrt(2 - sqrt(2))*atan((2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/4 - a*sqrt(sqrt(2) + 2)*atan((-2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/4 + a*sqrt(sqrt(2) + 2)*atan((2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_131():
    f = exp(acoth(a*x)/4)/x**3
    F = a**2*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(9)/8)/2 + a**2*(1 - 1/(a*x))**(sympy.S(7)/8)*(1 + 1/(a*x))**(sympy.S(1)/8)/8 + a**2*sqrt(2 - sqrt(2))*log(-(1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(2 - sqrt(2))/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/64 - a**2*sqrt(2 - sqrt(2))*log((1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(2 - sqrt(2))/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/64 + a**2*sqrt(sqrt(2) + 2)*log(-(1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/64 - a**2*sqrt(sqrt(2) + 2)*log((1 - 1/(a*x))**(sympy.S(1)/8)*sqrt(sqrt(2) + 2)/(1 + 1/(a*x))**(sympy.S(1)/8) + (1 - 1/(a*x))**(sympy.S(1)/4)/(1 + 1/(a*x))**(sympy.S(1)/4) + 1)/64 - a**2*sqrt(2 - sqrt(2))*atan((-2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/32 + a**2*sqrt(2 - sqrt(2))*atan((2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/32 - a**2*sqrt(sqrt(2) + 2)*atan((-2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/32 + a**2*sqrt(sqrt(2) + 2)*atan((2*(1 - 1/(a*x))**(sympy.S(1)/8)/(1 + 1/(a*x))**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/32
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_132():
    f = x**m*exp(4*acoth(a*x))
    F = -4*x**(m + 1)*hyper((1, m + 1), (m + 2,), a*x) + 4*x**(m + 1)/(-a*x + 1) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_133():
    f = x**m*exp(3*acoth(a*x))
    F = -3*x**(m + 1)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/2), (sympy.S.Half - m/2,), 1/(a**2*x**2))/(m + 1) + 4*x**(m + 1)*hyper((sympy.S(3)/2, -m/2 + sympy.S(-1)/2), (sympy.S.Half - m/2,), 1/(a**2*x**2))/(m + 1) - x**m*hyper((sympy.S.Half, -m/2), (1 - m/2,), 1/(a**2*x**2))/(a*m) + 4*x**m*hyper((sympy.S(3)/2, -m/2), (1 - m/2,), 1/(a**2*x**2))/(a*m)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_134():
    f = x**m*exp(2*acoth(a*x))
    F = -2*x**(m + 1)*hyper((1, m + 1), (m + 2,), a*x)/(m + 1) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_135():
    f = x**m*exp(acoth(a*x))
    F = x**(m + 1)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/2), (sympy.S.Half - m/2,), 1/(a**2*x**2))/(m + 1) + x**m*hyper((sympy.S.Half, -m/2), (1 - m/2,), 1/(a**2*x**2))/(a*m)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_136():
    f = x**m*exp(-acoth(a*x))
    F = x**(m + 1)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/2), (sympy.S.Half - m/2,), 1/(a**2*x**2))/(m + 1) - x**m*hyper((sympy.S.Half, -m/2), (1 - m/2,), 1/(a**2*x**2))/(a*m)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_137():
    f = x**m*exp(-2*acoth(a*x))
    F = -2*x**(m + 1)*hyper((1, m + 1), (m + 2,), -a*x)/(m + 1) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_138():
    f = x**m*exp(-3*acoth(a*x))
    F = -3*x**(m + 1)*hyper((sympy.S.Half, -m/2 + sympy.S(-1)/2), (sympy.S.Half - m/2,), 1/(a**2*x**2))/(m + 1) + 4*x**(m + 1)*hyper((sympy.S(3)/2, -m/2 + sympy.S(-1)/2), (sympy.S.Half - m/2,), 1/(a**2*x**2))/(m + 1) + x**m*hyper((sympy.S.Half, -m/2), (1 - m/2,), 1/(a**2*x**2))/(a*m) - 4*x**m*hyper((sympy.S(3)/2, -m/2), (1 - m/2,), 1/(a**2*x**2))/(a*m)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_139():
    f = x**m*exp(5*acoth(a*x)/2)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-5)/4, sympy.S(5)/4, -m, -1/(a*x), 1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_140():
    f = x**m*exp(3*acoth(a*x)/2)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-3)/4, sympy.S(3)/4, -m, -1/(a*x), 1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_141():
    f = x**m*exp(acoth(a*x)/2)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-1)/4, sympy.S(1)/4, -m, -1/(a*x), 1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_142():
    f = x**m*exp(-acoth(a*x)/2)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-1)/4, sympy.S(1)/4, -m, 1/(a*x), -1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_143():
    f = x**m*exp(-3*acoth(a*x)/2)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-3)/4, sympy.S(3)/4, -m, 1/(a*x), -1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_144():
    f = x**m*exp(-5*acoth(a*x)/2)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-5)/4, sympy.S(5)/4, -m, 1/(a*x), -1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_145():
    f = x**m*exp(2*acoth(x)/3)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-1)/3, sympy.S(1)/3, -m, -1/x, 1/x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_146():
    f = x**m*exp(acoth(x)/3)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-1)/6, sympy.S(1)/6, -m, -1/x, 1/x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_147():
    f = x**m*exp(acoth(a*x)/4)
    F = x**(m + 1)*appellf1(-m - 1, sympy.S(-1)/8, sympy.S(1)/8, -m, -1/(a*x), 1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_148():
    f = x**m*exp(n*acoth(a*x))
    F = x**(m + 1)*appellf1(-m - 1, -n/2, n/2, -m, -1/(a*x), 1/(a*x))/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_149():
    f = x**2*exp(n*acoth(a*x))
    F = x**3*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/3 + n*x**2*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/(6*a) + (1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 - 1)*(2*n**2 + 4)*hyper((2, 1 - n/2), (2 - n/2,), (a - 1/x)/(a + 1/x))/(3*a**3*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_150():
    f = x*exp(n*acoth(a*x))
    F = x**2*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/2 + 2*n*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (a - 1/x)/(a + 1/x))/(a**2*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_151():
    f = exp(n*acoth(a*x))
    F = 4*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (a - 1/x)/(a + 1/x))/(a*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_152():
    f = exp(n*acoth(a*x))/x
    F = 2**(n/2 + 1)*hyper((-n/2, -n/2), (1 - n/2,), (a - 1/x)/(2*a))/(n*(1 - 1/(a*x))**(n/2)) - 2*(1 + 1/(a*x))**(n/2)*hyper((1, -n/2), (1 - n/2,), (a - 1/x)/(a + 1/x))/(n*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_153():
    f = exp(n*acoth(a*x))/x**2
    F = 2**(n/2 + 1)*a*(1 - 1/(a*x))**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), (a - 1/x)/(2*a))/(2 - n)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_154():
    f = exp(n*acoth(a*x))/x**3
    F = 2**(n/2)*a**2*n*(1 - 1/(a*x))**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), (a - 1/x)/(2*a))/(2 - n) + a**2*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_155():
    f = exp(n*acoth(a*x))/x**4
    F = 2**(n/2)*a**3*(1 - 1/(a*x))**(1 - n/2)*(n**2 + 2)*hyper((-n/2, 1 - n/2), (2 - n/2,), (a - 1/x)/(2*a))/(6 - 3*n) + a**3*n*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/6 + a**2*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/(3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_156():
    f = exp(n*acoth(a*x))/x**5
    F = 2**(n/2 - 2)*a**4*n*(1 - 1/(a*x))**(1 - n/2)*(n**2 + 8)*hyper((-n/2, 1 - n/2), (2 - n/2,), (a - 1/x)/(2*a))/(6 - 3*n) + a**3*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)*(a*(n**2 + 6) + 2*n/x)/24 + a**2*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/(4*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_157():
    f = (-a*c*x + c)**p*exp(acoth(a*x))
    F = x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))*(-a*c*x + c)**p/(p + 1) + ((a - 1/x)/(a + 1/x))**(sympy.S.Half - p)*sqrt(1 + 1/(a*x))*(-a*c*x + c)**p*hyper((-p, sympy.S.Half - p), (1 - p,), 2/(x*(a + 1/x)))/(a*p*sqrt(1 - 1/(a*x))*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_158():
    f = (-a*c*x + c)**4*exp(acoth(a*x))
    F = a**4*c**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/5 - 3*a**3*c**4*x**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/4 + 17*a**2*c**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/15 - 7*a*c**4*x**2*sqrt(1 - 1/(a**2*x**2))/8 + 7*c**4*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_159():
    f = (-a*c*x + c)**3*exp(acoth(a*x))
    F = -a**3*c**3*x**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/4 + 2*a**2*c**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/3 - 5*a*c**3*x**2*sqrt(1 - 1/(a**2*x**2))/8 + 5*c**3*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_160():
    f = (-a*c*x + c)**2*exp(acoth(a*x))
    F = a**2*c**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/3 - a*c**2*x**2*sqrt(1 - 1/(a**2*x**2))/2 + c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_161():
    f = (-a*c*x + c)*exp(acoth(a*x))
    F = -a*c*x**2*sqrt(1 - 1/(a**2*x**2))/2 + c*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_162():
    f = exp(acoth(a*x))/(-a*c*x + c)
    F = -atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c) + (2*a + 2/x)/(a**2*c*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_163():
    f = exp(acoth(a*x))/(-a*c*x + c)**2
    F = -a**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(3*c**2*(a - 1/x)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_164():
    f = exp(acoth(a*x))/(-a*c*x + c)**3
    F = a**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(5*c**3*(a - 1/x)**4) - 4*a**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(15*c**3*(a - 1/x)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_165():
    f = exp(acoth(a*x))/(-a*c*x + c)**4
    F = -a**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(7*c**4*(a - 1/x)**5) + 12*a**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(35*c**4*(a - 1/x)**4) - 23*a**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(105*c**4*(a - 1/x)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_166():
    f = exp(acoth(a*x))/(-a*c*x + c)**5
    F = a**5*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(9*c**5*(a - 1/x)**6) - 8*a**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(21*c**5*(a - 1/x)**5) + 47*a**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(105*c**5*(a - 1/x)**4) - 58*a**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(315*c**5*(a - 1/x)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_167():
    f = (-a*c*x + c)**p*exp(2*acoth(a*x))
    F = 2*(-a*c*x + c)**p/(a*p) - (-a*c*x + c)**(p + 1)/(a*c*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_168():
    f = (-a*c*x + c)**5*exp(2*acoth(a*x))
    F = -c**5*(-a*x + 1)**6/(6*a) + 2*c**5*(-a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_169():
    f = (-a*c*x + c)**4*exp(2*acoth(a*x))
    F = -c**4*(-a*x + 1)**5/(5*a) + c**4*(-a*x + 1)**4/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_170():
    f = (-a*c*x + c)**3*exp(2*acoth(a*x))
    F = -c**3*(-a*x + 1)**4/(4*a) + 2*c**3*(-a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_171():
    f = (-a*c*x + c)**2*exp(2*acoth(a*x))
    F = a**2*c**2*x**3/3 - c**2*x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_172():
    f = exp(2*acoth(a*x))/(-a*c*x + c)
    F = -log(-a*x + 1)/(a*c) - 2/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_173():
    f = exp(2*acoth(a*x))/(-a*c*x + c)**2
    F = -x/(c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_174():
    f = exp(2*acoth(a*x))/(-a*c*x + c)**3
    F = 1/(2*a*c**3*(-a*x + 1)**2) - 2/(3*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_175():
    f = exp(2*acoth(a*x))/(-a*c*x + c)**4
    F = 1/(3*a*c**4*(-a*x + 1)**3) - 1/(2*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_176():
    f = (-a*c*x + c)**p*exp(3*acoth(a*x))
    F = x*(1 + 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**p/(sqrt(1 - 1/(a*x))*(p + 1)) + 3*sqrt(1 + 1/(a*x))*(-a*c*x + c)**p/(a*p*sqrt(1 - 1/(a*x))*(p + 1)) - 3*((a - 1/x)/(a + 1/x))**(sympy.S(3)/2 - p)*sqrt(1 + 1/(a*x))*(-a*c*x + c)**p*hyper((1 - p, sympy.S(3)/2 - p), (2 - p,), 2/(x*(a + 1/x)))/(a**2*p*x*(1 - p**2)*(1 - 1/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_177():
    f = (-a*c*x + c)**4*exp(3*acoth(a*x))
    F = a**4*c**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/5 - a**3*c**4*x**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/4 + 3*a*c**4*x**2*sqrt(1 - 1/(a**2*x**2))/8 - 3*c**4*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_178():
    f = (-a*c*x + c)**3*exp(3*acoth(a*x))
    F = -a**3*c**3*x**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/4 + 3*a*c**3*x**2*sqrt(1 - 1/(a**2*x**2))/8 - 3*c**3*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_179():
    f = (-a*c*x + c)**2*exp(3*acoth(a*x))
    F = a**2*c**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/3 + a*c**2*x**2*sqrt(1 - 1/(a**2*x**2))/2 - c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_180():
    f = (-a*c*x + c)*exp(3*acoth(a*x))
    F = -a*c*x**2*sqrt(1 - 1/(a**2*x**2))/2 - 2*c*x*sqrt(1 - 1/(a**2*x**2)) - 3*c*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_181():
    f = exp(3*acoth(a*x))/(-a*c*x + c)
    F = -atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c) + (8*a + 8/x)/(3*a**2*c*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) + 4/(3*a**2*c*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_182():
    f = exp(3*acoth(a*x))/(-a*c*x + c)**2
    F = -a**4*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(5*c**2*(a - 1/x)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_183():
    f = exp(3*acoth(a*x))/(-a*c*x + c)**3
    F = a**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(7*c**3*(a - 1/x)**6) - 6*a**4*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(35*c**3*(a - 1/x)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_184():
    f = exp(3*acoth(a*x))/(-a*c*x + c)**4
    F = -47*(a + 1/x)**5/(315*a**6*c**4*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) + 16*(a + 1/x)**6/(63*a**7*c**4*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - (a + 1/x)**7/(9*a**8*c**4*(1 - 1/(a**2*x**2))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_185():
    f = exp(3*acoth(a*x))/(-a*c*x + c)**5
    F = -152*(a + 1/x)**5/(1155*a**6*c**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) + 79*(a + 1/x)**6/(231*a**7*c**5*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 10*(a + 1/x)**7/(33*a**8*c**5*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) + (a + 1/x)**8/(11*a**9*c**5*(1 - 1/(a**2*x**2))**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_186():
    f = (-a*c*x + c)**p*exp(4*acoth(a*x))
    F = 4*c*(-a*c*x + c)**(p - 1)/(a*(1 - p)) + 4*(-a*c*x + c)**p/(a*p) - (-a*c*x + c)**(p + 1)/(a*c*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_187():
    f = (-a*c*x + c)**5*exp(4*acoth(a*x))
    F = -c**5*(-a*x + 1)**6/(6*a) + 4*c**5*(-a*x + 1)**5/(5*a) - c**5*(-a*x + 1)**4/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_188():
    f = (-a*c*x + c)**4*exp(4*acoth(a*x))
    F = a**4*c**4*x**5/5 - 2*a**2*c**4*x**3/3 + c**4*x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_189():
    f = (-a*c*x + c)**3*exp(4*acoth(a*x))
    F = -c**3*(a*x + 1)**4/(4*a) + 2*c**3*(a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_190():
    f = (-a*c*x + c)**2*exp(4*acoth(a*x))
    F = c**2*(a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_191():
    f = (-a*c*x + c)*exp(4*acoth(a*x))
    F = -a*c*x**2/2 - 3*c*x - 4*c*log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_192():
    f = exp(4*acoth(a*x))/(-a*c*x + c)
    F = -log(-a*x + 1)/(a*c) - 4/(a*c*(-a*x + 1)) + 2/(a*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_193():
    f = exp(4*acoth(a*x))/(-a*c*x + c)**2
    F = (a*x + 1)**3/(6*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_194():
    f = exp(4*acoth(a*x))/(-a*c*x + c)**3
    F = 1/(2*a*c**3*(-a*x + 1)**2) - 4/(3*a*c**3*(-a*x + 1)**3) + 1/(a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_195():
    f = exp(4*acoth(a*x))/(-a*c*x + c)**4
    F = 1/(3*a*c**4*(-a*x + 1)**3) - 1/(a*c**4*(-a*x + 1)**4) + 4/(5*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_196():
    f = (-a*c*x + c)**p*exp(-acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(-p + sympy.S(-1)/2)*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))*(-a*c*x + c)**p*hyper((-p - 1, -p + sympy.S(-1)/2), (-p,), 2/(x*(a + 1/x)))/(p + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_197():
    f = (-a*c*x + c)**3*exp(-acoth(a*x))
    F = -a**3*c**3*x**4*sqrt(1 - 1/(a**2*x**2))/4 + 4*a**2*c**3*x**3*sqrt(1 - 1/(a**2*x**2))/3 - 27*a*c**3*x**2*sqrt(1 - 1/(a**2*x**2))/8 + 20*c**3*x*sqrt(1 - 1/(a**2*x**2))/3 - 35*c**3*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_198():
    f = (-a*c*x + c)**2*exp(-acoth(a*x))
    F = a**2*c**2*x**3*sqrt(1 - 1/(a**2*x**2))/3 - 3*a*c**2*x**2*sqrt(1 - 1/(a**2*x**2))/2 + 11*c**2*x*sqrt(1 - 1/(a**2*x**2))/3 - 5*c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_199():
    f = (-a*c*x + c)*exp(-acoth(a*x))
    F = -a*c*x**2*sqrt(1 - 1/(a**2*x**2))/2 + 2*c*x*sqrt(1 - 1/(a**2*x**2)) - 3*c*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_200():
    f = exp(-acoth(a*x))/(-a*c*x + c)
    F = -atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_201():
    f = exp(-acoth(a*x))/(-a*c*x + c)**2
    F = -sqrt(1 - 1/(a**2*x**2))/(c**2*(a - 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_202():
    f = exp(-acoth(a*x))/(-a*c*x + c)**3
    F = a*sqrt(1 - 1/(a**2*x**2))/(3*c**3*(a - 1/x)**2) - 2*sqrt(1 - 1/(a**2*x**2))/(3*c**3*(a - 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_203():
    f = exp(-acoth(a*x))/(-a*c*x + c)**4
    F = -a**2*sqrt(1 - 1/(a**2*x**2))/(5*c**4*(a - 1/x)**3) + 8*a*sqrt(1 - 1/(a**2*x**2))/(15*c**4*(a - 1/x)**2) - 7*sqrt(1 - 1/(a**2*x**2))/(15*c**4*(a - 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_204():
    f = exp(-acoth(a*x))/(-a*c*x + c)**5
    F = a**3*sqrt(1 - 1/(a**2*x**2))/(7*c**5*(a - 1/x)**4) - 18*a**2*sqrt(1 - 1/(a**2*x**2))/(35*c**5*(a - 1/x)**3) + 23*a*sqrt(1 - 1/(a**2*x**2))/(35*c**5*(a - 1/x)**2) - 12*sqrt(1 - 1/(a**2*x**2))/(35*c**5*(a - 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_205():
    f = (-a*c*x + c)**p*exp(-2*acoth(a*x))
    F = (-a*c*x + c)**(p + 2)*hyper((1, p + 2), (p + 3,), -a*x/2 + sympy.S.Half)/(2*a*c**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_206():
    f = (-a*c*x + c)**4*exp(-2*acoth(a*x))
    F = 16*c**4*x - c**4*(-a*x + 1)**5/(5*a) - c**4*(-a*x + 1)**4/(2*a) - 4*c**4*(-a*x + 1)**3/(3*a) - 4*c**4*(-a*x + 1)**2/a - 32*c**4*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_207():
    f = (-a*c*x + c)**3*exp(-2*acoth(a*x))
    F = 8*c**3*x - c**3*(-a*x + 1)**4/(4*a) - 2*c**3*(-a*x + 1)**3/(3*a) - 2*c**3*(-a*x + 1)**2/a - 16*c**3*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_208():
    f = (-a*c*x + c)**2*exp(-2*acoth(a*x))
    F = 4*c**2*x - c**2*(-a*x + 1)**3/(3*a) - c**2*(-a*x + 1)**2/a - 8*c**2*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_209():
    f = (-a*c*x + c)*exp(-2*acoth(a*x))
    F = -a*c*x**2/2 + 3*c*x - 4*c*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_210():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)
    F = -log(a*x + 1)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_211():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**2
    F = -atanh(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_212():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**3
    F = -atanh(a*x)/(2*a*c**3) - 1/(2*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_213():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**4
    F = -atanh(a*x)/(4*a*c**4) - 1/(4*a*c**4*(-a*x + 1)) - 1/(4*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_214():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**5
    F = -atanh(a*x)/(8*a*c**5) - 1/(8*a*c**5*(-a*x + 1)) - 1/(8*a*c**5*(-a*x + 1)**2) - 1/(6*a*c**5*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_215():
    f = (-a*c*x + c)**p*exp(-3*acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(-p + sympy.S(-3)/2)*(1 - 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**p*hyper((-p + sympy.S(-3)/2, -p - 1), (-p,), 2/(x*(a + 1/x)))/(sqrt(1 + 1/(a*x))*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_216():
    f = (-a*c*x + c)**3*exp(-3*acoth(a*x))
    F = -a**3*c**3*x**4*sqrt(1 - 1/(a**2*x**2))/4 + 2*a**2*c**3*x**3*sqrt(1 - 1/(a**2*x**2)) - 67*a*c**3*x**2*sqrt(1 - 1/(a**2*x**2))/8 + 30*c**3*x*sqrt(1 - 1/(a**2*x**2)) - 315*c**3*atanh(sqrt(1 - 1/(a**2*x**2)))/(8*a) + 32*c**3*(a - 1/x)/(a**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_217():
    f = (-a*c*x + c)**2*exp(-3*acoth(a*x))
    F = a**2*c**2*x**3*sqrt(1 - 1/(a**2*x**2))/3 - 5*a*c**2*x**2*sqrt(1 - 1/(a**2*x**2))/2 + 35*c**2*x*sqrt(1 - 1/(a**2*x**2))/3 - 35*c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a) + 16*c**2*(a - 1/x)/(a**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_218():
    f = (-a*c*x + c)*exp(-3*acoth(a*x))
    F = -a*c*x**2*sqrt(1 - 1/(a**2*x**2))/2 + 4*c*x*sqrt(1 - 1/(a**2*x**2)) - 15*c*atanh(sqrt(1 - 1/(a**2*x**2)))/(2*a) + 8*c*(a - 1/x)/(a**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_219():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)
    F = -atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c) + (2*a - 2/x)/(a**2*c*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_220():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**2
    F = (a - 1/x)/(a**2*c**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_221():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**3
    F = 1/(a*c**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_222():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**4
    F = 2/(3*a*c**4*sqrt(1 - 1/(a**2*x**2))) - 1/(3*a**2*c**4*x**2*sqrt(1 - 1/(a**2*x**2))*(a - 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_223():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**5
    F = (5*a + 2/x)/(5*a**2*c**5*sqrt(1 - 1/(a**2*x**2))) - (4*a + 4/x)/(5*a**2*c**5*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) + (a + 1/x)**2/(5*a**3*c**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_224():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**6
    F = (35*a + 13/x)/(35*a**2*c**6*sqrt(1 - 1/(a**2*x**2))) - (46*a + 46/x)/(35*a**2*c**6*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) + 24*(a + 1/x)**2/(35*a**3*c**6*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) - (a + 1/x)**3/(7*a**4*c**6*(1 - 1/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_225():
    f = (-a*c*x + c)**(sympy.S(9)/2)*exp(acoth(a*x))
    F = 128*(1 + 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(9)/2)/(231*a**2*x*(1 - 1/(a*x))**(sympy.S(9)/2)) - 768*(1 + 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(9)/2)/(385*a**3*x**2*(1 - 1/(a*x))**(sympy.S(9)/2)) + 2*x*(1 + 1/(a*x))**(sympy.S(3)/2)*(a - 1/x)**4*(-a*c*x + c)**(sympy.S(9)/2)/(11*a**4*(1 - 1/(a*x))**(sympy.S(9)/2)) - 32*(1 + 1/(a*x))**(sympy.S(3)/2)*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(9)/2)/(99*a**4*(1 - 1/(a*x))**(sympy.S(9)/2)) + 9088*(1 + 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(9)/2)/(3465*a**4*x**3*(1 - 1/(a*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_226():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(acoth(a*x))
    F = -8*(1 + 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(7)/2)/(21*a*(1 - 1/(a*x))**(sympy.S(7)/2)) + 48*(1 + 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(7)/2)/(35*a**2*x*(1 - 1/(a*x))**(sympy.S(7)/2)) + 2*x*(1 + 1/(a*x))**(sympy.S(3)/2)*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(7)/2)/(9*a**3*(1 - 1/(a*x))**(sympy.S(7)/2)) - 568*(1 + 1/(a*x))**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(7)/2)/(315*a**3*x**2*(1 - 1/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_227():
    f = sqrt(-a*c*x + c)*exp(acoth(a*x))
    F = 2*(a*x + 1)*sqrt(-a*c*x + c)*exp(acoth(a*x))/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_228():
    f = exp(acoth(a*x))/sqrt(-a*c*x + c)
    F = 2*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/sqrt(-a*c*x + c) - 2*sqrt(2)*sqrt(1 - 1/(a*x))*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(sqrt(a)*sqrt(-a*c*x + c)*sqrt(1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_229():
    f = exp(acoth(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -sqrt(2)*sqrt(a)*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(2*(-a*c*x + c)**(sympy.S(3)/2)*(1/x)**(sympy.S(3)/2)) - a*x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_230():
    f = exp(acoth(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = sqrt(2)*a**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(16*(-a*c*x + c)**(sympy.S(5)/2)*(1/x)**(sympy.S(5)/2)) - a**3*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(4*(a - 1/x)**2*(-a*c*x + c)**(sympy.S(5)/2)) + a**2*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((8*a - 8/x)*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_231():
    f = exp(acoth(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = -sqrt(2)*a**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(64*(-a*c*x + c)**(sympy.S(7)/2)*(1/x)**(sympy.S(7)/2)) + a**4*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(16*(a - 1/x)**2*(-a*c*x + c)**(sympy.S(7)/2)) - a**4*x**2*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(6*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(7)/2)) - a**3*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))/((32*a - 32/x)*(-a*c*x + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_232():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(2*acoth(a*x))
    F = 4*(-a*c*x + c)**(sympy.S(7)/2)/(7*a) - 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_233():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(2*acoth(a*x))
    F = 4*(-a*c*x + c)**(sympy.S(5)/2)/(5*a) - 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_234():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(2*acoth(a*x))
    F = 4*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) - 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_235():
    f = sqrt(-a*c*x + c)*exp(2*acoth(a*x))
    F = 4*sqrt(-a*c*x + c)/a - 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_236():
    f = exp(2*acoth(a*x))/sqrt(-a*c*x + c)
    F = -4/(a*sqrt(-a*c*x + c)) - 2*sqrt(-a*c*x + c)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_237():
    f = exp(2*acoth(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -4/(3*a*(-a*c*x + c)**(sympy.S(3)/2)) + 2/(a*c*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_238():
    f = exp(2*acoth(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = -4/(5*a*(-a*c*x + c)**(sympy.S(5)/2)) + 2/(3*a*c*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_239():
    f = exp(2*acoth(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = -4/(7*a*(-a*c*x + c)**(sympy.S(7)/2)) + 2/(5*a*c*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_240():
    f = (-a*c*x + c)**(sympy.S(9)/2)*exp(3*acoth(a*x))
    F = -8*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(9)/2)/(33*a*(1 - 1/(a*x))**(sympy.S(9)/2)) + 16*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(9)/2)/(21*a**2*x*(1 - 1/(a*x))**(sympy.S(9)/2)) + 2*x*(1 + 1/(a*x))**(sympy.S(5)/2)*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(9)/2)/(11*a**3*(1 - 1/(a*x))**(sympy.S(9)/2)) - 856*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(9)/2)/(1155*a**3*x**2*(1 - 1/(a*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_241():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(3*acoth(a*x))
    F = 2*x*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(7)/2)/(9*(1 - 1/(a*x))**(sympy.S(7)/2)) - 44*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(7)/2)/(63*a*(1 - 1/(a*x))**(sympy.S(7)/2)) + 214*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(7)/2)/(315*a**2*x*(1 - 1/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_242():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(3*acoth(a*x))
    F = 2*x*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(5)/2)/(7*(1 - 1/(a*x))**(sympy.S(5)/2)) - 18*(1 + 1/(a*x))**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(5)/2)/(35*a*(1 - 1/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_243():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(3*acoth(a*x))
    F = 2*(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)*exp(3*acoth(a*x))/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_244():
    f = sqrt(-a*c*x + c)*exp(3*acoth(a*x))
    F = 2*x*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(3*sqrt(1 - 1/(a*x))) + 4*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(a*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(a**(sympy.S(3)/2)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_245():
    f = exp(3*acoth(a*x))/sqrt(-a*c*x + c)
    F = 2*a*x*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/((a - 1/x)*sqrt(-a*c*x + c)) - 6*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/((a - 1/x)*sqrt(-a*c*x + c)) - 3*sqrt(2)*sqrt(1 - 1/(a*x))*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(sqrt(a)*sqrt(-a*c*x + c)*sqrt(1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_246():
    f = exp(3*acoth(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -3*sqrt(2)*sqrt(a)*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(8*(-a*c*x + c)**(sympy.S(3)/2)*(1/x)**(sympy.S(3)/2)) - a**2*x*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(2*(a - 1/x)**2*(-a*c*x + c)**(sympy.S(3)/2)) - 3*a*x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/((4*a - 4/x)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_247():
    f = exp(3*acoth(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = sqrt(2)*a**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(32*(-a*c*x + c)**(sympy.S(5)/2)*(1/x)**(sympy.S(5)/2)) - a**4*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(6*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(5)/2)) + a**3*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(24*(a - 1/x)**2*(-a*c*x + c)**(sympy.S(5)/2)) + a**2*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((16*a - 16/x)*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_248():
    f = exp(3*acoth(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = -3*sqrt(2)*a**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(512*(-a*c*x + c)**(sympy.S(7)/2)*(1/x)**(sympy.S(7)/2)) + a**5*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(32*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(7)/2)) - a**5*x**2*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(8*(a - 1/x)**4*(-a*c*x + c)**(sympy.S(7)/2)) - a**4*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(128*(a - 1/x)**2*(-a*c*x + c)**(sympy.S(7)/2)) - 3*a**3*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))/((256*a - 256/x)*(-a*c*x + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_249():
    f = exp(-acoth(a*x))/sqrt(-a*c*x + c)
    F = (2*a*x + 2)*exp(-acoth(a*x))/(a*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_250():
    f = exp(-acoth(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -sqrt(2)*sqrt(a)*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/((-a*c*x + c)**(sympy.S(3)/2)*(1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_251():
    f = exp(-acoth(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = sqrt(2)*a**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(4*(-a*c*x + c)**(sympy.S(5)/2)*(1/x)**(sympy.S(5)/2)) - a**2*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((2*a - 2/x)*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_252():
    f = exp(-acoth(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = -3*sqrt(2)*a**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(32*(-a*c*x + c)**(sympy.S(7)/2)*(1/x)**(sympy.S(7)/2)) + 3*a**3*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))/((16*a - 16/x)*(-a*c*x + c)**(sympy.S(7)/2)) - a**3*x**2*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))/(4*(a - 1/x)**2*(-a*c*x + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_253():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(-2*acoth(a*x))
    F = 32*sqrt(2)*c**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a - 32*c**3*sqrt(-a*c*x + c)/a - 16*c**2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) - 8*c*(-a*c*x + c)**(sympy.S(5)/2)/(5*a) - 4*(-a*c*x + c)**(sympy.S(7)/2)/(7*a) - 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_254():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(-2*acoth(a*x))
    F = 16*sqrt(2)*c**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a - 16*c**2*sqrt(-a*c*x + c)/a - 8*c*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) - 4*(-a*c*x + c)**(sympy.S(5)/2)/(5*a) - 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_255():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(-2*acoth(a*x))
    F = 8*sqrt(2)*c**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a - 8*c*sqrt(-a*c*x + c)/a - 4*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) - 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_256():
    f = sqrt(-a*c*x + c)*exp(-2*acoth(a*x))
    F = 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a - 4*sqrt(-a*c*x + c)/a - 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_257():
    f = exp(-2*acoth(a*x))/sqrt(-a*c*x + c)
    F = -2*sqrt(-a*c*x + c)/(a*c) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_258():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_259():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = -1/(a*c**2*sqrt(-a*c*x + c)) + sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(2*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_260():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = -1/(3*a*c**2*(-a*c*x + c)**(sympy.S(3)/2)) - 1/(2*a*c**3*sqrt(-a*c*x + c)) + sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(4*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_261():
    f = exp(-2*acoth(a*x))/(-a*c*x + c)**(sympy.S(9)/2)
    F = -1/(5*a*c**2*(-a*c*x + c)**(sympy.S(5)/2)) - 1/(6*a*c**3*(-a*c*x + c)**(sympy.S(3)/2)) - 1/(4*a*c**4*sqrt(-a*c*x + c)) + sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(8*a*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_262():
    f = (-a*c*x + c)**(sympy.S(9)/2)*exp(-3*acoth(a*x))
    F = 4096*(-a*c*x + c)**(sympy.S(9)/2)/(231*a**4*x**3*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x))) - 40960*(-a*c*x + c)**(sympy.S(9)/2)/(231*a**5*x**4*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x))) + 2*x*(a - 1/x)**6*(-a*c*x + c)**(sympy.S(9)/2)/(11*a**6*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x))) - 16*(a - 1/x)**5*(-a*c*x + c)**(sympy.S(9)/2)/(33*a**6*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x))) + 320*(a - 1/x)**4*(-a*c*x + c)**(sympy.S(9)/2)/(231*a**6*x*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x))) - 1024*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(9)/2)/(231*a**6*x**2*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x))) - 94208*(-a*c*x + c)**(sympy.S(9)/2)/(231*a**6*x**5*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_263():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(-3*acoth(a*x))
    F = -512*(-a*c*x + c)**(sympy.S(7)/2)/(63*a**3*x**2*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))) + 5120*(-a*c*x + c)**(sympy.S(7)/2)/(63*a**4*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))) + 2*x*(a - 1/x)**5*(-a*c*x + c)**(sympy.S(7)/2)/(9*a**5*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))) - 40*(a - 1/x)**4*(-a*c*x + c)**(sympy.S(7)/2)/(63*a**5*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))) + 128*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(7)/2)/(63*a**5*x*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))) + 11776*(-a*c*x + c)**(sympy.S(7)/2)/(63*a**5*x**4*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_264():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(-3*acoth(a*x))
    F = 128*(-a*c*x + c)**(sympy.S(5)/2)/(35*a**2*x*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))) - 256*(-a*c*x + c)**(sympy.S(5)/2)/(7*a**3*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))) + 2*x*(a - 1/x)**4*(-a*c*x + c)**(sympy.S(5)/2)/(7*a**4*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))) - 32*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(5)/2)/(35*a**4*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))) - 2944*(-a*c*x + c)**(sympy.S(5)/2)/(35*a**4*x**3*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_265():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(-3*acoth(a*x))
    F = -8*(-a*c*x + c)**(sympy.S(3)/2)/(5*a*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))) + 16*(-a*c*x + c)**(sympy.S(3)/2)/(a**2*x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))) + 2*x*(a - 1/x)**3*(-a*c*x + c)**(sympy.S(3)/2)/(5*a**3*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))) + 184*(-a*c*x + c)**(sympy.S(3)/2)/(5*a**3*x**2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_266():
    f = sqrt(-a*c*x + c)*exp(-3*acoth(a*x))
    F = 2*x*sqrt(-a*c*x + c)/(3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 20*sqrt(-a*c*x + c)/(3*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 46*sqrt(-a*c*x + c)/(3*a**2*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_267():
    f = exp(-3*acoth(a*x))/sqrt(-a*c*x + c)
    F = 2*x*sqrt(1 - 1/(a*x))/(sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)) + 6*sqrt(1 - 1/(a*x))/(a*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_268():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -(2*a*x + 2)*exp(-3*acoth(a*x))/(a*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_269():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = -sqrt(2)*a**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(2*(-a*c*x + c)**(sympy.S(5)/2)*(1/x)**(sympy.S(5)/2)) + a*x**2*(1 - 1/(a*x))**(sympy.S(5)/2)/(sqrt(1 + 1/(a*x))*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_270():
    f = exp(-3*acoth(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = 3*sqrt(2)*a**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(8*(-a*c*x + c)**(sympy.S(7)/2)*(1/x)**(sympy.S(7)/2)) - 3*a**2*x**3*(1 - 1/(a*x))**(sympy.S(7)/2)/(4*sqrt(1 + 1/(a*x))*(-a*c*x + c)**(sympy.S(7)/2)) - a**2*x**2*(1 - 1/(a*x))**(sympy.S(7)/2)/(sqrt(1 + 1/(a*x))*(2*a - 2/x)*(-a*c*x + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_271():
    f = x*(x + 1)*exp(acoth(x))
    F = x**3*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(5)/2)/3 + x**2*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(3)/2)/3 + x*sqrt((x - 1)/x)*sqrt(1 + 1/x) + atanh(sqrt((x - 1)/x)*sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_272():
    f = (x + 1)*exp(acoth(x))
    F = x**2*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(3)/2)/2 + 3*x*sqrt((x - 1)/x)*sqrt(1 + 1/x)/2 + 3*atanh(sqrt((x - 1)/x)*sqrt(1 + 1/x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_273():
    f = x*(1 - x)*exp(acoth(x))
    F = -x**3*(1 - 1/x**2)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_274():
    f = (1 - x)*exp(acoth(x))
    F = -x**2*sqrt(1 - 1/x**2)/2 + atanh(sqrt(1 - 1/x**2))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_275():
    f = x*(x + 1)**2*exp(acoth(x))
    F = x**4*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(7)/2)/4 + x**3*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(5)/2)/4 + 5*x**2*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(3)/2)/8 + 15*x*sqrt((x - 1)/x)*sqrt(1 + 1/x)/8 + 15*atanh(sqrt((x - 1)/x)*sqrt(1 + 1/x))/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_276():
    f = (x + 1)**2*exp(acoth(x))
    F = x**3*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(5)/2)/3 + 5*x**2*sqrt((x - 1)/x)*(1 + 1/x)**(sympy.S(3)/2)/6 + 5*x*sqrt((x - 1)/x)*sqrt(1 + 1/x)/2 + 5*atanh(sqrt((x - 1)/x)*sqrt(1 + 1/x))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_277():
    f = x*(1 - x)**2*exp(acoth(x))
    F = x**4*(1 - 1/x**2)**(sympy.S(3)/2)/4 - x**3*(1 - 1/x**2)**(sympy.S(3)/2)/3 + x**2*sqrt(1 - 1/x**2)/8 - atanh(sqrt(1 - 1/x**2))/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_278():
    f = (1 - x)**2*exp(acoth(x))
    F = x**3*(1 - 1/x**2)**(sympy.S(3)/2)/3 - x**2*sqrt(1 - 1/x**2)/2 + atanh(sqrt(1 - 1/x**2))/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_279():
    f = x*exp(acoth(x))/(x + 1)
    F = x*sqrt((x - 1)/x)*sqrt(1 + 1/x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_280():
    f = exp(acoth(x))/(x + 1)
    F = atanh(sqrt((x - 1)/x)*sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_281():
    f = x*exp(acoth(x))/(1 - x)
    F = -x*sqrt(1 - 1/x**2) - 2*atanh(sqrt(1 - 1/x**2)) + (2 + 2/x)/sqrt(1 - 1/x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_282():
    f = exp(acoth(x))/(1 - x)
    F = -atanh(sqrt(1 - 1/x**2)) + (2 + 2/x)/sqrt(1 - 1/x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_283():
    f = x*exp(acoth(x))/(x + 1)**2
    F = -sqrt((x - 1)/x)/sqrt(1 + 1/x) + atanh(sqrt((x - 1)/x)*sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_284():
    f = exp(acoth(x))/(x + 1)**2
    F = sqrt((x - 1)/x)/sqrt(1 + 1/x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_285():
    f = x*exp(acoth(x))/(1 - x)**2
    F = atanh(sqrt(1 - 1/x**2)) - (3 + 5/x)/(3*sqrt(1 - 1/x**2)) - (4 + 4/x)/(3*(1 - 1/x**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_286():
    f = exp(acoth(x))/(1 - x)**2
    F = -(1 - 1/x**2)**(sympy.S(3)/2)/(3*(1 - 1/x)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_287():
    f = x**m*sqrt(-a*c*x + c)*exp(acoth(a*x))
    F = 2*x**(m + 1)*sqrt(-a*c*x + c)*hyper((sympy.S(-1)/2, -m + sympy.S(-3)/2), (-m + sympy.S(-1)/2,), -1/(a*x))/(sqrt(1 - 1/(a*x))*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_288():
    f = x**2*sqrt(-a*c*x + c)*exp(acoth(a*x))
    F = 2*x**3*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(7*sqrt(1 - 1/(a*x))) - 8*x**2*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(35*a*sqrt(1 - 1/(a*x))) + 16*x*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(105*a**2*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_289():
    f = x*sqrt(-a*c*x + c)*exp(acoth(a*x))
    F = 2*x**2*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(5*sqrt(1 - 1/(a*x))) - 4*x*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(15*a*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_290():
    f = sqrt(-a*c*x + c)*exp(acoth(a*x))
    F = 2*(a*x + 1)*sqrt(-a*c*x + c)*exp(acoth(a*x))/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_291():
    f = sqrt(-a*c*x + c)*exp(acoth(a*x))/x
    F = 2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/sqrt(1 - 1/(a*x)) - 2*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(sqrt(a)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_292():
    f = sqrt(-a*c*x + c)*exp(acoth(a*x))/x**2
    F = -sqrt(a)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/sqrt(1 - 1/(a*x)) - sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(x*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_293():
    f = x**3*sqrt(-a*c*x + c)*exp(2*acoth(a*x))
    F = 4*sqrt(-a*c*x + c)/a**4 - 14*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**4*c) + 18*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**4*c**2) - 10*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**4*c**3) + 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a**4*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_294():
    f = x**2*sqrt(-a*c*x + c)*exp(2*acoth(a*x))
    F = 4*sqrt(-a*c*x + c)/a**3 - 10*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**3*c) + 8*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**3*c**2) - 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**3*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_295():
    f = x*sqrt(-a*c*x + c)*exp(2*acoth(a*x))
    F = 4*sqrt(-a*c*x + c)/a**2 - 2*(-a*c*x + c)**(sympy.S(3)/2)/(a**2*c) + 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_296():
    f = sqrt(-a*c*x + c)*exp(2*acoth(a*x))
    F = 4*sqrt(-a*c*x + c)/a - 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_297():
    f = sqrt(-a*c*x + c)*exp(2*acoth(a*x))/x
    F = 2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) + 2*sqrt(-a*c*x + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_298():
    f = sqrt(-a*c*x + c)*exp(2*acoth(a*x))/x**2
    F = 3*a*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) + sqrt(-a*c*x + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_299():
    f = sqrt(-a*c*x + c)*exp(2*acoth(a*x))/x**3
    F = 7*a**2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/4 + 7*a*sqrt(-a*c*x + c)/(4*x) + sqrt(-a*c*x + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_300():
    f = sqrt(-a*c*x + c)*exp(2*acoth(a*x))/x**4
    F = 11*a**3*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/8 + 11*a**2*sqrt(-a*c*x + c)/(8*x) + 11*a*sqrt(-a*c*x + c)/(12*x**2) + sqrt(-a*c*x + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_301():
    f = sqrt(-a*c*x + c)*exp(2*acoth(a*x))/x**5
    F = 75*a**4*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/64 + 75*a**3*sqrt(-a*c*x + c)/(64*x) + 25*a**2*sqrt(-a*c*x + c)/(32*x**2) + 5*a*sqrt(-a*c*x + c)/(8*x**3) + sqrt(-a*c*x + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_302():
    f = x**3*sqrt(-a*c*x + c)*exp(3*acoth(a*x))
    F = 2*x**4*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(9*sqrt(1 - 1/(a*x))) + 38*x**3*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(63*a*sqrt(1 - 1/(a*x))) + 92*x**2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(105*a**2*sqrt(1 - 1/(a*x))) + 472*x*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(315*a**3*sqrt(1 - 1/(a*x))) + 1576*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(315*a**4*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(a**(sympy.S(9)/2)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_303():
    f = x**2*sqrt(-a*c*x + c)*exp(3*acoth(a*x))
    F = 2*x**3*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(7*sqrt(1 - 1/(a*x))) + 6*x**2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(7*a*sqrt(1 - 1/(a*x))) + 32*x*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(21*a**2*sqrt(1 - 1/(a*x))) + 104*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(21*a**3*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(a**(sympy.S(7)/2)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_304():
    f = x*sqrt(-a*c*x + c)*exp(3*acoth(a*x))
    F = 2*x**2*(1 + 1/(a*x))**(sympy.S(5)/2)*sqrt(-a*c*x + c)/(5*sqrt(1 - 1/(a*x))) + 2*x*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(3*a*sqrt(1 - 1/(a*x))) + 4*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(a**2*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(a**(sympy.S(5)/2)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_305():
    f = sqrt(-a*c*x + c)*exp(3*acoth(a*x))
    F = 2*x*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(3*sqrt(1 - 1/(a*x))) + 4*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(a*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(a**(sympy.S(3)/2)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_306():
    f = sqrt(-a*c*x + c)*exp(3*acoth(a*x))/x
    F = 2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/sqrt(1 - 1/(a*x)) + 2*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(sqrt(a)*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/(sqrt(a)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_307():
    f = sqrt(-a*c*x + c)*exp(3*acoth(a*x))/x**2
    F = 5*sqrt(a)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/sqrt(1 - 1/(a*x)) - 4*sqrt(2)*sqrt(a)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/sqrt(1 - 1/(a*x)) + sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(x*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_308():
    f = sqrt(-a*c*x + c)*exp(3*acoth(a*x))/x**3
    F = 23*a**(sympy.S(3)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(4*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*a**(sympy.S(3)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/sqrt(1 - 1/(a*x)) + a*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(2*x*sqrt(1 - 1/(a*x))) + 7*a*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(4*x*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_309():
    f = sqrt(-a*c*x + c)*exp(3*acoth(a*x))/x**4
    F = 45*a**(sympy.S(5)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(8*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*a**(sympy.S(5)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/sqrt(1 - 1/(a*x)) + 3*a**2*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(4*x*sqrt(1 - 1/(a*x))) + 13*a**2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(8*x*sqrt(1 - 1/(a*x))) + a*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(3*x**2*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_310():
    f = sqrt(-a*c*x + c)*exp(3*acoth(a*x))/x**5
    F = 363*a**(sympy.S(7)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(64*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*a**(sympy.S(7)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*atanh(sqrt(2)*sqrt(1/x)/(sqrt(a)*sqrt(1 + 1/(a*x))))/sqrt(1 - 1/(a*x)) + 21*a**3*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(32*x*sqrt(1 - 1/(a*x))) + 107*a**3*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(64*x*sqrt(1 - 1/(a*x))) + 11*a**2*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(24*x**2*sqrt(1 - 1/(a*x))) + a*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(4*x**3*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_311():
    f = x*(x + 1)**(sympy.S(3)/2)*exp(acoth(x))
    F = 2*x**2*sqrt(-(1 - x)/x)*(x + 1)**(sympy.S(3)/2)/(7*(1 + 1/x)**(sympy.S(3)/2)) + 8*x*sqrt(-(1 - x)/x)*(x + 1)**(sympy.S(3)/2)/(7*(1 + 1/x)**(sympy.S(3)/2)) + 46*sqrt(-(1 - x)/x)*(x + 1)**(sympy.S(3)/2)/(21*(1 + 1/x)**(sympy.S(3)/2)) + 92*sqrt(-(1 - x)/x)*(x + 1)**(sympy.S(3)/2)/(21*x*(1 + 1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_312():
    f = (x + 1)**(sympy.S(3)/2)*exp(acoth(x))
    F = 2*x*sqrt(-(1 - x)/x)*(x + 1)**(sympy.S(3)/2)/(5*(1 + 1/x)**(sympy.S(3)/2)) + 28*sqrt(-(1 - x)/x)*(x + 1)**(sympy.S(3)/2)/(15*(1 + 1/x)**(sympy.S(3)/2)) + 86*sqrt(-(1 - x)/x)*(x + 1)**(sympy.S(3)/2)/(15*x*(1 + 1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_313():
    f = x*(1 - x)**(sympy.S(3)/2)*exp(acoth(x))
    F = 2*x**2*(1 + 1/x)**(sympy.S(3)/2)*(1 - x)**(sympy.S(3)/2)/(7*(1 - 1/x)**(sympy.S(3)/2)) - 22*x*(1 + 1/x)**(sympy.S(3)/2)*(1 - x)**(sympy.S(3)/2)/(35*(1 - 1/x)**(sympy.S(3)/2)) + 44*(1 + 1/x)**(sympy.S(3)/2)*(1 - x)**(sympy.S(3)/2)/(105*(1 - 1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_314():
    f = (1 - x)**(sympy.S(3)/2)*exp(acoth(x))
    F = 2*x*(1 + 1/x)**(sympy.S(3)/2)*(1 - x)**(sympy.S(3)/2)/(5*(1 - 1/x)**(sympy.S(3)/2)) - 14*(1 + 1/x)**(sympy.S(3)/2)*(1 - x)**(sympy.S(3)/2)/(15*(1 - 1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_315():
    f = x*sqrt(x + 1)*exp(acoth(x))
    F = 2*x**2*sqrt(-(1 - x)/x)*sqrt(x + 1)/(5*sqrt(1 + 1/x)) + 6*x*sqrt(-(1 - x)/x)*sqrt(x + 1)/(5*sqrt(1 + 1/x)) + 12*sqrt(-(1 - x)/x)*sqrt(x + 1)/(5*sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_316():
    f = sqrt(x + 1)*exp(acoth(x))
    F = 2*x*sqrt(-(1 - x)/x)*sqrt(x + 1)/(3*sqrt(1 + 1/x)) + 10*sqrt(-(1 - x)/x)*sqrt(x + 1)/(3*sqrt(1 + 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_317():
    f = x*sqrt(1 - x)*exp(acoth(x))
    F = 2*x**2*(1 + 1/x)**(sympy.S(3)/2)*sqrt(1 - x)/(5*sqrt(1 - 1/x)) - 4*x*(1 + 1/x)**(sympy.S(3)/2)*sqrt(1 - x)/(15*sqrt(1 - 1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_318():
    f = sqrt(1 - x)*exp(acoth(x))
    F = 2*sqrt(1 - x)*(x + 1)*exp(acoth(x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_319():
    f = x*exp(acoth(x))/sqrt(x + 1)
    F = 2*x**2*sqrt(-(1 - x)/x)*sqrt(1 + 1/x)/(3*sqrt(x + 1)) + 4*x*sqrt(-(1 - x)/x)*sqrt(1 + 1/x)/(3*sqrt(x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_320():
    f = exp(acoth(x))/sqrt(x + 1)
    F = 2*x*sqrt(-(1 - x)/x)*sqrt(1 + 1/x)/sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_321():
    f = x*exp(acoth(x))/sqrt(1 - x)
    F = 2*x**2*sqrt(1 - 1/x)*(1 + 1/x)**(sympy.S(3)/2)/(3*sqrt(1 - x)) + 2*x*sqrt(1 - 1/x)*sqrt(1 + 1/x)/sqrt(1 - x) - 2*sqrt(2)*sqrt(1 - 1/x)*atanh(sqrt(2)*sqrt(1/x)/sqrt(1 + 1/x))/(sqrt(1 - x)*sqrt(1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_322():
    f = exp(acoth(x))/sqrt(1 - x)
    F = 2*x*sqrt(1 - 1/x)*sqrt(1 + 1/x)/sqrt(1 - x) - 2*sqrt(2)*sqrt(1 - 1/x)*atanh(sqrt(2)*sqrt(1/x)/sqrt(1 + 1/x))/(sqrt(1 - x)*sqrt(1/x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_323():
    f = x*exp(acoth(x))/(x + 1)**(sympy.S(3)/2)
    F = 2*x**2*sqrt(-(1 - x)/x)*(1 + 1/x)**(sympy.S(3)/2)/(x + 1)**(sympy.S(3)/2) + sqrt(2)*(1 + 1/x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(1/x)/sqrt(-(1 - x)/x))/((x + 1)**(sympy.S(3)/2)*(1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_324():
    f = exp(acoth(x))/(x + 1)**(sympy.S(3)/2)
    F = -sqrt(2)*(1 + 1/x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(1/x)/sqrt(-(1 - x)/x))/((x + 1)**(sympy.S(3)/2)*(1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_325():
    f = x*exp(acoth(x))/(1 - x)**(sympy.S(3)/2)
    F = 5*x**2*(1 - 1/x)**(sympy.S(3)/2)*sqrt(1 + 1/x)/(2*(1 - x)**(sympy.S(3)/2)) - x**2*sqrt(1 - 1/x)*(1 + 1/x)**(sympy.S(3)/2)/(2*(1 - x)**(sympy.S(3)/2)) - 5*sqrt(2)*(1 - 1/x)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1/x)/sqrt(1 + 1/x))/(2*(1 - x)**(sympy.S(3)/2)*(1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_326():
    f = exp(acoth(x))/(1 - x)**(sympy.S(3)/2)
    F = -x*sqrt(1 - 1/x)*sqrt(1 + 1/x)/(1 - x)**(sympy.S(3)/2) - sqrt(2)*(1 - 1/x)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1/x)/sqrt(1 + 1/x))/(2*(1 - x)**(sympy.S(3)/2)*(1/x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_327():
    f = x**m*sqrt(-a*c*x + c)*exp(-acoth(a*x))
    F = 2*x**(m + 1)*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(sqrt(1 - 1/(a*x))*(2*m + 3)) - x**m*(8*m + 10)*sqrt(-a*c*x + c)*hyper((sympy.S.Half, -m + sympy.S(-1)/2), (sympy.S.Half - m,), -1/(a*x))/(a*sqrt(1 - 1/(a*x))*(2*m + 1)*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_328():
    f = sqrt(-a*c*x + c)*exp(-acoth(a*x))/x
    F = 2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/sqrt(1 - 1/(a*x)) + 2*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(sqrt(a)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_329():
    f = sqrt(-a*c*x + c)*exp(-acoth(a*x))/x**2
    F = -3*sqrt(a)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/sqrt(1 - 1/(a*x)) + sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(x*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_330():
    f = x**3*sqrt(-a*c*x + c)*exp(-2*acoth(a*x))
    F = -4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a**4 + 4*sqrt(-a*c*x + c)/a**4 + 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**4*c) + 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**4*c**2) - 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**4*c**3) + 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a**4*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_331():
    f = x**2*sqrt(-a*c*x + c)*exp(-2*acoth(a*x))
    F = 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a**3 - 4*sqrt(-a*c*x + c)/a**3 - 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**3*c) - 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**3*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_332():
    f = x*sqrt(-a*c*x + c)*exp(-2*acoth(a*x))
    F = -4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a**2 + 4*sqrt(-a*c*x + c)/a**2 + 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**2*c) + 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_333():
    f = sqrt(-a*c*x + c)*exp(-2*acoth(a*x))
    F = 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a - 4*sqrt(-a*c*x + c)/a - 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_334():
    f = sqrt(-a*c*x + c)*exp(-2*acoth(a*x))/x
    F = 2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) - 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) + 2*sqrt(-a*c*x + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_335():
    f = sqrt(-a*c*x + c)*exp(-2*acoth(a*x))/x**2
    F = -5*a*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) + 4*sqrt(2)*a*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) + sqrt(-a*c*x + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_336():
    f = sqrt(-a*c*x + c)*exp(-2*acoth(a*x))/x**3
    F = 23*a**2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/4 - 4*sqrt(2)*a**2*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) - 9*a*sqrt(-a*c*x + c)/(4*x) + sqrt(-a*c*x + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_337():
    f = sqrt(-a*c*x + c)*exp(-2*acoth(a*x))/x**4
    F = -45*a**3*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/8 + 4*sqrt(2)*a**3*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) + 19*a**2*sqrt(-a*c*x + c)/(8*x) - 13*a*sqrt(-a*c*x + c)/(12*x**2) + sqrt(-a*c*x + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_338():
    f = sqrt(-a*c*x + c)*exp(-2*acoth(a*x))/x**5
    F = 363*a**4*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/64 - 4*sqrt(2)*a**4*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) - 149*a**3*sqrt(-a*c*x + c)/(64*x) + 107*a**2*sqrt(-a*c*x + c)/(96*x**2) - 17*a*sqrt(-a*c*x + c)/(24*x**3) + sqrt(-a*c*x + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_339():
    f = x**3*sqrt(-a*c*x + c)*exp(-3*acoth(a*x))
    F = 2*x**4*sqrt(-a*c*x + c)/(9*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 8*x**3*sqrt(-a*c*x + c)/(9*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 164*x**2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(15*a**2*sqrt(1 - 1/(a*x))) - 82*x**2*sqrt(-a*c*x + c)/(9*a**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 656*x*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(45*a**3*sqrt(1 - 1/(a*x))) + 1312*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(45*a**4*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_340():
    f = x**2*sqrt(-a*c*x + c)*exp(-3*acoth(a*x))
    F = 2*x**3*sqrt(-a*c*x + c)/(7*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 44*x**2*sqrt(-a*c*x + c)/(35*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 1336*x*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(105*a**2*sqrt(1 - 1/(a*x))) - 334*x*sqrt(-a*c*x + c)/(35*a**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 2672*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(105*a**3*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_341():
    f = x*sqrt(-a*c*x + c)*exp(-3*acoth(a*x))
    F = 2*x**2*sqrt(-a*c*x + c)/(5*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 32*x*sqrt(-a*c*x + c)/(15*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 316*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(15*a**2*sqrt(1 - 1/(a*x))) - 158*sqrt(-a*c*x + c)/(15*a**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_342():
    f = sqrt(-a*c*x + c)*exp(-3*acoth(a*x))
    F = 2*x*sqrt(-a*c*x + c)/(3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 20*sqrt(-a*c*x + c)/(3*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 46*sqrt(-a*c*x + c)/(3*a**2*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_343():
    f = sqrt(-a*c*x + c)*exp(-3*acoth(a*x))/x
    F = 2*sqrt(-a*c*x + c)/(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 10*sqrt(-a*c*x + c)/(a*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 2*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(sqrt(a)*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_344():
    f = sqrt(-a*c*x + c)*exp(-3*acoth(a*x))/x**2
    F = 7*sqrt(a)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/sqrt(1 - 1/(a*x)) - sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(x*sqrt(1 - 1/(a*x))) - 8*sqrt(-a*c*x + c)/(x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_345():
    f = sqrt(-a*c*x + c)*exp(-3*acoth(a*x))/x**3
    F = -47*a**(sympy.S(3)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(4*sqrt(1 - 1/(a*x))) + 47*a*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(4*x*sqrt(1 - 1/(a*x))) - sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(2*x**2*sqrt(1 - 1/(a*x))) - 8*sqrt(-a*c*x + c)/(x**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_346():
    f = sqrt(-a*c*x + c)*exp(-3*acoth(a*x))/x**4
    F = 119*a**(sympy.S(5)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(8*sqrt(1 - 1/(a*x))) - 119*a**2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(8*x*sqrt(1 - 1/(a*x))) + 119*a*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(12*x**2*sqrt(1 - 1/(a*x))) - sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(3*x**3*sqrt(1 - 1/(a*x))) - 8*sqrt(-a*c*x + c)/(x**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_347():
    f = sqrt(-a*c*x + c)*exp(-3*acoth(a*x))/x**5
    F = -1115*a**(sympy.S(7)/2)*sqrt(-a*c*x + c)*sqrt(1/x)*asinh(sqrt(1/x)/sqrt(a))/(64*sqrt(1 - 1/(a*x))) + 1115*a**3*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(64*x*sqrt(1 - 1/(a*x))) - 1115*a**2*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(96*x**2*sqrt(1 - 1/(a*x))) + 223*a*sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(24*x**3*sqrt(1 - 1/(a*x))) - sqrt(1 + 1/(a*x))*sqrt(-a*c*x + c)/(4*x**4*sqrt(1 - 1/(a*x))) - 8*sqrt(-a*c*x + c)/(x**4*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_348():
    f = (-a*c*x + c)**(n/2 + 2)*exp(n*acoth(a*x))
    F = x*(1 - 1/(a*x))**(-n/2 - 2)*(1 + 1/(a*x))**(n/2 + 1)*(n + 8)*(-a*c*x + c)**(n/2 + 2)/(n + 6) - x*(1 - 1/(a*x))**(-n/2 - 2)*(1 + 1/(a*x))**(n/2 + 1)*(a - 1/x)*(-a*c*x + c)**(n/2 + 2)/a - (1 - 1/(a*x))**(-n/2 - 2)*(1 + 1/(a*x))**(n/2 + 1)*(-a*c*x + c)**(n/2 + 2)*(n**2 + 14*n + 56)/(a*(n + 4)*(n + 6)) + (1 - 1/(a*x))**(-n/2 - 2)*(1 + 1/(a*x))**(n/2 + 1)*(-a*c*x + c)**(n/2 + 2)*(2*n**2 + 28*n + 112)/(a**2*x*(n + 6)*(n**2 + 6*n + 8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_349():
    f = (-a*c*x + c)**(n/2 + 1)*exp(n*acoth(a*x))
    F = 2*x*(1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 + 1)*(-a*c*x + c)**(n/2 + 1)/(n + 4) - (1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 + 1)*(2*n + 12)*(-a*c*x + c)**(n/2 + 1)/(a*(n + 2)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_350():
    f = (-a*c*x + c)**(n/2)*exp(n*acoth(a*x))
    F = 2*(a*x + 1)*(-a*c*x + c)**(n/2)*exp(n*acoth(a*x))/(a*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_351():
    f = (-a*c*x + c)**(n/2 - 1)*exp(n*acoth(a*x))
    F = 2*x*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2)*(-a*c*x + c)**(n/2 - 1)*hyper((1, -n/2), (1 - n/2,), 2/(x*(a + 1/x)))/n
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_352():
    f = (-a*c*x + c)**(n/2 - 2)*exp(n*acoth(a*x))
    F = -2*x*(1 - 1/(a*x))**(2 - n/2)*(1 + 1/(a*x))**(n/2 - 1)*(-a*c*x + c)**(n/2 - 2)*hyper((2, 1 - n/2), (2 - n/2,), 2/(x*(a + 1/x)))/(2 - n)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_353():
    f = (-a*c*x + c)**p*exp(n*acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(n/2 - p)*(1 + 1/(a*x))**(n/2 + 1)*(-a*c*x + c)**p*hyper((-p - 1, n/2 - p), (-p,), 2/(x*(a + 1/x)))/((1 - 1/(a*x))**(n/2)*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_354():
    f = (-a*c*x + c)**3*exp(n*acoth(a*x))
    F = -32*c**3*(1 - 1/(a*x))**(4 - n/2)*(1 + 1/(a*x))**(n/2 - 4)*hyper((5, 4 - n/2), (5 - n/2,), (a - 1/x)/(a + 1/x))/(a*(8 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_355():
    f = (-a*c*x + c)**2*exp(n*acoth(a*x))
    F = 16*c**2*(1 - 1/(a*x))**(3 - n/2)*(1 + 1/(a*x))**(n/2 - 3)*hyper((4, 3 - n/2), (4 - n/2,), (a - 1/x)/(a + 1/x))/(a*(6 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_356():
    f = (-a*c*x + c)*exp(n*acoth(a*x))
    F = -8*c*(1 - 1/(a*x))**(2 - n/2)*(1 + 1/(a*x))**(n/2 - 2)*hyper((3, 2 - n/2), (3 - n/2,), (a - 1/x)/(a + 1/x))/(a*(4 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_357():
    f = exp(n*acoth(a*x))/(-a*c*x + c)
    F = 2*(1 + 1/(a*x))**(n/2)*hyper((1, -n/2), (1 - n/2,), (a - 1/x)/(a + 1/x))/(a*c*n*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_358():
    f = exp(n*acoth(a*x))/(-a*c*x + c)**2
    F = -(1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 + 1)/(a*c**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_359():
    f = exp(n*acoth(a*x))/(-a*c*x + c)**3
    F = (1 - 1/(a*x))**(-n/2 - 2)*(1 + 1/(a*x))**(n/2 + 1)/(a*c**3*(n + 4)) - (1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 + 1)*(n + 3)/(a*c**3*(n + 2)*(n + 4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_360():
    f = exp(n*acoth(a*x))/(-a*c*x + c)**4
    F = (1 - 1/(a*x))**(-n/2 - 3)*(1 + 1/(a*x))**(n/2 + 1)*(n + 5)/(a*c**4*(n + 6)) - (1 - 1/(a*x))**(-n/2 - 2)*(1 + 1/(a*x))**(n/2 + 1)*(n**2 + 8*n + 14)/(a*c**4*(n + 4)*(n + 6)) - (1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 + 1)*(n**2 + 8*n + 14)/(a*c**4*(n + 6)*(n**2 + 6*n + 8)) - (1 - 1/(a*x))**(-n/2 - 3)*(1 + 1/(a*x))**(n/2 + 1)/(a**2*c**4*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_361():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(n*acoth(a*x))
    F = 2*x*((a - 1/x)/(a + 1/x))**(n/2 + sympy.S(-5)/2)*(1 + 1/(a*x))**(n/2 + 1)*(-a*c*x + c)**(sympy.S(5)/2)*hyper((sympy.S(-7)/2, n/2 + sympy.S(-5)/2), (sympy.S(-5)/2,), 2/(x*(a + 1/x)))/(7*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_362():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(n*acoth(a*x))
    F = 2*x*((a - 1/x)/(a + 1/x))**(n/2 + sympy.S(-3)/2)*(1 + 1/(a*x))**(n/2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)*hyper((sympy.S(-5)/2, n/2 + sympy.S(-3)/2), (sympy.S(-3)/2,), 2/(x*(a + 1/x)))/(5*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_363():
    f = sqrt(-a*c*x + c)*exp(n*acoth(a*x))
    F = 2*x*((a - 1/x)/(a + 1/x))**(n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + 1)*sqrt(-a*c*x + c)*hyper((sympy.S(-3)/2, n/2 + sympy.S(-1)/2), (sympy.S(-1)/2,), 2/(x*(a + 1/x)))/(3*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_364():
    f = exp(n*acoth(a*x))/sqrt(-a*c*x + c)
    F = 2*x*((a - 1/x)/(a + 1/x))**(n/2 + sympy.S.Half)*(1 + 1/(a*x))**(n/2 + 1)*hyper((sympy.S(-1)/2, n/2 + sympy.S.Half), (sympy.S.Half,), 2/(x*(a + 1/x)))/((1 - 1/(a*x))**(n/2)*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_365():
    f = exp(n*acoth(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -2*x*((a - 1/x)/(a + 1/x))**(n/2 + sympy.S(3)/2)*(1 + 1/(a*x))**(n/2 + 1)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (sympy.S(3)/2,), 2/(x*(a + 1/x)))/((1 - 1/(a*x))**(n/2)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_366():
    f = exp(n*acoth(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = a*x**2*((a - 1/x)/(a + 1/x))**(n/2 + sympy.S(3)/2)*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (sympy.S(3)/2,), 2/(x*(a + 1/x)))/((n + 3)*(-a*c*x + c)**(sympy.S(5)/2)) - a*x**2*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/((n + 3)*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_367():
    f = exp(n*acoth(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = -3*a**2*x**3*((a - 1/x)/(a + 1/x))**(n/2 + sympy.S(3)/2)*(1 - 1/(a*x))**(2 - n/2)*(1 + 1/(a*x))**(n/2 + 1)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (sympy.S(3)/2,), 2/(x*(a + 1/x)))/((-a*c*x + c)**(sympy.S(7)/2)*(2*n**2 + 16*n + 30)) + 3*a**2*x**3*(1 - 1/(a*x))**(2 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/((-a*c*x + c)**(sympy.S(7)/2)*(2*n**2 + 16*n + 30)) - a*x**2*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 + 1)/((n + 5)*(-a*c*x + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_368():
    f = (c - c/(a*x))**4*exp(acoth(a*x))
    F = c**4*x*(1 - 1/(a**2*x**2))**(sympy.S(3)/2) - c**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(3*a) - c**4*acsc(a*x)/(2*a) - 3*c**4*atanh(sqrt(1 - 1/(a**2*x**2)))/a + c**4*sqrt(1 - 1/(a**2*x**2))*(6*a - 1/x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_369():
    f = (c - c/(a*x))**3*exp(acoth(a*x))
    F = c**3*x*(1 - 1/(a**2*x**2))**(sympy.S(3)/2) + c**3*acsc(a*x)/(2*a) - 2*c**3*atanh(sqrt(1 - 1/(a**2*x**2)))/a + c**3*sqrt(1 - 1/(a**2*x**2))*(4*a + 1/x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_370():
    f = (c - c/(a*x))**2*exp(acoth(a*x))
    F = c**2*x*sqrt(1 - 1/(a**2*x**2))*(a + 1/x)/a + c**2*acsc(a*x)/a - c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_371():
    f = (c - c/(a*x))*exp(acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2)) + c*acsc(a*x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_372():
    f = exp(acoth(a*x))/(c - c/(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2))/c + 2*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c) - (2*a + 2/x)/(a**2*c*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_373():
    f = exp(acoth(a*x))/(c - c/(a*x))**2
    F = x*sqrt(1 - 1/(a**2*x**2))/c**2 + 3*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**2) - (9*a + 11/x)/(3*a**2*c**2*sqrt(1 - 1/(a**2*x**2))) - (4*a + 4/x)/(3*a**2*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_374():
    f = exp(acoth(a*x))/(c - c/(a*x))**3
    F = x*sqrt(1 - 1/(a**2*x**2))/c**3 + 4*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**3) - (60*a + 79/x)/(15*a**2*c**3*sqrt(1 - 1/(a**2*x**2))) - (20*a + 32/x)/(15*a**2*c**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - (8*a + 8/x)/(5*a**2*c**3*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_375():
    f = exp(acoth(a*x))/(c - c/(a*x))**4
    F = x*sqrt(1 - 1/(a**2*x**2))/c**4 + 5*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**4) - (525*a + 719/x)/(105*a**2*c**4*sqrt(1 - 1/(a**2*x**2))) - (175*a + 307/x)/(105*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - (28*a + 68/x)/(35*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) - (16*a + 16/x)/(7*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_376():
    f = (c - c/(a*x))**5*exp(2*acoth(a*x))
    F = c**5*x - 3*c**5*log(x)/a - 2*c**5/(a**2*x) - c**5/(a**3*x**2) + c**5/(a**4*x**3) - c**5/(4*a**5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_377():
    f = (c - c/(a*x))**4*exp(2*acoth(a*x))
    F = c**4*x - 2*c**4*log(x)/a - c**4/(a**3*x**2) + c**4/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_378():
    f = (c - c/(a*x))**3*exp(2*acoth(a*x))
    F = c**3*x - c**3*log(x)/a + c**3/(a**2*x) - c**3/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_379():
    f = (c - c/(a*x))**2*exp(2*acoth(a*x))
    F = c**2*x + c**2/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_380():
    f = (c - c/(a*x))*exp(2*acoth(a*x))
    F = c*x + c*log(x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_381():
    f = exp(2*acoth(a*x))/(c - c/(a*x))
    F = x/c + 3*log(-a*x + 1)/(a*c) + 2/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_382():
    f = exp(2*acoth(a*x))/(c - c/(a*x))**2
    F = x/c**2 + 4*log(-a*x + 1)/(a*c**2) + 5/(a*c**2*(-a*x + 1)) - 1/(a*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_383():
    f = exp(2*acoth(a*x))/(c - c/(a*x))**3
    F = x/c**3 + 5*log(-a*x + 1)/(a*c**3) + 9/(a*c**3*(-a*x + 1)) - 7/(2*a*c**3*(-a*x + 1)**2) + 2/(3*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_384():
    f = exp(2*acoth(a*x))/(c - c/(a*x))**4
    F = x/c**4 + 6*log(-a*x + 1)/(a*c**4) + 14/(a*c**4*(-a*x + 1)) - 8/(a*c**4*(-a*x + 1)**2) + 3/(a*c**4*(-a*x + 1)**3) - 1/(2*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_385():
    f = (c - c/(a*x))**4*exp(3*acoth(a*x))
    F = c**4*x*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(3*a + 1/x)/(3*a) + 3*c**4*acsc(a*x)/(2*a) - c**4*atanh(sqrt(1 - 1/(a**2*x**2)))/a + c**4*sqrt(1 - 1/(a**2*x**2))*(2*a + 3/x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_386():
    f = (c - c/(a*x))**3*exp(3*acoth(a*x))
    F = c**3*x*(1 - 1/(a**2*x**2))**(sympy.S(3)/2) + 3*c**3*acsc(a*x)/(2*a) + 3*c**3*sqrt(1 - 1/(a**2*x**2))/(2*a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_387():
    f = (c - c/(a*x))**2*exp(3*acoth(a*x))
    F = c**2*x*sqrt(1 - 1/(a**2*x**2))*(a - 1/x)/a + c**2*acsc(a*x)/a + c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_388():
    f = (c - c/(a*x))*exp(3*acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2)) - c*acsc(a*x)/a + 2*c*atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_389():
    f = exp(3*acoth(a*x))/(c - c/(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2))/c + 4*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c) - (12*a + 16/x)/(3*a**2*c*sqrt(1 - 1/(a**2*x**2))) - (8*a + 8/x)/(3*a**2*c*(1 - 1/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_390():
    f = exp(3*acoth(a*x))/(c - c/(a*x))**2
    F = x*sqrt(1 - 1/(a**2*x**2))/c**2 + 5*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**2) - (75*a + 103/x)/(15*a**2*c**2*sqrt(1 - 1/(a**2*x**2))) - (20*a + 44/x)/(15*a**2*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - (16*a + 16/x)/(5*a**2*c**2*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_391():
    f = exp(3*acoth(a*x))/(c - c/(a*x))**3
    F = x*sqrt(1 - 1/(a**2*x**2))/c**3 + 6*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**3) - (42*a + 59/x)/(7*a**2*c**3*sqrt(1 - 1/(a**2*x**2))) - (14*a + 26/x)/(7*a**2*c**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - (32*a + 32/x)/(7*a**2*c**3*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 16/(7*a**2*c**3*x*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_392():
    f = exp(3*acoth(a*x))/(c - c/(a*x))**4
    F = x*sqrt(1 - 1/(a**2*x**2))/c**4 + 7*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**4) - (2205*a + 3149/x)/(315*a**2*c**4*sqrt(1 - 1/(a**2*x**2))) - (735*a + 1417/x)/(315*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - (168*a + 328/x)/(105*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) + (144*a - 80/x)/(63*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - (64*a + 64/x)/(9*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_393():
    f = (c - c/(a*x))**5*exp(4*acoth(a*x))
    F = c**5*x - c**5*log(x)/a + 2*c**5/(a**2*x) - c**5/(a**3*x**2) - c**5/(3*a**4*x**3) + c**5/(4*a**5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_394():
    f = (c - c/(a*x))**4*exp(4*acoth(a*x))
    F = c**4*x + 2*c**4/(a**2*x) - c**4/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_395():
    f = (c - c/(a*x))**3*exp(4*acoth(a*x))
    F = c**3*x + c**3*log(x)/a + c**3/(a**2*x) + c**3/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_396():
    f = (c - c/(a*x))**2*exp(4*acoth(a*x))
    F = c**2*x + 2*c**2*log(x)/a - c**2/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_397():
    f = (c - c/(a*x))*exp(4*acoth(a*x))
    F = c*x - c*log(x)/a + 4*c*log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_398():
    f = exp(4*acoth(a*x))/(c - c/(a*x))
    F = x/c + 5*log(-a*x + 1)/(a*c) + 8/(a*c*(-a*x + 1)) - 2/(a*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_399():
    f = exp(4*acoth(a*x))/(c - c/(a*x))**2
    F = x/c**2 + 6*log(-a*x + 1)/(a*c**2) + 13/(a*c**2*(-a*x + 1)) - 6/(a*c**2*(-a*x + 1)**2) + 4/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_400():
    f = exp(4*acoth(a*x))/(c - c/(a*x))**3
    F = x/c**3 + 7*log(-a*x + 1)/(a*c**3) + 19/(a*c**3*(-a*x + 1)) - 25/(2*a*c**3*(-a*x + 1)**2) + 16/(3*a*c**3*(-a*x + 1)**3) - 1/(a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_401():
    f = exp(4*acoth(a*x))/(c - c/(a*x))**4
    F = x/c**4 + 8*log(-a*x + 1)/(a*c**4) + 26/(a*c**4*(-a*x + 1)) - 22/(a*c**4*(-a*x + 1)**2) + 41/(3*a*c**4*(-a*x + 1)**3) - 5/(a*c**4*(-a*x + 1)**4) + 4/(5*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_402():
    f = (c - c/(a*x))**4*exp(-acoth(a*x))
    F = c**4*x*sqrt(1 - 1/(a**2*x**2)) - 32*c**4*sqrt(1 - 1/(a**2*x**2))/(3*a) - 25*c**4*acsc(a*x)/(2*a) - 5*c**4*atanh(sqrt(1 - 1/(a**2*x**2)))/a + 5*c**4*sqrt(1 - 1/(a**2*x**2))/(2*a**2*x) - c**4*sqrt(1 - 1/(a**2*x**2))/(3*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_403():
    f = (c - c/(a*x))**3*exp(-acoth(a*x))
    F = c**3*x*sqrt(1 - 1/(a**2*x**2)) - 4*c**3*sqrt(1 - 1/(a**2*x**2))/a - 13*c**3*acsc(a*x)/(2*a) - 4*c**3*atanh(sqrt(1 - 1/(a**2*x**2)))/a + c**3*sqrt(1 - 1/(a**2*x**2))/(2*a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_404():
    f = (c - c/(a*x))**2*exp(-acoth(a*x))
    F = c**2*x*sqrt(1 - 1/(a**2*x**2)) - c**2*sqrt(1 - 1/(a**2*x**2))/a - 3*c**2*acsc(a*x)/a - 3*c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_405():
    f = (c - c/(a*x))*exp(-acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2)) - c*acsc(a*x)/a - 2*c*atanh(sqrt(1 - 1/(a**2*x**2)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_406():
    f = exp(-acoth(a*x))/(c - c/(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2))/c
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_407():
    f = exp(-acoth(a*x))/(c - c/(a*x))**2
    F = -a*x*sqrt(1 - 1/(a**2*x**2))/(c**2*(a - 1/x)) + 2*x*sqrt(1 - 1/(a**2*x**2))/c**2 + atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_408():
    f = exp(-acoth(a*x))/(c - c/(a*x))**3
    F = x*sqrt(1 - 1/(a**2*x**2))/c**3 + 2*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**3) - (6*a + 7/x)/(3*a**2*c**3*sqrt(1 - 1/(a**2*x**2))) - (2*a + 2/x)/(3*a**2*c**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_409():
    f = exp(-acoth(a*x))/(c - c/(a*x))**4
    F = x*sqrt(1 - 1/(a**2*x**2))/c**4 + 3*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**4) - (15*a + 19/x)/(5*a**2*c**4*sqrt(1 - 1/(a**2*x**2))) - (5*a + 7/x)/(5*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - (4*a + 4/x)/(5*a**2*c**4*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_410():
    f = (c - c/(a*x))**4*exp(-2*acoth(a*x))
    F = c**4*x + 26*c**4*log(x)/a - 32*c**4*log(a*x + 1)/a + 16*c**4/(a**2*x) - 3*c**4/(a**3*x**2) + c**4/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_411():
    f = (c - c/(a*x))**3*exp(-2*acoth(a*x))
    F = c**3*x + 11*c**3*log(x)/a - 16*c**3*log(a*x + 1)/a + 5*c**3/(a**2*x) - c**3/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_412():
    f = (c - c/(a*x))**2*exp(-2*acoth(a*x))
    F = c**2*x + 4*c**2*log(x)/a - 8*c**2*log(a*x + 1)/a + c**2/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_413():
    f = (c - c/(a*x))*exp(-2*acoth(a*x))
    F = c*x + c*log(x)/a - 4*c*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_414():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))
    F = x/c - log(a*x + 1)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_415():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))**2
    F = x/c**2 - atanh(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_416():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))**3
    F = x/c**3 + 5*log(-a*x + 1)/(4*a*c**3) - log(a*x + 1)/(4*a*c**3) + 1/(2*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_417():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))**4
    F = x/c**4 + 17*log(-a*x + 1)/(8*a*c**4) - log(a*x + 1)/(8*a*c**4) + 7/(4*a*c**4*(-a*x + 1)) - 1/(4*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_418():
    f = (c - c/(a*x))**4*exp(-3*acoth(a*x))
    F = c**4*x*sqrt(1 - 1/(a**2*x**2)) + 68*c**4*sqrt(1 - 1/(a**2*x**2))/(3*a) + 91*c**4*acsc(a*x)/(2*a) - 7*c**4*atanh(sqrt(1 - 1/(a**2*x**2)))/a + 64*c**4*(a - 1/x)/(a**2*sqrt(1 - 1/(a**2*x**2))) - 7*c**4*sqrt(1 - 1/(a**2*x**2))/(2*a**2*x) + c**4*sqrt(1 - 1/(a**2*x**2))/(3*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_419():
    f = (c - c/(a*x))**3*exp(-3*acoth(a*x))
    F = c**3*x*sqrt(1 - 1/(a**2*x**2)) + 6*c**3*sqrt(1 - 1/(a**2*x**2))/a + 33*c**3*acsc(a*x)/(2*a) - 6*c**3*atanh(sqrt(1 - 1/(a**2*x**2)))/a + 32*c**3*(a - 1/x)/(a**2*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(1 - 1/(a**2*x**2))/(2*a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_420():
    f = (c - c/(a*x))**2*exp(-3*acoth(a*x))
    F = c**2*x*sqrt(1 - 1/(a**2*x**2)) + c**2*sqrt(1 - 1/(a**2*x**2))/a + 5*c**2*acsc(a*x)/a - 5*c**2*atanh(sqrt(1 - 1/(a**2*x**2)))/a + 16*c**2*(a - 1/x)/(a**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_421():
    f = (c - c/(a*x))*exp(-3*acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2)) + c*acsc(a*x)/a - 4*c*atanh(sqrt(1 - 1/(a**2*x**2)))/a + 8*c*(a - 1/x)/(a**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_422():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2))/c - 2*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c) + (2*a - 2/x)/(a**2*c*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_423():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))**2
    F = 2*x*sqrt(1 - 1/(a**2*x**2))/c**2 - x*(a - 1/x)/(a*c**2*sqrt(1 - 1/(a**2*x**2))) - atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_424():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))**3
    F = x/(c**3*sqrt(1 - 1/(a**2*x**2))) - 2/(a**2*c**3*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_425():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))**4
    F = -a*x/(3*c**4*sqrt(1 - 1/(a**2*x**2))*(a - 1/x)) + 8*x*sqrt(1 - 1/(a**2*x**2))/(3*c**4) - x*(4*a + 3/x)/(3*a*c**4*sqrt(1 - 1/(a**2*x**2))) + atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_426():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))**5
    F = x*sqrt(1 - 1/(a**2*x**2))/c**5 + 2*atanh(sqrt(1 - 1/(a**2*x**2)))/(a*c**5) - (30*a + 41/x)/(15*a**2*c**5*sqrt(1 - 1/(a**2*x**2))) - (10*a + 13/x)/(15*a**2*c**5*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - (2*a + 2/x)/(5*a**2*c**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_427():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(acoth(a*x))
    F = c**3*x*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(c - c/(a*x))**(sympy.S(3)/2) - c**(sympy.S(3)/2)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/a + c**2*sqrt(1 - 1/(a**2*x**2))/(a*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_428():
    f = sqrt(c - c/(a*x))*exp(acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)) + sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_429():
    f = exp(acoth(a*x))/sqrt(c - c/(a*x))
    F = x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/sqrt(c - c/(a*x)) - 2*sqrt(2)*sqrt(1 - 1/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a*sqrt(c - c/(a*x))) + 3*sqrt(1 - 1/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(a*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_430():
    f = exp(acoth(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = a*x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)*(c - c/(a*x))**(sympy.S(3)/2)) - 2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)*(c - c/(a*x))**(sympy.S(3)/2)) - 7*sqrt(2)*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(2*a*(c - c/(a*x))**(sympy.S(3)/2)) + 5*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_431():
    f = exp(acoth(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = a**2*x*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)**2*(c - c/(a*x))**(sympy.S(5)/2)) - 3*a*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/(2*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(5)/2)) - 23*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((8*a - 8/x)*(c - c/(a*x))**(sympy.S(5)/2)) - 79*sqrt(2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(16*a*(c - c/(a*x))**(sympy.S(5)/2)) + 7*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_432():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(2*acoth(a*x))
    F = x*(c - c/(a*x))**(sympy.S(9)/2) - 5*c**(sympy.S(9)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + 5*c**4*sqrt(c - c/(a*x))/a + 5*c**3*(c - c/(a*x))**(sympy.S(3)/2)/(3*a) + c**2*(c - c/(a*x))**(sympy.S(5)/2)/a + 5*c*(c - c/(a*x))**(sympy.S(7)/2)/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_433():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(2*acoth(a*x))
    F = x*(c - c/(a*x))**(sympy.S(7)/2) - 3*c**(sympy.S(7)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + 3*c**3*sqrt(c - c/(a*x))/a + c**2*(c - c/(a*x))**(sympy.S(3)/2)/a + 3*c*(c - c/(a*x))**(sympy.S(5)/2)/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_434():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(2*acoth(a*x))
    F = x*(c - c/(a*x))**(sympy.S(5)/2) - c**(sympy.S(5)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + c**2*sqrt(c - c/(a*x))/a + c*(c - c/(a*x))**(sympy.S(3)/2)/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_435():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(2*acoth(a*x))
    F = x*(c - c/(a*x))**(sympy.S(3)/2) + c**(sympy.S(3)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - c*sqrt(c - c/(a*x))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_436():
    f = sqrt(c - c/(a*x))*exp(2*acoth(a*x))
    F = x*sqrt(c - c/(a*x)) + 3*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_437():
    f = exp(2*acoth(a*x))/sqrt(c - c/(a*x))
    F = x/sqrt(c - c/(a*x)) - 5/(a*sqrt(c - c/(a*x))) + 5*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_438():
    f = exp(2*acoth(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = x/(c - c/(a*x))**(sympy.S(3)/2) - 7/(3*a*(c - c/(a*x))**(sympy.S(3)/2)) - 7/(a*c*sqrt(c - c/(a*x))) + 7*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_439():
    f = exp(2*acoth(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = x/(c - c/(a*x))**(sympy.S(5)/2) - 9/(5*a*(c - c/(a*x))**(sympy.S(5)/2)) - 3/(a*c*(c - c/(a*x))**(sympy.S(3)/2)) - 9/(a*c**2*sqrt(c - c/(a*x))) + 9*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_440():
    f = exp(2*acoth(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = x/(c - c/(a*x))**(sympy.S(7)/2) - 11/(7*a*(c - c/(a*x))**(sympy.S(7)/2)) - 11/(5*a*c*(c - c/(a*x))**(sympy.S(5)/2)) - 11/(3*a*c**2*(c - c/(a*x))**(sympy.S(3)/2)) - 11/(a*c**3*sqrt(c - c/(a*x))) + 11*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_441():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(3*acoth(a*x))
    F = 3*sqrt(1 + 1/(a*x))*(c - c/(a*x))**(sympy.S(9)/2)/(a*(1 - 1/(a*x))**(sympy.S(9)/2)) - 3*(c - c/(a*x))**(sympy.S(9)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(9)/2)) + (1 + 1/(a*x))**(sympy.S(3)/2)*(84*a - 51/x)*(c - c/(a*x))**(sympy.S(9)/2)/(35*a**2*(1 - 1/(a*x))**(sympy.S(9)/2)) + x*(1 + 1/(a*x))**(sympy.S(3)/2)*(a - 1/x)**3*(c - c/(a*x))**(sympy.S(9)/2)/(a**3*(1 - 1/(a*x))**(sympy.S(9)/2)) + 9*(1 + 1/(a*x))**(sympy.S(3)/2)*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(9)/2)/(7*a**3*(1 - 1/(a*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_442():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(3*acoth(a*x))
    F = x*(1 + 1/(a*x))**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(7)/2)/(1 - 1/(a*x))**(sympy.S(7)/2) - 2*(1 + 1/(a*x))**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(7)/2)/(5*a*(1 - 1/(a*x))**(sympy.S(7)/2)) + (1 + 1/(a*x))**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(7)/2)/(3*a*(1 - 1/(a*x))**(sympy.S(7)/2)) + sqrt(1 + 1/(a*x))*(c - c/(a*x))**(sympy.S(7)/2)/(a*(1 - 1/(a*x))**(sympy.S(7)/2)) - (c - c/(a*x))**(sympy.S(7)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_443():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(3*acoth(a*x))
    F = c**5*x*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(c - c/(a*x))**(sympy.S(5)/2) + c**(sympy.S(5)/2)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/a - c**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(3*a*(c - c/(a*x))**(sympy.S(3)/2)) - c**3*sqrt(1 - 1/(a**2*x**2))/(a*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_444():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(3*acoth(a*x))
    F = c**3*x*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(c - c/(a*x))**(sympy.S(3)/2) + 3*c**(sympy.S(3)/2)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/a - 3*c**2*sqrt(1 - 1/(a**2*x**2))/(a*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_445():
    f = sqrt(c - c/(a*x))*exp(3*acoth(a*x))
    F = x*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/sqrt(1 - 1/(a*x)) - 4*sqrt(2)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a*sqrt(1 - 1/(a*x))) + 5*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(a*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_446():
    f = exp(3*acoth(a*x))/sqrt(c - c/(a*x))
    F = a*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/((a - 1/x)*sqrt(c - c/(a*x))) - 3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/((a - 1/x)*sqrt(c - c/(a*x))) - 5*sqrt(2)*sqrt(1 - 1/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a*sqrt(c - c/(a*x))) + 7*sqrt(1 - 1/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(a*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_447():
    f = exp(3*acoth(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = a**2*x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)**2*(c - c/(a*x))**(sympy.S(3)/2)) - 2*a*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)**2*(c - c/(a*x))**(sympy.S(3)/2)) - 15*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/((4*a - 4/x)*(c - c/(a*x))**(sympy.S(3)/2)) - 51*sqrt(2)*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(8*a*(c - c/(a*x))**(sympy.S(3)/2)) + 9*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_448():
    f = exp(3*acoth(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = a**3*x*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)**3*(c - c/(a*x))**(sympy.S(5)/2)) - 5*a**2*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/(3*(a - 1/x)**3*(c - c/(a*x))**(sympy.S(5)/2)) - 29*a*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/(12*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(5)/2)) - 73*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((16*a - 16/x)*(c - c/(a*x))**(sympy.S(5)/2)) - 249*sqrt(2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(32*a*(c - c/(a*x))**(sympy.S(5)/2)) + 11*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_449():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(-acoth(a*x))
    F = -9*(c - c/(a*x))**(sympy.S(7)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(7)/2)) - sqrt(1 + 1/(a*x))*(80*a - 7/x)*(c - c/(a*x))**(sympy.S(7)/2)/(5*a**2*(1 - 1/(a*x))**(sympy.S(7)/2)) + x*sqrt(1 + 1/(a*x))*(a - 1/x)**3*(c - c/(a*x))**(sympy.S(7)/2)/(a**3*(1 - 1/(a*x))**(sympy.S(7)/2)) + 3*sqrt(1 + 1/(a*x))*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(7)/2)/(5*a**3*(1 - 1/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_450():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(-acoth(a*x))
    F = -7*(c - c/(a*x))**(sympy.S(5)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(5)/2)) + x*sqrt(1 + 1/(a*x))*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(5)/2)/(a**2*(1 - 1/(a*x))**(sympy.S(5)/2)) - sqrt(1 + 1/(a*x))*(16*a + 1/x)*(c - c/(a*x))**(sympy.S(5)/2)/(3*a**2*(1 - 1/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_451():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(-acoth(a*x))
    F = x*sqrt(1 + 1/(a*x))*(c - c/(a*x))**(sympy.S(3)/2)/(1 - 1/(a*x))**(sympy.S(3)/2) - 2*sqrt(1 + 1/(a*x))*(c - c/(a*x))**(sympy.S(3)/2)/(a*(1 - 1/(a*x))**(sympy.S(3)/2)) - 5*(c - c/(a*x))**(sympy.S(3)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_452():
    f = sqrt(c - c/(a*x))*exp(-acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)) - 3*sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_453():
    f = exp(-acoth(a*x))/sqrt(c - c/(a*x))
    F = x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)) - atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_454():
    f = exp(-acoth(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/(c - c/(a*x))**(sympy.S(3)/2) - sqrt(2)*(1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a*(c - c/(a*x))**(sympy.S(3)/2)) + (1 - 1/(a*x))**(sympy.S(3)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_455():
    f = exp(-acoth(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = a*x*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)*(c - c/(a*x))**(sympy.S(5)/2)) - 3*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/((2*a - 2/x)*(c - c/(a*x))**(sympy.S(5)/2)) - 9*sqrt(2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(4*a*(c - c/(a*x))**(sympy.S(5)/2)) + 3*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_456():
    f = exp(-acoth(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = a**2*x*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))/((a - 1/x)**2*(c - c/(a*x))**(sympy.S(7)/2)) - 5*a*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))/(4*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(7)/2)) - 35*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))/((16*a - 16/x)*(c - c/(a*x))**(sympy.S(7)/2)) - 115*sqrt(2)*(1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(32*a*(c - c/(a*x))**(sympy.S(7)/2)) + 5*(1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_457():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(-2*acoth(a*x))
    F = x*(c - c/(a*x))**(sympy.S(7)/2) - 11*c**(sympy.S(7)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + 32*sqrt(2)*c**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a - 21*c**3*sqrt(c - c/(a*x))/a - 5*c**2*(c - c/(a*x))**(sympy.S(3)/2)/(3*a) + 3*c*(c - c/(a*x))**(sympy.S(5)/2)/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_458():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(-2*acoth(a*x))
    F = x*(c - c/(a*x))**(sympy.S(5)/2) - 9*c**(sympy.S(5)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + 16*sqrt(2)*c**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a - 7*c**2*sqrt(c - c/(a*x))/a + c*(c - c/(a*x))**(sympy.S(3)/2)/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_459():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(-2*acoth(a*x))
    F = x*(c - c/(a*x))**(sympy.S(3)/2) - 7*c**(sympy.S(3)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + 8*sqrt(2)*c**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a - c*sqrt(c - c/(a*x))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_460():
    f = sqrt(c - c/(a*x))*exp(-2*acoth(a*x))
    F = x*sqrt(c - c/(a*x)) - 5*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_461():
    f = exp(-2*acoth(a*x))/sqrt(c - c/(a*x))
    F = x*sqrt(c - c/(a*x))/c - 3*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*sqrt(c)) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_462():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = x*sqrt(c - c/(a*x))/c**2 - atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_463():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = x/(c**2*sqrt(c - c/(a*x))) - 2/(a*c**2*sqrt(c - c/(a*x))) + atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(5)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(2*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_464():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = x/(c**2*(c - c/(a*x))**(sympy.S(3)/2)) - 4/(3*a*c**2*(c - c/(a*x))**(sympy.S(3)/2)) - 7/(2*a*c**3*sqrt(c - c/(a*x))) + 3*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(7)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(4*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_465():
    f = exp(-2*acoth(a*x))/(c - c/(a*x))**(sympy.S(9)/2)
    F = x/(c**2*(c - c/(a*x))**(sympy.S(5)/2)) - 6/(5*a*c**2*(c - c/(a*x))**(sympy.S(5)/2)) - 11/(6*a*c**3*(c - c/(a*x))**(sympy.S(3)/2)) - 21/(4*a*c**4*sqrt(c - c/(a*x))) + 5*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(9)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(8*a*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_466():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(-3*acoth(a*x))
    F = -15*(c - c/(a*x))**(sympy.S(9)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(9)/2)) + sqrt(1 + 1/(a*x))*(1520*a - 325/x)*(c - c/(a*x))**(sympy.S(9)/2)/(7*a**2*(1 - 1/(a*x))**(sympy.S(9)/2)) + 135*sqrt(1 + 1/(a*x))*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(9)/2)/(7*a**3*(1 - 1/(a*x))**(sympy.S(9)/2)) + 65*sqrt(1 + 1/(a*x))*(a - 1/x)**3*(c - c/(a*x))**(sympy.S(9)/2)/(7*a**4*(1 - 1/(a*x))**(sympy.S(9)/2)) + x*(a - 1/x)**5*(c - c/(a*x))**(sympy.S(9)/2)/(a**5*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x))) + 10*(a - 1/x)**4*(c - c/(a*x))**(sympy.S(9)/2)/(a**5*(1 - 1/(a*x))**(sympy.S(9)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_467():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(-3*acoth(a*x))
    F = -13*(c - c/(a*x))**(sympy.S(7)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(7)/2)) + sqrt(1 + 1/(a*x))*(1360*a - 311/x)*(c - c/(a*x))**(sympy.S(7)/2)/(15*a**2*(1 - 1/(a*x))**(sympy.S(7)/2)) + 47*sqrt(1 + 1/(a*x))*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(7)/2)/(5*a**3*(1 - 1/(a*x))**(sympy.S(7)/2)) + x*(a - 1/x)**4*(c - c/(a*x))**(sympy.S(7)/2)/(a**4*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))) + 10*(a - 1/x)**3*(c - c/(a*x))**(sympy.S(7)/2)/(a**4*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_468():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(-3*acoth(a*x))
    F = -11*(c - c/(a*x))**(sympy.S(5)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(5)/2)) + sqrt(1 + 1/(a*x))*(112*a - 29/x)*(c - c/(a*x))**(sympy.S(5)/2)/(3*a**2*(1 - 1/(a*x))**(sympy.S(5)/2)) + x*(a - 1/x)**3*(c - c/(a*x))**(sympy.S(5)/2)/(a**3*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))) + 10*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(5)/2)/(a**3*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_469():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(-3*acoth(a*x))
    F = -9*(c - c/(a*x))**(sympy.S(3)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(1 - 1/(a*x))**(sympy.S(3)/2)) + x*(a - 1/x)**2*(c - c/(a*x))**(sympy.S(3)/2)/(a**2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))) + (21*a + 1/x)*(c - c/(a*x))**(sympy.S(3)/2)/(a**2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_470():
    f = sqrt(c - c/(a*x))*exp(-3*acoth(a*x))
    F = x*sqrt(c - c/(a*x))/(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 7*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(a*sqrt(1 - 1/(a*x))) + 9*sqrt(c - c/(a*x))/(a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_471():
    f = exp(-3*acoth(a*x))/sqrt(c - c/(a*x))
    F = x*sqrt(c - c/(a*x))/(c*sqrt(1 - 1/(a**2*x**2))) + 5*sqrt(c - c/(a*x))/(a*c*sqrt(1 - 1/(a**2*x**2))) - 5*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_472():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = 3*x*sqrt(1 - 1/(a**2*x**2))/(c*sqrt(c - c/(a*x))) - 2*x*sqrt(c - c/(a*x))/(c**2*sqrt(1 - 1/(a**2*x**2))) - 3*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_473():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = x*(1 - 1/(a*x))**(sympy.S(5)/2)/(sqrt(1 + 1/(a*x))*(c - c/(a*x))**(sympy.S(5)/2)) - sqrt(2)*(1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(2*a*(c - c/(a*x))**(sympy.S(5)/2)) - (1 - 1/(a*x))**(sympy.S(5)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(5)/2)) + 2*(1 - 1/(a*x))**(sympy.S(5)/2)/(a*sqrt(1 + 1/(a*x))*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_474():
    f = exp(-3*acoth(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = a*x*(1 - 1/(a*x))**(sympy.S(7)/2)/(sqrt(1 + 1/(a*x))*(a - 1/x)*(c - c/(a*x))**(sympy.S(7)/2)) - 3*(1 - 1/(a*x))**(sympy.S(7)/2)/(sqrt(1 + 1/(a*x))*(2*a - 2/x)*(c - c/(a*x))**(sympy.S(7)/2)) - 11*sqrt(2)*(1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(8*a*(c - c/(a*x))**(sympy.S(7)/2)) + (1 - 1/(a*x))**(sympy.S(7)/2)*atanh(sqrt(1 + 1/(a*x)))/(a*(c - c/(a*x))**(sympy.S(7)/2)) + 7*(1 - 1/(a*x))**(sympy.S(7)/2)/(4*a*sqrt(1 + 1/(a*x))*(c - c/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_475():
    f = x**m*sqrt(c - c/(a*x))*exp(acoth(a*x))
    F = x**(m + 1)*sqrt(c - c/(a*x))*hyper((sympy.S(-1)/2, -m - 1), (-m,), -1/(a*x))/(sqrt(1 - 1/(a*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_476():
    f = x**2*sqrt(c - c/(a*x))*exp(acoth(a*x))
    F = c*x**3*sqrt(1 - 1/(a**2*x**2))/(3*sqrt(c - c/(a*x))) + c*x**2*sqrt(1 - 1/(a**2*x**2))/(12*a*sqrt(c - c/(a*x))) - c*x*sqrt(1 - 1/(a**2*x**2))/(8*a**2*sqrt(c - c/(a*x))) + sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_477():
    f = x*sqrt(c - c/(a*x))*exp(acoth(a*x))
    F = c*x**2*sqrt(1 - 1/(a**2*x**2))/(2*sqrt(c - c/(a*x))) + c*x*sqrt(1 - 1/(a**2*x**2))/(4*a*sqrt(c - c/(a*x))) - sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_478():
    f = sqrt(c - c/(a*x))*exp(acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)) + sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_479():
    f = sqrt(c - c/(a*x))*exp(acoth(a*x))/x
    F = 2*sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x))) - 2*c*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_480():
    f = sqrt(c - c/(a*x))*exp(acoth(a*x))/x**2
    F = -2*a*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(3*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_481():
    f = sqrt(c - c/(a*x))*exp(acoth(a*x))/x**3
    F = -2*a**2*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(15*(c - c/(a*x))**(sympy.S(3)/2)) + 2*a**2*c*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(5*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_482():
    f = sqrt(c - c/(a*x))*exp(acoth(a*x))/x**4
    F = 8*a**3*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(105*(c - c/(a*x))**(sympy.S(3)/2)) - 8*a**3*c*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(35*sqrt(c - c/(a*x))) - 2*a*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(7*x**2*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_483():
    f = sqrt(c - c/(a*x))*exp(acoth(a*x))/x**5
    F = -16*a**4*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(315*(c - c/(a*x))**(sympy.S(3)/2)) + 16*a**4*c*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(105*sqrt(c - c/(a*x))) + 4*a**2*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(21*x**2*(c - c/(a*x))**(sympy.S(3)/2)) - 2*a*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(9*x**3*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_484():
    f = x**3*sqrt(c - c/(a*x))*exp(2*acoth(a*x))
    F = x**4*sqrt(c - c/(a*x))/4 + 5*x**3*sqrt(c - c/(a*x))/(8*a) + 25*x**2*sqrt(c - c/(a*x))/(32*a**2) + 75*x*sqrt(c - c/(a*x))/(64*a**3) + 75*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_485():
    f = x**2*sqrt(c - c/(a*x))*exp(2*acoth(a*x))
    F = x**3*sqrt(c - c/(a*x))/3 + 11*x**2*sqrt(c - c/(a*x))/(12*a) + 11*x*sqrt(c - c/(a*x))/(8*a**2) + 11*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_486():
    f = x*sqrt(c - c/(a*x))*exp(2*acoth(a*x))
    F = x**2*sqrt(c - c/(a*x))/2 + 7*x*sqrt(c - c/(a*x))/(4*a) + 7*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_487():
    f = sqrt(c - c/(a*x))*exp(2*acoth(a*x))
    F = x*sqrt(c - c/(a*x)) + 3*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_488():
    f = sqrt(c - c/(a*x))*exp(2*acoth(a*x))/x
    F = 2*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c)) + 2*sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_489():
    f = sqrt(c - c/(a*x))*exp(2*acoth(a*x))/x**2
    F = 4*a*sqrt(c - c/(a*x)) - 2*a*(c - c/(a*x))**(sympy.S(3)/2)/(3*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_490():
    f = sqrt(c - c/(a*x))*exp(2*acoth(a*x))/x**3
    F = 4*a**2*sqrt(c - c/(a*x)) - 2*a**2*(c - c/(a*x))**(sympy.S(3)/2)/c + 2*a**2*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_491():
    f = sqrt(c - c/(a*x))*exp(2*acoth(a*x))/x**4
    F = 4*a**3*sqrt(c - c/(a*x)) - 10*a**3*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) + 8*a**3*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2) - 2*a**3*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_492():
    f = sqrt(c - c/(a*x))*exp(2*acoth(a*x))/x**5
    F = 4*a**4*sqrt(c - c/(a*x)) - 14*a**4*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) + 18*a**4*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2) - 10*a**4*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3) + 2*a**4*(c - c/(a*x))**(sympy.S(9)/2)/(9*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_493():
    f = x**3*sqrt(c - c/(a*x))*exp(3*acoth(a*x))
    F = x**4*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(4*sqrt(1 - 1/(a*x))) + 17*x**3*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(24*a*sqrt(1 - 1/(a*x))) + 107*x**2*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(96*a**2*sqrt(1 - 1/(a*x))) + 149*x*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(64*a**3*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a**4*sqrt(1 - 1/(a*x))) + 363*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(64*a**4*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_494():
    f = x**2*sqrt(c - c/(a*x))*exp(3*acoth(a*x))
    F = x**3*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(3*sqrt(1 - 1/(a*x))) + 13*x**2*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(12*a*sqrt(1 - 1/(a*x))) + 19*x*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(8*a**2*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a**3*sqrt(1 - 1/(a*x))) + 45*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(8*a**3*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_495():
    f = x*sqrt(c - c/(a*x))*exp(3*acoth(a*x))
    F = x**2*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(2*sqrt(1 - 1/(a*x))) + 9*x*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(4*a*sqrt(1 - 1/(a*x))) - 4*sqrt(2)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a**2*sqrt(1 - 1/(a*x))) + 23*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(4*a**2*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_496():
    f = sqrt(c - c/(a*x))*exp(3*acoth(a*x))
    F = x*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/sqrt(1 - 1/(a*x)) - 4*sqrt(2)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/(a*sqrt(1 - 1/(a*x))) + 5*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(a*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_497():
    f = sqrt(c - c/(a*x))*exp(3*acoth(a*x))/x
    F = 2*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/sqrt(1 - 1/(a*x)) - 4*sqrt(2)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/sqrt(1 - 1/(a*x)) + 2*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/sqrt(1 - 1/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_498():
    f = sqrt(c - c/(a*x))*exp(3*acoth(a*x))/x**2
    F = -4*sqrt(2)*a*sqrt(c)*atanh(sqrt(2)*sqrt(c)*sqrt(1 - 1/(a**2*x**2))/(2*sqrt(c - c/(a*x)))) + 2*a*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(3*(c - c/(a*x))**(sympy.S(3)/2)) + 4*a*c*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_499():
    f = sqrt(c - c/(a*x))*exp(3*acoth(a*x))/x**3
    F = -4*sqrt(2)*a**2*sqrt(c)*atanh(sqrt(2)*sqrt(c)*sqrt(1 - 1/(a**2*x**2))/(2*sqrt(c - c/(a*x)))) + 2*a**2*c**3*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(5*(c - c/(a*x))**(sympy.S(5)/2)) + 2*a**2*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(3*(c - c/(a*x))**(sympy.S(3)/2)) + 4*a**2*c*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_500():
    f = sqrt(c - c/(a*x))*exp(3*acoth(a*x))/x**4
    F = -4*sqrt(2)*a**3*sqrt(c)*atanh(sqrt(2)*sqrt(c)*sqrt(1 - 1/(a**2*x**2))/(2*sqrt(c - c/(a*x)))) + 4*a**3*c**3*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(7*(c - c/(a*x))**(sympy.S(5)/2)) - 2*a**3*c**2*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(7*(c - c/(a*x))**(sympy.S(3)/2)) + 2*a**3*c**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(3*(c - c/(a*x))**(sympy.S(3)/2)) + 4*a**3*c*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_501():
    f = sqrt(c - c/(a*x))*exp(3*acoth(a*x))/x**5
    F = 2*a**4*(1 + 1/(a*x))**(sympy.S(9)/2)*sqrt(c - c/(a*x))/(9*sqrt(1 - 1/(a*x))) - 2*a**4*(1 + 1/(a*x))**(sympy.S(7)/2)*sqrt(c - c/(a*x))/(7*sqrt(1 - 1/(a*x))) + 2*a**4*(1 + 1/(a*x))**(sympy.S(5)/2)*sqrt(c - c/(a*x))/(5*sqrt(1 - 1/(a*x))) + 2*a**4*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(c - c/(a*x))/(3*sqrt(1 - 1/(a*x))) + 4*a**4*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/sqrt(1 - 1/(a*x)) - 4*sqrt(2)*a**4*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(1 + 1/(a*x))/2)/sqrt(1 - 1/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_502():
    f = x**m*sqrt(c - c/(a*x))*exp(-acoth(a*x))
    F = x**(m + 1)*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/(sqrt(1 - 1/(a*x))*(m + 1)) - x**m*sqrt(c - c/(a*x))*(4*m + 3)*hyper((sympy.S.Half, -m), (1 - m,), -1/(a*x))/(2*a*m*sqrt(1 - 1/(a*x))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_503():
    f = x**2*sqrt(c - c/(a*x))*exp(-acoth(a*x))
    F = c*x**3*sqrt(1 - 1/(a**2*x**2))/(3*sqrt(c - c/(a*x))) - 11*c*x**2*sqrt(1 - 1/(a**2*x**2))/(12*a*sqrt(c - c/(a*x))) + 11*c*x*sqrt(1 - 1/(a**2*x**2))/(8*a**2*sqrt(c - c/(a*x))) - 11*sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_504():
    f = x*sqrt(c - c/(a*x))*exp(-acoth(a*x))
    F = c*x**2*sqrt(1 - 1/(a**2*x**2))/(2*sqrt(c - c/(a*x))) - 7*c*x*sqrt(1 - 1/(a**2*x**2))/(4*a*sqrt(c - c/(a*x))) + 7*sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_505():
    f = sqrt(c - c/(a*x))*exp(-acoth(a*x))
    F = c*x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)) - 3*sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_506():
    f = sqrt(c - c/(a*x))*exp(-acoth(a*x))/x
    F = 2*sqrt(c)*atanh(sqrt(c)*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x))) + 2*c*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_507():
    f = sqrt(c - c/(a*x))*exp(-acoth(a*x))/x**2
    F = -8*a*c*sqrt(1 - 1/(a**2*x**2))/(3*sqrt(c - c/(a*x))) - 2*a*sqrt(1 - 1/(a**2*x**2))*sqrt(c - c/(a*x))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_508():
    f = sqrt(c - c/(a*x))*exp(-acoth(a*x))/x**3
    F = 8*a**2*c*sqrt(1 - 1/(a**2*x**2))/(5*sqrt(c - c/(a*x))) + 2*a**2*sqrt(1 - 1/(a**2*x**2))*sqrt(c - c/(a*x))/5 + 2*a**2*sqrt(1 - 1/(a**2*x**2))*(c - c/(a*x))**(sympy.S(3)/2)/(5*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_509():
    f = sqrt(c - c/(a*x))*exp(-acoth(a*x))/x**4
    F = -104*a**3*c*sqrt(1 - 1/(a**2*x**2))/(105*sqrt(c - c/(a*x))) - 104*a**3*sqrt(1 - 1/(a**2*x**2))*sqrt(c - c/(a*x))/105 - 26*a*c*sqrt(1 - 1/(a**2*x**2))/(35*x**2*sqrt(c - c/(a*x))) + 2*c*sqrt(1 - 1/(a**2*x**2))/(7*x**3*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_510():
    f = x**3*sqrt(c - c/(a*x))*exp(-2*acoth(a*x))
    F = x**4*sqrt(c - c/(a*x))/4 - 17*x**3*sqrt(c - c/(a*x))/(24*a) + 107*x**2*sqrt(c - c/(a*x))/(96*a**2) - 149*x*sqrt(c - c/(a*x))/(64*a**3) + 363*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(64*a**4) - 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_511():
    f = x**2*sqrt(c - c/(a*x))*exp(-2*acoth(a*x))
    F = x**3*sqrt(c - c/(a*x))/3 - 13*x**2*sqrt(c - c/(a*x))/(12*a) + 19*x*sqrt(c - c/(a*x))/(8*a**2) - 45*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(8*a**3) + 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_512():
    f = x*sqrt(c - c/(a*x))*exp(-2*acoth(a*x))
    F = x**2*sqrt(c - c/(a*x))/2 - 9*x*sqrt(c - c/(a*x))/(4*a) + 23*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(4*a**2) - 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_513():
    f = sqrt(c - c/(a*x))*exp(-2*acoth(a*x))
    F = x*sqrt(c - c/(a*x)) - 5*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_514():
    f = sqrt(c - c/(a*x))*exp(-2*acoth(a*x))/x
    F = 2*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c)) - 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) + 2*sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_515():
    f = sqrt(c - c/(a*x))*exp(-2*acoth(a*x))/x**2
    F = 4*sqrt(2)*a*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) - 4*a*sqrt(c - c/(a*x)) - 2*a*(c - c/(a*x))**(sympy.S(3)/2)/(3*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_516():
    f = sqrt(c - c/(a*x))*exp(-2*acoth(a*x))/x**3
    F = -4*sqrt(2)*a**2*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) + 4*a**2*sqrt(c - c/(a*x)) + 2*a**2*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) + 2*a**2*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_517():
    f = sqrt(c - c/(a*x))*exp(-2*acoth(a*x))/x**4
    F = 4*sqrt(2)*a**3*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) - 4*a**3*sqrt(c - c/(a*x)) - 2*a**3*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) - 2*a**3*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_518():
    f = sqrt(c - c/(a*x))*exp(-2*acoth(a*x))/x**5
    F = -4*sqrt(2)*a**4*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) + 4*a**4*sqrt(c - c/(a*x)) + 2*a**4*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) + 2*a**4*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2) - 2*a**4*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3) + 2*a**4*(c - c/(a*x))**(sympy.S(9)/2)/(9*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_519():
    f = x**3*sqrt(c - c/(a*x))*exp(-3*acoth(a*x))
    F = x**4*sqrt(c - c/(a*x))/(4*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 25*x**3*sqrt(c - c/(a*x))/(24*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 223*x**2*sqrt(c - c/(a*x))/(96*a**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 1115*x*sqrt(c - c/(a*x))/(192*a**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 1115*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(64*a**4*sqrt(1 - 1/(a*x))) - 1115*sqrt(c - c/(a*x))/(64*a**4*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_520():
    f = x**2*sqrt(c - c/(a*x))*exp(-3*acoth(a*x))
    F = x**3*sqrt(c - c/(a*x))/(3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 19*x**2*sqrt(c - c/(a*x))/(12*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 119*x*sqrt(c - c/(a*x))/(24*a**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 119*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(8*a**3*sqrt(1 - 1/(a*x))) + 119*sqrt(c - c/(a*x))/(8*a**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_521():
    f = x*sqrt(c - c/(a*x))*exp(-3*acoth(a*x))
    F = x**2*sqrt(c - c/(a*x))/(2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 13*x*sqrt(c - c/(a*x))/(4*a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) + 47*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(4*a**2*sqrt(1 - 1/(a*x))) - 47*sqrt(c - c/(a*x))/(4*a**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_522():
    f = sqrt(c - c/(a*x))*exp(-3*acoth(a*x))
    F = x*sqrt(c - c/(a*x))/(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 7*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/(a*sqrt(1 - 1/(a*x))) + 9*sqrt(c - c/(a*x))/(a*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_523():
    f = sqrt(c - c/(a*x))*exp(-3*acoth(a*x))/x
    F = -2*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/sqrt(1 - 1/(a*x)) + 2*sqrt(c - c/(a*x))*atanh(sqrt(1 + 1/(a*x)))/sqrt(1 - 1/(a*x)) - 8*sqrt(c - c/(a*x))/(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_524():
    f = sqrt(c - c/(a*x))*exp(-3*acoth(a*x))/x**2
    F = 64*a*sqrt(c - c/(a*x))/(3*sqrt(1 - 1/(a**2*x**2))) - 16*a*(c - c/(a*x))**(sympy.S(3)/2)/(3*c*sqrt(1 - 1/(a**2*x**2))) - 2*a*(c - c/(a*x))**(sympy.S(5)/2)/(3*c**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_525():
    f = sqrt(c - c/(a*x))*exp(-3*acoth(a*x))/x**3
    F = -224*a**2*c*sqrt(1 - 1/(a**2*x**2))/(15*sqrt(c - c/(a*x))) - 56*a**2*sqrt(1 - 1/(a**2*x**2))*sqrt(c - c/(a*x))/15 - 7*a**2*sqrt(1 - 1/(a**2*x**2))*(c - c/(a*x))**(sympy.S(3)/2)/(5*c) - a**2*(c - c/(a*x))**(sympy.S(7)/2)/(c**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_526():
    f = sqrt(c - c/(a*x))*exp(-3*acoth(a*x))/x**4
    F = 1888*a**3*c*sqrt(1 - 1/(a**2*x**2))/(105*sqrt(c - c/(a*x))) + 472*a**3*sqrt(1 - 1/(a**2*x**2))*sqrt(c - c/(a*x))/105 + 59*a**3*sqrt(1 - 1/(a**2*x**2))*(c - c/(a*x))**(sympy.S(3)/2)/(35*c) + 2*a**3*sqrt(1 - 1/(a**2*x**2))*(c - c/(a*x))**(sympy.S(5)/2)/(7*c**2) + a**3*(c - c/(a*x))**(sympy.S(7)/2)/(c**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_527():
    f = sqrt(c - c/(a*x))*exp(-3*acoth(a*x))/x**5
    F = -2*a**4*(1 + 1/(a*x))**(sympy.S(9)/2)*sqrt(c - c/(a*x))/(9*sqrt(1 - 1/(a*x))) + 2*a**4*(1 + 1/(a*x))**(sympy.S(7)/2)*sqrt(c - c/(a*x))/sqrt(1 - 1/(a*x)) - 38*a**4*(1 + 1/(a*x))**(sympy.S(5)/2)*sqrt(c - c/(a*x))/(5*sqrt(1 - 1/(a*x))) + 50*a**4*(1 + 1/(a*x))**(sympy.S(3)/2)*sqrt(c - c/(a*x))/(3*sqrt(1 - 1/(a*x))) - 32*a**4*sqrt(1 + 1/(a*x))*sqrt(c - c/(a*x))/sqrt(1 - 1/(a*x)) - 8*a**4*sqrt(c - c/(a*x))/(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_528():
    f = (c - c/(a*x))*exp(n*acoth(a*x))
    F = -2**(n/2)*c*(1 - 1/(a*x))**(1 - n/2)*hyper((1 - n/2, 1 - n/2), (2 - n/2,), (a - 1/x)/(2*a))/(a*(2 - n)) + c*x*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2) - 2*c*(1 - n)*(1 + 1/(a*x))**(n/2)*hyper((1, n/2), (n/2 + 1,), (a + 1/x)/(a - 1/x))/(a*n*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_529():
    f = exp(n*acoth(a*x))/(c - c/(a*x))
    F = x*(1 + 1/(a*x))**(n/2 + 1)/(c*(1 - 1/(a*x))**(n/2)) - (1 + 1/(a*x))**(n/2)*(2*n + 2)*hyper((1, -n/2), (1 - n/2,), (a - 1/x)/(a + 1/x))/(a*c*n*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_530():
    f = exp(n*acoth(a*x))/(c - c/(a*x))**2
    F = x*(1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 + 1)/c**2 - (1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 + 1)*(n + 3)/(a*c**2*(n + 2)) - (1 + 1/(a*x))**(n/2)*(2*n + 4)*hyper((1, -n/2), (1 - n/2,), (a - 1/x)/(a + 1/x))/(a*c**2*n*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_531():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(n*acoth(a*x))
    F = -2**(sympy.S(5)/2 - n/2)*(1 + 1/(a*x))**(n/2 + 1)*(c - c/(a*x))**(sympy.S(3)/2)*appellf1(n/2 + 1, 2, n/2 + sympy.S(-3)/2, n/2 + 2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*(1 - 1/(a*x))**(sympy.S(3)/2)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_532():
    f = sqrt(c - c/(a*x))*exp(n*acoth(a*x))
    F = -2**(sympy.S(3)/2 - n/2)*(1 + 1/(a*x))**(n/2 + 1)*sqrt(c - c/(a*x))*appellf1(n/2 + 1, 2, n/2 + sympy.S(-1)/2, n/2 + 2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*sqrt(1 - 1/(a*x))*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_533():
    f = exp(n*acoth(a*x))/sqrt(c - c/(a*x))
    F = -2**(sympy.S.Half - n/2)*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(n/2 + 1)*appellf1(n/2 + 1, 2, n/2 + sympy.S.Half, n/2 + 2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*sqrt(c - c/(a*x))*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_534():
    f = exp(n*acoth(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = -2**(-n/2 + sympy.S(-1)/2)*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(n/2 + 1)*appellf1(n/2 + 1, 2, n/2 + sympy.S(3)/2, n/2 + 2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*(c - c/(a*x))**(sympy.S(3)/2)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_535():
    f = (c - c/(a*x))**p*exp(n*acoth(a*x))
    F = -2**(-n/2 + p + 1)*(1 + 1/(a*x))**(n/2 + 1)*(c - c/(a*x))**p*appellf1(n/2 + 1, 2, n/2 - p, n/2 + 2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*(1 - 1/(a*x))**p*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_536():
    f = (c - c/(a*x))**p*exp(2*p*acoth(a*x))
    F = -(1 + 1/(a*x))**(p + 1)*(c - c/(a*x))**p*hyper((2, p + 1), (p + 2,), 1 + 1/(a*x))/(a*(1 - 1/(a*x))**p*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_537():
    f = (c - c/(a*x))**p*exp(-2*p*acoth(a*x))
    F = -4**p*(1 + 1/(a*x))**(1 - p)*(c - c/(a*x))**p*appellf1(1 - p, 2, -2*p, 2 - p, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*(1 - p)*(1 - 1/(a*x))**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_538():
    f = (c - c/(a*x))**p*exp(2*acoth(a*x))
    F = x*(c - c/(a*x))**p + (2 - p)*(c - c/(a*x))**p*hyper((1, p), (p + 1,), 1 - 1/(a*x))/(a*p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_539():
    f = (c - c/(a*x))**p*exp(acoth(a*x))
    F = -2**(p + sympy.S.Half)*(1 + 1/(a*x))**(sympy.S(3)/2)*(c - c/(a*x))**p*appellf1(sympy.S(3)/2, 2, sympy.S.Half - p, sympy.S(5)/2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(3*a*(1 - 1/(a*x))**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_540():
    f = (c - c/(a*x))**p*exp(-acoth(a*x))
    F = -2**(p + sympy.S(3)/2)*sqrt(1 + 1/(a*x))*(c - c/(a*x))**p*appellf1(sympy.S.Half, 2, -p + sympy.S(-1)/2, sympy.S(3)/2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*(1 - 1/(a*x))**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_541():
    f = (c - c/(a*x))**p*exp(-2*acoth(a*x))
    F = x*(c - c/(a*x))**(p + 2)/c**2 - (c - c/(a*x))**(p + 2)*hyper((1, p + 2), (p + 3,), 1 - 1/(a*x))/(a*c**2) + (c - c/(a*x))**(p + 2)*hyper((1, p + 2), (p + 3,), (a - 1/x)/(2*a))/(2*a*c**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_542():
    f = (-a**2*c*x**2 + c)**4*exp(acoth(a*x))
    F = a**8*c**4*x**9*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(11)/2)/9 - 7*a**7*c**4*x**8*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(11)/2)/72 + 5*a**6*c**4*x**7*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(11)/2)/72 - 5*a**5*c**4*x**6*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(11)/2)/144 + a**4*c**4*x**5*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/144 + a**3*c**4*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/64 + 7*a**2*c**4*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/192 + 35*a*c**4*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/384 + 35*c**4*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/128 + 35*c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_543():
    f = (-a**2*c*x**2 + c)**3*exp(acoth(a*x))
    F = -a**6*c**3*x**7*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/7 + 5*a**5*c**3*x**6*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/42 - a**4*c**3*x**5*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/14 + a**3*c**3*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/56 + a**2*c**3*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/24 + 5*a*c**3*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/48 + 5*c**3*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/16 + 5*c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_544():
    f = (-a**2*c*x**2 + c)**2*exp(acoth(a*x))
    F = a**4*c**2*x**5*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/5 - 3*a**3*c**2*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/20 + a**2*c**2*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/20 + a*c**2*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/8 + 3*c**2*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/8 + 3*c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_545():
    f = (-a**2*c*x**2 + c)*exp(acoth(a*x))
    F = -a**2*c*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/3 + a*c*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/6 + c*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/2 + c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_546():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)
    F = exp(acoth(a*x))/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_547():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = -(-2*a*x + 1)*exp(acoth(a*x))/(3*a*c**2*(-a**2*x**2 + 1)) + 2*exp(acoth(a*x))/(3*a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_548():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = -(-4*a*x + 1)*exp(acoth(a*x))/(15*a*c**3*(-a**2*x**2 + 1)**2) - 4*(-2*a*x + 1)*exp(acoth(a*x))/(15*a*c**3*(-a**2*x**2 + 1)) + 8*exp(acoth(a*x))/(15*a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_549():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = -(-6*a*x + 1)*exp(acoth(a*x))/(35*a*c**4*(-a**2*x**2 + 1)**3) - 2*(-4*a*x + 1)*exp(acoth(a*x))/(35*a*c**4*(-a**2*x**2 + 1)**2) - 8*(-2*a*x + 1)*exp(acoth(a*x))/(35*a*c**4*(-a**2*x**2 + 1)) + 16*exp(acoth(a*x))/(35*a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_550():
    f = (-a**2*c*x**2 + c)**5*exp(2*acoth(a*x))
    F = -c**5*(a*x + 1)**11/(11*a) + 4*c**5*(a*x + 1)**10/(5*a) - 8*c**5*(a*x + 1)**9/(3*a) + 4*c**5*(a*x + 1)**8/a - 16*c**5*(a*x + 1)**7/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_551():
    f = (-a**2*c*x**2 + c)**4*exp(2*acoth(a*x))
    F = c**4*(a*x + 1)**9/(9*a) - 3*c**4*(a*x + 1)**8/(4*a) + 12*c**4*(a*x + 1)**7/(7*a) - 4*c**4*(a*x + 1)**6/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_552():
    f = (-a**2*c*x**2 + c)**3*exp(2*acoth(a*x))
    F = -c**3*(a*x + 1)**7/(7*a) + 2*c**3*(a*x + 1)**6/(3*a) - 4*c**3*(a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_553():
    f = (-a**2*c*x**2 + c)**2*exp(2*acoth(a*x))
    F = c**2*(a*x + 1)**5/(5*a) - c**2*(a*x + 1)**4/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_554():
    f = (-a**2*c*x**2 + c)*exp(2*acoth(a*x))
    F = -c*(a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_555():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)
    F = -1/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_556():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = -atanh(a*x)/(4*a*c**2) - 1/(4*a*c**2*(-a*x + 1)) - 1/(4*a*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_557():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = -atanh(a*x)/(4*a*c**3) + 1/(16*a*c**3*(a*x + 1)) - 3/(16*a*c**3*(-a*x + 1)) - 1/(8*a*c**3*(-a*x + 1)**2) - 1/(12*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_558():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = -15*atanh(a*x)/(64*a*c**4) + 5/(64*a*c**4*(a*x + 1)) + 1/(64*a*c**4*(a*x + 1)**2) - 5/(32*a*c**4*(-a*x + 1)) - 3/(32*a*c**4*(-a*x + 1)**2) - 1/(16*a*c**4*(-a*x + 1)**3) - 1/(32*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_559():
    f = (-a**2*c*x**2 + c)**4*exp(3*acoth(a*x))
    F = a**8*c**4*x**9*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(13)/2)/9 - 5*a**7*c**4*x**8*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(13)/2)/72 + 5*a**6*c**4*x**7*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(13)/2)/168 - 5*a**5*c**4*x**6*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(11)/2)/1008 - 11*a**4*c**4*x**5*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/1008 - 11*a**3*c**4*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/448 - 11*a**2*c**4*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/192 - 55*a*c**4*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/384 - 55*c**4*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/128 - 55*c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_560():
    f = (-a**2*c*x**2 + c)**3*exp(3*acoth(a*x))
    F = -a**6*c**3*x**7*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(11)/2)/7 + a**5*c**3*x**6*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(11)/2)/14 - a**4*c**3*x**5*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/70 - 9*a**3*c**3*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/280 - 3*a**2*c**3*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/40 - 3*a*c**3*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/16 - 9*c**3*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/16 - 9*c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_561():
    f = (-a**2*c*x**2 + c)**2*exp(3*acoth(a*x))
    F = a**4*c**2*x**5*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/5 - a**3*c**2*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/20 - 7*a**2*c**2*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/60 - 7*a*c**2*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/24 - 7*c**2*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/8 - 7*c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_562():
    f = (-a**2*c*x**2 + c)*exp(3*acoth(a*x))
    F = -a**2*c*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/3 - 5*a*c*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/6 - 5*c*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/2 - 5*c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_563():
    f = exp(3*acoth(a*x))/(-a**2*c*x**2 + c)
    F = exp(3*acoth(a*x))/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_564():
    f = exp(3*acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = (-2*a*x + 3)*exp(3*acoth(a*x))/(5*a*c**2*(-a**2*x**2 + 1)) - 2*exp(3*acoth(a*x))/(15*a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_565():
    f = exp(3*acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = -(-4*a*x + 3)*exp(3*acoth(a*x))/(7*a*c**3*(-a**2*x**2 + 1)**2) + 12*(-2*a*x + 3)*exp(3*acoth(a*x))/(35*a*c**3*(-a**2*x**2 + 1)) - 8*exp(3*acoth(a*x))/(35*a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_566():
    f = exp(3*acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = -10*(-4*a*x + 3)*exp(3*acoth(a*x))/(63*a*c**4*(-a**2*x**2 + 1)**2) - (-2*a*x + 1)*exp(3*acoth(a*x))/(9*a*c**4*(-a**2*x**2 + 1)**3) + 8*(-2*a*x + 3)*exp(3*acoth(a*x))/(21*a*c**4*(-a**2*x**2 + 1)) - 16*exp(3*acoth(a*x))/(63*a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_567():
    f = (-a**2*c*x**2 + c)**5*exp(4*acoth(a*x))
    F = -c**5*(a*x + 1)**11/(11*a) + 3*c**5*(a*x + 1)**10/(5*a) - 4*c**5*(a*x + 1)**9/(3*a) + c**5*(a*x + 1)**8/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_568():
    f = (-a**2*c*x**2 + c)**4*exp(4*acoth(a*x))
    F = c**4*(a*x + 1)**9/(9*a) - c**4*(a*x + 1)**8/(2*a) + 4*c**4*(a*x + 1)**7/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_569():
    f = (-a**2*c*x**2 + c)**3*exp(4*acoth(a*x))
    F = -c**3*(a*x + 1)**7/(7*a) + c**3*(a*x + 1)**6/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_570():
    f = (-a**2*c*x**2 + c)**2*exp(4*acoth(a*x))
    F = c**2*(a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_571():
    f = (-a**2*c*x**2 + c)*exp(4*acoth(a*x))
    F = -4*c*x - c*(a*x + 1)**3/(3*a) - c*(a*x + 1)**2/a - 8*c*log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_572():
    f = exp(4*acoth(a*x))/(-a**2*c*x**2 + c)
    F = x/(c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_573():
    f = exp(4*acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = 1/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_574():
    f = exp(4*acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = atanh(a*x)/(16*a*c**3) + 1/(16*a*c**3*(-a*x + 1)) + 1/(16*a*c**3*(-a*x + 1)**2) + 1/(12*a*c**3*(-a*x + 1)**3) + 1/(8*a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_575():
    f = exp(4*acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = 3*atanh(a*x)/(32*a*c**4) - 1/(64*a*c**4*(a*x + 1)) + 5/(64*a*c**4*(-a*x + 1)) + 1/(16*a*c**4*(-a*x + 1)**2) + 1/(16*a*c**4*(-a*x + 1)**3) + 1/(16*a*c**4*(-a*x + 1)**4) + 1/(20*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_576():
    f = (-a**2*c*x**2 + c)**4*exp(-acoth(a*x))
    F = a**8*c**4*x**9*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/9 - a**7*c**4*x**8*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/8 + a**6*c**4*x**7*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/8 - 5*a**5*c**4*x**6*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/48 + a**4*c**4*x**5*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/16 - a**3*c**4*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/64 - 7*a**2*c**4*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/192 - 35*a*c**4*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/384 - 35*c**4*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/128 - 35*c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_577():
    f = (-a**2*c*x**2 + c)**3*exp(-acoth(a*x))
    F = -a**6*c**3*x**7*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/7 + a**5*c**3*x**6*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/6 - a**4*c**3*x**5*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/6 + a**3*c**3*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/8 - a**2*c**3*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/24 - 5*a*c**3*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/48 - 5*c**3*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/16 - 5*c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_578():
    f = (-a**2*c*x**2 + c)**2*exp(-acoth(a*x))
    F = a**4*c**2*x**5*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/5 - a**3*c**2*x**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/4 + a**2*c**2*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/4 - a*c**2*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/8 - 3*c**2*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/8 - 3*c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_579():
    f = (-a**2*c*x**2 + c)*exp(-acoth(a*x))
    F = -a**2*c*x**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/3 + a*c*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/2 - c*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/2 - c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_580():
    f = exp(-acoth(a*x))/(-a**2*c*x**2 + c)
    F = -exp(-acoth(a*x))/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_581():
    f = exp(-acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = (2*a*x + 1)*exp(-acoth(a*x))/(3*a*c**2*(-a**2*x**2 + 1)) - 2*exp(-acoth(a*x))/(3*a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_582():
    f = exp(-acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = (4*a*x + 1)*exp(-acoth(a*x))/(15*a*c**3*(-a**2*x**2 + 1)**2) + (8*a*x + 4)*exp(-acoth(a*x))/(15*a*c**3*(-a**2*x**2 + 1)) - 8*exp(-acoth(a*x))/(15*a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_583():
    f = exp(-acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = (6*a*x + 1)*exp(-acoth(a*x))/(35*a*c**4*(-a**2*x**2 + 1)**3) + (8*a*x + 2)*exp(-acoth(a*x))/(35*a*c**4*(-a**2*x**2 + 1)**2) + (16*a*x + 8)*exp(-acoth(a*x))/(35*a*c**4*(-a**2*x**2 + 1)) - 16*exp(-acoth(a*x))/(35*a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_584():
    f = (-a**2*c*x**2 + c)**4*exp(-2*acoth(a*x))
    F = -c**4*(-a*x + 1)**9/(9*a) + 3*c**4*(-a*x + 1)**8/(4*a) - 12*c**4*(-a*x + 1)**7/(7*a) + 4*c**4*(-a*x + 1)**6/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_585():
    f = (-a**2*c*x**2 + c)**3*exp(-2*acoth(a*x))
    F = c**3*(-a*x + 1)**7/(7*a) - 2*c**3*(-a*x + 1)**6/(3*a) + 4*c**3*(-a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_586():
    f = (-a**2*c*x**2 + c)**2*exp(-2*acoth(a*x))
    F = -c**2*(-a*x + 1)**5/(5*a) + c**2*(-a*x + 1)**4/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_587():
    f = (-a**2*c*x**2 + c)*exp(-2*acoth(a*x))
    F = c*(-a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_588():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)
    F = 1/(a*c*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_589():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = -atanh(a*x)/(4*a*c**2) + 1/(4*a*c**2*(a*x + 1)) + 1/(4*a*c**2*(a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_590():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = -atanh(a*x)/(4*a*c**3) + 3/(16*a*c**3*(a*x + 1)) + 1/(8*a*c**3*(a*x + 1)**2) + 1/(12*a*c**3*(a*x + 1)**3) - 1/(16*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_591():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = -15*atanh(a*x)/(64*a*c**4) + 5/(32*a*c**4*(a*x + 1)) + 3/(32*a*c**4*(a*x + 1)**2) + 1/(16*a*c**4*(a*x + 1)**3) + 1/(32*a*c**4*(a*x + 1)**4) - 5/(64*a*c**4*(-a*x + 1)) - 1/(64*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_592():
    f = (-a**2*c*x**2 + c)**4*exp(-3*acoth(a*x))
    F = a**8*c**4*x**9*(1 - 1/(a*x))**(sympy.S(11)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/9 - 11*a**7*c**4*x**8*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/72 + 11*a**6*c**4*x**7*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/56 - 11*a**5*c**4*x**6*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/48 + 11*a**4*c**4*x**5*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/48 - 11*a**3*c**4*x**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/64 + 11*a**2*c**4*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/192 + 55*a*c**4*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/384 + 55*c**4*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/128 + 55*c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_593():
    f = (-a**2*c*x**2 + c)**3*exp(-3*acoth(a*x))
    F = -a**6*c**3*x**7*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/7 + 3*a**5*c**3*x**6*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/14 - 3*a**4*c**3*x**5*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/10 + 3*a**3*c**3*x**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/8 - 3*a**2*c**3*x**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/8 + 3*a*c**3*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/16 + 9*c**3*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/16 + 9*c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_594():
    f = (-a**2*c*x**2 + c)**2*exp(-3*acoth(a*x))
    F = a**4*c**2*x**5*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/5 - 7*a**3*c**2*x**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/20 + 7*a**2*c**2*x**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/12 - 7*a*c**2*x**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/8 + 7*c**2*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/8 + 7*c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_595():
    f = (-a**2*c*x**2 + c)*exp(-3*acoth(a*x))
    F = -a**2*c*x**3*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/3 + 5*a*c*x**2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/6 - 5*c*x*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/2 + 5*c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_596():
    f = exp(-3*acoth(a*x))/(-a**2*c*x**2 + c)
    F = -exp(-3*acoth(a*x))/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_597():
    f = exp(-3*acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = -(2*a*x + 3)*exp(-3*acoth(a*x))/(5*a*c**2*(-a**2*x**2 + 1)) + 2*exp(-3*acoth(a*x))/(15*a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_598():
    f = exp(-3*acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = (4*a*x + 3)*exp(-3*acoth(a*x))/(7*a*c**3*(-a**2*x**2 + 1)**2) - (24*a*x + 36)*exp(-3*acoth(a*x))/(35*a*c**3*(-a**2*x**2 + 1)) + 8*exp(-3*acoth(a*x))/(35*a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_599():
    f = exp(-3*acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = (2*a*x + 1)*exp(-3*acoth(a*x))/(9*a*c**4*(-a**2*x**2 + 1)**3) - (16*a*x + 24)*exp(-3*acoth(a*x))/(21*a*c**4*(-a**2*x**2 + 1)) + (40*a*x + 30)*exp(-3*acoth(a*x))/(63*a*c**4*(-a**2*x**2 + 1)**2) + 16*exp(-3*acoth(a*x))/(63*a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_600():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(acoth(a*x))
    F = (a*x + 1)**10*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(10*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 8*(a*x + 1)**9*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(9*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) + 3*(a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 32*(a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(7*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) + 8*(a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(3*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_601():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(acoth(a*x))
    F = (a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(8*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 6*(a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(7*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) + 2*(a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 8*(a*x + 1)**5*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(5*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_602():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(acoth(a*x))
    F = (a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) - 4*(a*x + 1)**5*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(5*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) + (a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_603():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(acoth(a*x))
    F = (a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - 2*(a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(3*a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_604():
    f = sqrt(-a**2*c*x**2 + c)*exp(acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_605():
    f = exp(acoth(a*x))/sqrt(-a**2*c*x**2 + c)
    F = x*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*c*x**2 + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_606():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*atanh(a*x)/(2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_607():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -3*a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-4*a*x + 4)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_608():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = 5*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)*atanh(a*x)/(16*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(32*(a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + 3*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(32*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(24*(-a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + 3*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((-16*a*x + 16)*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_609():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(2*acoth(a*x))
    F = -77*c**4*x*sqrt(-a**2*c*x**2 + c)/256 - 77*c**3*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/384 - 77*c**2*x*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/480 - 11*c*x*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/80 - 77*c**(sympy.S(9)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(256*a) + (a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(10*a) + 11*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(90*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_610():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(2*acoth(a*x))
    F = -45*c**3*x*sqrt(-a**2*c*x**2 + c)/128 - 15*c**2*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/64 - 3*c*x*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/16 - 45*c**(sympy.S(7)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(128*a) + (a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(8*a) + 9*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(56*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_611():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*acoth(a*x))
    F = -7*c**2*x*sqrt(-a**2*c*x**2 + c)/16 - 7*c*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/24 - 7*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(16*a) + (a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a) + 7*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(30*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_612():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*acoth(a*x))
    F = -5*c*x*sqrt(-a**2*c*x**2 + c)/8 - 5*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a) + (a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a) + 5*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(12*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_613():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))
    F = -3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(2*a) + (a*x + 1)*sqrt(-a**2*c*x**2 + c)/(2*a) + 3*sqrt(-a**2*c*x**2 + c)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_614():
    f = exp(2*acoth(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -(2*a*x + 2)/(a*sqrt(-a**2*c*x**2 + c)) + atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_615():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -x/(3*c*sqrt(-a**2*c*x**2 + c)) - (2*a*x + 2)/(3*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_616():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -x/(5*c*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 2*x/(5*c**2*sqrt(-a**2*c*x**2 + c)) - (2*a*x + 2)/(5*a*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_617():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -x/(7*c*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 4*x/(21*c**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 8*x/(21*c**3*sqrt(-a**2*c*x**2 + c)) - (2*a*x + 2)/(7*a*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_618():
    f = exp(2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(9)/2)
    F = -x/(9*c*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - 2*x/(15*c**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 8*x/(45*c**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 16*x/(45*c**4*sqrt(-a**2*c*x**2 + c)) - (2*a*x + 2)/(9*a*(-a**2*c*x**2 + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_619():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(3*acoth(a*x))
    F = (a*x + 1)**10*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(10*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 2*(a*x + 1)**9*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(3*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) + 3*(a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(2*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 8*(a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(7*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_620():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(3*acoth(a*x))
    F = (a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(8*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 4*(a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(7*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) + 2*(a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(3*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_621():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(3*acoth(a*x))
    F = (a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) - 2*(a*x + 1)**5*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(5*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_622():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(3*acoth(a*x))
    F = (a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_623():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**2*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_624():
    f = exp(3*acoth(a*x))/sqrt(-a**2*c*x**2 + c)
    F = x*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*c*x**2 + c) + 2*x*sqrt(1 - 1/(a**2*x**2))/((-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_625():
    f = exp(3*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_626():
    f = exp(3*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(6*(-a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_627():
    f = exp(3*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -5*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)*atanh(a*x)/(32*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((32*a*x + 32)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - 3*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(32*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(12*(-a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(16*(-a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((-8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_628():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(-acoth(a*x))
    F = (-a*x + 1)**10*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(10*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 8*(-a*x + 1)**9*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(9*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) + 3*(-a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 32*(-a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(7*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) + 8*(-a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(3*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_629():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(-acoth(a*x))
    F = (-a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(8*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 6*(-a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(7*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) + 2*(-a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 8*(-a*x + 1)**5*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(5*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_630():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(-acoth(a*x))
    F = (-a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) - 4*(-a*x + 1)**5*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(5*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) + (-a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_631():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-acoth(a*x))
    F = (-a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)) - 2*(-a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(3*a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_632():
    f = sqrt(-a**2*c*x**2 + c)*exp(-acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) - sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_633():
    f = exp(-acoth(a*x))/sqrt(-a**2*c*x**2 + c)
    F = x*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*c*x**2 + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_634():
    f = exp(-acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*atanh(a*x)/(2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_635():
    f = exp(-acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = 3*a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((4*a*x + 4)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_636():
    f = exp(-acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -5*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)*atanh(a*x)/(16*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + 3*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((16*a*x + 16)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + 3*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(32*(a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(24*(a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(32*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((-8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_637():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(-2*acoth(a*x))
    F = -7*c**2*x*sqrt(-a**2*c*x**2 + c)/16 - 7*c*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/24 - 7*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(16*a) - (-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a) - 7*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(30*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_638():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-2*acoth(a*x))
    F = -5*c*x*sqrt(-a**2*c*x**2 + c)/8 - 5*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a) - (-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a) - 5*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(12*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_639():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))
    F = -3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(2*a) - (-a*x + 1)*sqrt(-a**2*c*x**2 + c)/(2*a) - 3*sqrt(-a**2*c*x**2 + c)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_640():
    f = exp(-2*acoth(a*x))/sqrt(-a**2*c*x**2 + c)
    F = (-2*a*x + 2)/(a*sqrt(-a**2*c*x**2 + c)) + atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_641():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -x/(3*c*sqrt(-a**2*c*x**2 + c)) + (-2*a*x + 2)/(3*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_642():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -x/(5*c*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 2*x/(5*c**2*sqrt(-a**2*c*x**2 + c)) + (-2*a*x + 2)/(5*a*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_643():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -x/(7*c*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 4*x/(21*c**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 8*x/(21*c**3*sqrt(-a**2*c*x**2 + c)) + (-2*a*x + 2)/(7*a*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_644():
    f = exp(-2*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(9)/2)
    F = -x/(9*c*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - 2*x/(15*c**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 8*x/(45*c**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 16*x/(45*c**4*sqrt(-a**2*c*x**2 + c)) + (-2*a*x + 2)/(9*a*(-a**2*c*x**2 + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_645():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(-3*acoth(a*x))
    F = (-a*x + 1)**10*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(10*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 2*(-a*x + 1)**9*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(3*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) + 3*(-a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(2*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2)) - 8*(-a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(9)/2)/(7*a**10*x**9*(1 - 1/(a**2*x**2))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_646():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(-3*acoth(a*x))
    F = (-a*x + 1)**8*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(8*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) - 4*(-a*x + 1)**7*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(7*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)) + 2*(-a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(3*a**8*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_647():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(-3*acoth(a*x))
    F = (-a*x + 1)**6*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)) - 2*(-a*x + 1)**5*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(5*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_648():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-3*acoth(a*x))
    F = (-a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_649():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) - 3*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**2*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_650():
    f = exp(-3*acoth(a*x))/sqrt(-a**2*c*x**2 + c)
    F = x*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*c*x**2 + c) + 2*x*sqrt(1 - 1/(a**2*x**2))/((a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_651():
    f = exp(-3*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*(a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_652():
    f = exp(-3*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(6*(a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_653():
    f = exp(-3*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = 5*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)*atanh(a*x)/(32*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - 3*a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(32*(a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(12*(a*x + 1)**3*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/(16*(a*x + 1)**4*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) + a**6*x**7*(1 - 1/(a**2*x**2))**(sympy.S(7)/2)/((-32*a*x + 32)*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_654():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(acoth(a*x))
    F = x**3*sqrt(-a**2*c*x**2 + c)/(4*sqrt(1 - 1/(a**2*x**2))) + x**2*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_655():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(acoth(a*x))
    F = x**2*sqrt(-a**2*c*x**2 + c)/(3*sqrt(1 - 1/(a**2*x**2))) + x*sqrt(-a**2*c*x**2 + c)/(2*a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_656():
    f = sqrt(-a**2*c*x**2 + c)*exp(acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_657():
    f = sqrt(-a**2*c*x**2 + c)*exp(acoth(a*x))/x
    F = sqrt(-a**2*c*x**2 + c)/sqrt(1 - 1/(a**2*x**2)) + sqrt(-a**2*c*x**2 + c)*log(x)/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_658():
    f = sqrt(-a**2*c*x**2 + c)*exp(acoth(a*x))/x**2
    F = sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) - sqrt(-a**2*c*x**2 + c)/(a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_659():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))
    F = x**4*sqrt(-a**2*c*x**2 + c)/5 + x**3*sqrt(-a**2*c*x**2 + c)/(2*a) + 3*x**2*sqrt(-a**2*c*x**2 + c)/(5*a**2) - 3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(4*a**4) + (15*a*x + 24)*sqrt(-a**2*c*x**2 + c)/(20*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_660():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))
    F = x**3*sqrt(-a**2*c*x**2 + c)/4 + 2*x**2*sqrt(-a**2*c*x**2 + c)/(3*a) - 7*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a**3) + (21*a*x + 32)*sqrt(-a**2*c*x**2 + c)/(24*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_661():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))
    F = x**2*sqrt(-a**2*c*x**2 + c)/3 - sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/a**2 + (3*a*x + 5)*sqrt(-a**2*c*x**2 + c)/(3*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_662():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))
    F = -3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(2*a) + (a*x + 1)*sqrt(-a**2*c*x**2 + c)/(2*a) + 3*sqrt(-a**2*c*x**2 + c)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_663():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))/x
    F = -2*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) + sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + sqrt(-a**2*c*x**2 + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_664():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))/x**2
    F = -a*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) + 2*a*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + sqrt(-a**2*c*x**2 + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_665():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))/x**3
    F = 3*a**2*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/2 + 2*a*sqrt(-a**2*c*x**2 + c)/x + sqrt(-a**2*c*x**2 + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_666():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))/x**4
    F = a**3*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + 5*a**2*sqrt(-a**2*c*x**2 + c)/(3*x) + a*sqrt(-a**2*c*x**2 + c)/x**2 + sqrt(-a**2*c*x**2 + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_667():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))/x**5
    F = 7*a**4*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/8 + 4*a**3*sqrt(-a**2*c*x**2 + c)/(3*x) + 7*a**2*sqrt(-a**2*c*x**2 + c)/(8*x**2) + 2*a*sqrt(-a**2*c*x**2 + c)/(3*x**3) + sqrt(-a**2*c*x**2 + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_668():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))
    F = x**4*sqrt(-a**2*c*x**2 + c)/(5*sqrt(1 - 1/(a**2*x**2))) + 3*x**3*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(1 - 1/(a**2*x**2))) + 4*x**2*sqrt(-a**2*c*x**2 + c)/(3*a**2*sqrt(1 - 1/(a**2*x**2))) + 2*x*sqrt(-a**2*c*x**2 + c)/(a**3*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)/(a**4*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**5*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_669():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))
    F = x**3*sqrt(-a**2*c*x**2 + c)/(4*sqrt(1 - 1/(a**2*x**2))) + x**2*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))) + 2*x*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)/(a**3*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**4*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_670():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))
    F = x**2*sqrt(-a**2*c*x**2 + c)/(3*sqrt(1 - 1/(a**2*x**2))) + 3*x*sqrt(-a**2*c*x**2 + c)/(2*a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**3*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_671():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**2*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_672():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))/x
    F = sqrt(-a**2*c*x**2 + c)/sqrt(1 - 1/(a**2*x**2)) - sqrt(-a**2*c*x**2 + c)*log(x)/(a*x*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_673():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))/x**2
    F = -3*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_674():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))/x**3
    F = -4*a*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*a*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(2*a*x**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_675():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))/x**4
    F = -4*a**2*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*a**2*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*a*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(-a**2*c*x**2 + c)/(2*x**3*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(3*a*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_676():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))/x**5
    F = -4*a**3*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*a**3*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*a**2*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(1 - 1/(a**2*x**2))) + 2*a*sqrt(-a**2*c*x**2 + c)/(x**3*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(x**4*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(4*a*x**5*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_677():
    f = x**4*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x**5*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(a*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + 7*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(-a*x + 1)/(4*a**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(a*x + 1)/(4*a**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*a**2*(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_678():
    f = x**3*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*c*x**2 + c)**(sympy.S(3)/2) + 5*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(-a*x + 1)/(4*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(a*x + 1)/(4*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*a*(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_679():
    f = x**2*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = 3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(-a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_680():
    f = x*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -a*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*atanh(a*x)/(2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_681():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*atanh(a*x)/(2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**2*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_682():
    f = exp(acoth(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = a**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(x)/(-a**2*c*x**2 + c)**(sympy.S(3)/2) - 3*a**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(-a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - a**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_683():
    f = exp(acoth(a*x))/(x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(x)/(-a**2*c*x**2 + c)**(sympy.S(3)/2) - 5*a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(-a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - a**3*x**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_684():
    f = exp(acoth(a*x))/(x**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = 2*a**5*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(x)/(-a**2*c*x**2 + c)**(sympy.S(3)/2) - 7*a**5*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(-a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - a**5*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*log(a*x + 1)/(4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**5*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - a**4*x**2*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*c*x**2 + c)**(sympy.S(3)/2) - a**3*x*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)/(2*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_685():
    f = x**5*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = x**6*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*c*x**2 + c)**(sympy.S(5)/2) + 23*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(-a*x + 1)/(16*a*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 7*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(a*x + 1)/(16*a*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*a*(a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(a*(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*a*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_686():
    f = x**4*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = 11*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(-a*x + 1)/(16*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + 5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(a*x + 1)/(16*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + 3*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-4*a*x + 4)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_687():
    f = x**3*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -3*a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_688():
    f = x**2*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = a**2*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**2*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**2*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**2*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-4*a*x + 4)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_689():
    f = x*exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = a**3*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**3*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**3*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_690():
    f = exp(acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -3*a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*atanh(a*x)/(8*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**4*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-4*a*x + 4)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_691():
    f = exp(acoth(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = -a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(x)/(-a**2*c*x**2 + c)**(sympy.S(5)/2) + 11*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(-a*x + 1)/(16*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + 5*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(a*x + 1)/(16*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-2*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_692():
    f = exp(acoth(a*x))/(x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = -a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(x)/(-a**2*c*x**2 + c)**(sympy.S(5)/2) + 23*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(-a*x + 1)/(16*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 7*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*log(a*x + 1)/(16*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((8*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 3*a**6*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/((-4*a*x + 4)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + a**5*x**4*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_693():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(-acoth(a*x))
    F = x**3*sqrt(-a**2*c*x**2 + c)/(4*sqrt(1 - 1/(a**2*x**2))) - x**2*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_694():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(-acoth(a*x))
    F = x**2*sqrt(-a**2*c*x**2 + c)/(3*sqrt(1 - 1/(a**2*x**2))) - x*sqrt(-a**2*c*x**2 + c)/(2*a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_695():
    f = sqrt(-a**2*c*x**2 + c)*exp(-acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) - sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_696():
    f = sqrt(-a**2*c*x**2 + c)*exp(-acoth(a*x))/x
    F = sqrt(-a**2*c*x**2 + c)/sqrt(1 - 1/(a**2*x**2)) - sqrt(-a**2*c*x**2 + c)*log(x)/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_697():
    f = sqrt(-a**2*c*x**2 + c)*exp(-acoth(a*x))/x**2
    F = sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_698():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))
    F = x**4*sqrt(-a**2*c*x**2 + c)/5 - x**3*sqrt(-a**2*c*x**2 + c)/(2*a) + 3*x**2*sqrt(-a**2*c*x**2 + c)/(5*a**2) + 3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(4*a**4) + (-15*a*x + 24)*sqrt(-a**2*c*x**2 + c)/(20*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_699():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))
    F = x**3*sqrt(-a**2*c*x**2 + c)/4 - 2*x**2*sqrt(-a**2*c*x**2 + c)/(3*a) - 7*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a**3) - (-21*a*x + 32)*sqrt(-a**2*c*x**2 + c)/(24*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_700():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))
    F = x**2*sqrt(-a**2*c*x**2 + c)/3 + sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/a**2 + (-3*a*x + 5)*sqrt(-a**2*c*x**2 + c)/(3*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_701():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))
    F = -3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(2*a) - (-a*x + 1)*sqrt(-a**2*c*x**2 + c)/(2*a) - 3*sqrt(-a**2*c*x**2 + c)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_702():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))/x
    F = 2*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) + sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + sqrt(-a**2*c*x**2 + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_703():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))/x**2
    F = -a*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) - 2*a*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + sqrt(-a**2*c*x**2 + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_704():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))/x**3
    F = 3*a**2*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/2 - 2*a*sqrt(-a**2*c*x**2 + c)/x + sqrt(-a**2*c*x**2 + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_705():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))/x**4
    F = -a**3*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + 5*a**2*sqrt(-a**2*c*x**2 + c)/(3*x) - a*sqrt(-a**2*c*x**2 + c)/x**2 + sqrt(-a**2*c*x**2 + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_706():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))/x**5
    F = 7*a**4*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/8 - 4*a**3*sqrt(-a**2*c*x**2 + c)/(3*x) + 7*a**2*sqrt(-a**2*c*x**2 + c)/(8*x**2) - 2*a*sqrt(-a**2*c*x**2 + c)/(3*x**3) + sqrt(-a**2*c*x**2 + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_707():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))
    F = x**4*sqrt(-a**2*c*x**2 + c)/(5*sqrt(1 - 1/(a**2*x**2))) - 3*x**3*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(1 - 1/(a**2*x**2))) + 4*x**2*sqrt(-a**2*c*x**2 + c)/(3*a**2*sqrt(1 - 1/(a**2*x**2))) - 2*x*sqrt(-a**2*c*x**2 + c)/(a**3*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)/(a**4*sqrt(1 - 1/(a**2*x**2))) - 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**5*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_708():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))
    F = x**3*sqrt(-a**2*c*x**2 + c)/(4*sqrt(1 - 1/(a**2*x**2))) - x**2*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))) + 2*x*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(1 - 1/(a**2*x**2))) - 4*sqrt(-a**2*c*x**2 + c)/(a**3*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**4*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_709():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))
    F = x**2*sqrt(-a**2*c*x**2 + c)/(3*sqrt(1 - 1/(a**2*x**2))) - 3*x*sqrt(-a**2*c*x**2 + c)/(2*a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(1 - 1/(a**2*x**2))) - 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**3*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_710():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))
    F = x*sqrt(-a**2*c*x**2 + c)/(2*sqrt(1 - 1/(a**2*x**2))) - 3*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**2*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_711():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))/x
    F = sqrt(-a**2*c*x**2 + c)/sqrt(1 - 1/(a**2*x**2)) + sqrt(-a**2*c*x**2 + c)*log(x)/(a*x*sqrt(1 - 1/(a**2*x**2))) - 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_712():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))/x**2
    F = -3*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) - sqrt(-a**2*c*x**2 + c)/(a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_713():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))/x**3
    F = 4*a*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) - 4*a*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(1 - 1/(a**2*x**2))) - sqrt(-a**2*c*x**2 + c)/(2*a*x**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_714():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))/x**4
    F = -4*a**2*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*a**2*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) - 4*a*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(-a**2*c*x**2 + c)/(2*x**3*sqrt(1 - 1/(a**2*x**2))) - sqrt(-a**2*c*x**2 + c)/(3*a*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_715():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))/x**5
    F = 4*a**3*sqrt(-a**2*c*x**2 + c)*log(x)/(x*sqrt(1 - 1/(a**2*x**2))) - 4*a**3*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(x*sqrt(1 - 1/(a**2*x**2))) + 4*a**2*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(1 - 1/(a**2*x**2))) - 2*a*sqrt(-a**2*c*x**2 + c)/(x**3*sqrt(1 - 1/(a**2*x**2))) + sqrt(-a**2*c*x**2 + c)/(x**4*sqrt(1 - 1/(a**2*x**2))) - sqrt(-a**2*c*x**2 + c)/(4*a*x**5*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_716():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(3*acoth(a*x))
    F = x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(sqrt(1 - 1/(a**2*x**2))*(m + 2)) - 4*x**m*sqrt(-a**2*c*x**2 + c)*hyper((1, m + 1), (m + 2,), a*x)/(a*sqrt(1 - 1/(a**2*x**2))*(m + 1)) + 3*x**m*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_717():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(2*acoth(a*x))
    F = -2*a*c*x**(m + 2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*sqrt(-a**2*c*x**2 + c)) - c*x**(m + 1)*(2*m + 3)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(m + 2)*sqrt(-a**2*c*x**2 + c)) + x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(m + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_718():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(acoth(a*x))
    F = x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(sqrt(1 - 1/(a**2*x**2))*(m + 2)) + x**m*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_719():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(-acoth(a*x))
    F = x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(sqrt(1 - 1/(a**2*x**2))*(m + 2)) - x**m*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_720():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(-2*acoth(a*x))
    F = 2*a*c*x**(m + 2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*sqrt(-a**2*c*x**2 + c)) - c*x**(m + 1)*(2*m + 3)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(m + 2)*sqrt(-a**2*c*x**2 + c)) + x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(m + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_721():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(-3*acoth(a*x))
    F = x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(sqrt(1 - 1/(a**2*x**2))*(m + 2)) + 4*x**m*sqrt(-a**2*c*x**2 + c)*hyper((1, m + 1), (m + 2,), -a*x)/(a*sqrt(1 - 1/(a**2*x**2))*(m + 1)) - 3*x**m*sqrt(-a**2*c*x**2 + c)/(a*sqrt(1 - 1/(a**2*x**2))*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_722():
    f = (-a**2*c*x**2 + c)**3*exp(n*acoth(a*x))
    F = -256*c**3*(1 - 1/(a*x))**(4 - n/2)*(1 + 1/(a*x))**(n/2 - 4)*hyper((8, 4 - n/2), (5 - n/2,), (a - 1/x)/(a + 1/x))/(a*(8 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_723():
    f = (-a**2*c*x**2 + c)**2*exp(n*acoth(a*x))
    F = 64*c**2*(1 - 1/(a*x))**(3 - n/2)*(1 + 1/(a*x))**(n/2 - 3)*hyper((6, 3 - n/2), (4 - n/2,), (a - 1/x)/(a + 1/x))/(a*(6 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_724():
    f = (-a**2*c*x**2 + c)*exp(n*acoth(a*x))
    F = -16*c*(1 - 1/(a*x))**(2 - n/2)*(1 + 1/(a*x))**(n/2 - 2)*hyper((4, 2 - n/2), (3 - n/2,), (a - 1/x)/(a + 1/x))/(a*(4 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_725():
    f = exp(n*acoth(a*x))
    F = 4*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (a - 1/x)/(a + 1/x))/(a*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_726():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)
    F = exp(n*acoth(a*x))/(a*c*n)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_727():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**2
    F = -(-2*a*x + n)*exp(n*acoth(a*x))/(a*c**2*(4 - n**2)*(-a**2*x**2 + 1)) + 2*exp(n*acoth(a*x))/(a*c**2*n*(4 - n**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_728():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**3
    F = -(-4*a*x + n)*exp(n*acoth(a*x))/(a*c**3*(16 - n**2)*(-a**2*x**2 + 1)**2) - 12*(-2*a*x + n)*exp(n*acoth(a*x))/(a*c**3*(4 - n**2)*(16 - n**2)*(-a**2*x**2 + 1)) + 24*exp(n*acoth(a*x))/(a*c**3*n*(n**4 - 20*n**2 + 64))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_729():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**4
    F = -(-6*a*x + n)*exp(n*acoth(a*x))/(a*c**4*(36 - n**2)*(-a**2*x**2 + 1)**3) - 30*(-4*a*x + n)*exp(n*acoth(a*x))/(a*c**4*(16 - n**2)*(36 - n**2)*(-a**2*x**2 + 1)**2) - 360*(-2*a*x + n)*exp(n*acoth(a*x))/(a*c**4*(4 - n**2)*(16 - n**2)*(36 - n**2)*(-a**2*x**2 + 1)) + 720*exp(n*acoth(a*x))/(a*c**4*n*(36 - n**2)*(n**4 - 20*n**2 + 64))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_730():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(n*acoth(a*x))
    F = 32*(1 - 1/(a*x))**(sympy.S(5)/2 - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-5)/2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)*hyper((5, sympy.S(5)/2 - n/2), (sympy.S(7)/2 - n/2,), (a - 1/x)/(a + 1/x))/(a**4*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(5 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_731():
    f = sqrt(-a**2*c*x**2 + c)*exp(n*acoth(a*x))
    F = 8*(1 - 1/(a*x))**(sympy.S(3)/2 - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*c*x**2 + c)*hyper((3, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), (a - 1/x)/(a + 1/x))/(a**2*x*sqrt(1 - 1/(a**2*x**2))*(3 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_732():
    f = exp(n*acoth(a*x))/sqrt(-a**2*c*x**2 + c)
    F = 2*x*sqrt(1 - 1/(a**2*x**2))*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (a - 1/x)/(a + 1/x))/((1 - n)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_733():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -(-a*x + n)*exp(n*acoth(a*x))/(a*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_734():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -(-3*a*x + n)*exp(n*acoth(a*x))/(a*c*(9 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 6*(-a*x + n)*exp(n*acoth(a*x))/(a*c**2*(1 - n**2)*(9 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_735():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -(-5*a*x + n)*exp(n*acoth(a*x))/(a*c*(25 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 20*(-3*a*x + n)*exp(n*acoth(a*x))/(a*c**2*(9 - n**2)*(25 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 120*(-a*x + n)*exp(n*acoth(a*x))/(a*c**3*(1 - n**2)*(9 - n**2)*(25 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_736():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(9)/2)
    F = -(-7*a*x + n)*exp(n*acoth(a*x))/(a*c*(49 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - 42*(-5*a*x + n)*exp(n*acoth(a*x))/(a*c**2*(25 - n**2)*(49 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 840*(-3*a*x + n)*exp(n*acoth(a*x))/(a*c**3*(9 - n**2)*(25 - n**2)*(49 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 5040*(-a*x + n)*exp(n*acoth(a*x))/(a*c**4*(1 - n**2)*(9 - n**2)*(25 - n**2)*(49 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_737():
    f = x**3*exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x**4*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)/(-a**2*c*x**2 + c)**(sympy.S(3)/2) - 2*n*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a + 1/x)/(a - 1/x))/(a*(1 - n)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*(n + 2)/(a*(n + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*(n**2 + 2*n + 2)/(a*(1 - n)*(n + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_738():
    f = x**2*exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -2*x*sqrt(1 - 1/(a**2*x**2))*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (a - 1/x)/(a + 1/x))/(a**2*c*(1 - n)*sqrt(-a**2*c*x**2 + c)) - (-a*x + n)*exp(n*acoth(a*x))/(a**3*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_739():
    f = x*exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = (-a*n*x + 1)*exp(n*acoth(a*x))/(a**2*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_740():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -(-a*x + n)*exp(n*acoth(a*x))/(a*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_741():
    f = exp(n*acoth(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = -2**(n/2 + sympy.S.Half)*a**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*hyper((sympy.S.Half - n/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (a - 1/x)/(2*a))/((1 - n)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - a**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)/((n + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**3*x**3*(1 - 1/(a**2*x**2))**(sympy.S(3)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)/((1 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_742():
    f = x**4*exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S(3)/2 - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)*(-n**3 - 2*n**2 + 7*n + 18)/((-a**2*c*x**2 + c)**(sympy.S(5)/2)*(n**4 - 10*n**2 + 9)) - x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-3)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)*(n + 6)/((n + 1)*(n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)*(n**2 + 6*n + 15)/((1 - n)*(n + 1)*(n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 2*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a + 1/x)/(a - 1/x))/((1 - n)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_743():
    f = x**3*exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -6*a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S(3)/2 - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((-a**2*c*x**2 + c)**(sympy.S(5)/2)*(n**4 - 10*n**2 + 9)) - a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-3)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 3*a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((-a**2*c*x**2 + c)**(sympy.S(5)/2)*(n**2 + 4*n + 3)) + 6*a*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((1 - n**2)*(n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_744():
    f = x**2*exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -(-3*a*x + n)*exp(n*acoth(a*x))/(a**3*c*(9 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + (3 - n**2)*(-a*x + n)*exp(n*acoth(a*x))/(a**3*c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_745():
    f = x*exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = (-a*n*x + 3)*exp(n*acoth(a*x))/(a**2*c*(9 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + 2*n*(-a*x + n)*exp(n*acoth(a*x))/(a**2*c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_746():
    f = exp(n*acoth(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -(-3*a*x + n)*exp(n*acoth(a*x))/(a*c*(9 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 6*(-a*x + n)*exp(n*acoth(a*x))/(a*c**2*(1 - n**2)*(9 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_747():
    f = exp(n*acoth(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = -2**(n/2 + sympy.S(5)/2)*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-3)/2)*hyper((-n/2 + sympy.S(-3)/2, -n/2 + sympy.S(-3)/2), (-n/2 + sympy.S(-1)/2,), (a - 1/x)/(2*a))/((n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 6*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S(3)/2 - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((-a**2*c*x**2 + c)**(sympy.S(5)/2)*(n**4 - 10*n**2 + 9)) - a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-3)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + 4*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-3)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)/((n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 6*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-3)/2)*(1 + 1/(a*x))**(n/2 + sympy.S.Half)/((n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + 4*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-3)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(3)/2)/((n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 3*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((-a**2*c*x**2 + c)**(sympy.S(5)/2)*(n**2 + 4*n + 3)) + 8*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)/((-a**2*c*x**2 + c)**(sympy.S(5)/2)*(n**2 + 4*n + 3)) - 6*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(-n/2 + sympy.S(-1)/2)*(1 + 1/(a*x))**(n/2 + sympy.S.Half)/((-a**2*c*x**2 + c)**(sympy.S(5)/2)*(n**2 + 4*n + 3)) + 6*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-3)/2)/((1 - n**2)*(n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 8*a**5*x**5*(1 - 1/(a**2*x**2))**(sympy.S(5)/2)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)/((1 - n**2)*(n + 3)*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_748():
    f = (-a**2*c*x**2 + c)**p*exp(n*acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(n/2 - p)*(1 - 1/(a*x))**(-n/2 + p)*(1 + 1/(a*x))**(n/2 + p + 1)*(-a**2*c*x**2 + c)**p*hyper((-2*p - 1, n/2 - p), (-2*p,), 2/(x*(a + 1/x)))/((1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_749():
    f = (-a**2*c*x**2 + c)**p*exp(2*p*acoth(a*x))
    F = x*(1 + 1/(a*x))**(2*p + 1)*(-a**2*c*x**2 + c)**p/((1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_750():
    f = (-a**2*c*x**2 + c)**p*exp(-2*p*acoth(a*x))
    F = x*(1 - 1/(a*x))**(2*p + 1)*(-a**2*c*x**2 + c)**p/((1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_751():
    f = (-a**2*c*x**2 + c)**p*exp(4*acoth(a*x))
    F = 2**(p + 2)*c*(a*x + 1)**(1 - p)*(-a**2*c*x**2 + c)**(p - 1)*hyper((p - 1, -p - 2), (p,), -a*x/2 + sympy.S.Half)/(a*(1 - p))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_752():
    f = (-a**2*c*x**2 + c)**p*exp(3*acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(sympy.S(3)/2 - p)*(1 - 1/(a*x))**(p + sympy.S(-3)/2)*(1 + 1/(a*x))**(p + sympy.S(5)/2)*(-a**2*c*x**2 + c)**p*hyper((sympy.S(3)/2 - p, -2*p - 1), (-2*p,), 2/(x*(a + 1/x)))/((1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_753():
    f = (-a**2*c*x**2 + c)**p*exp(2*acoth(a*x))
    F = 2**(p + 1)*(-a**2*c*x**2 + c)**p*hyper((p, -p - 1), (p + 1,), -a*x/2 + sympy.S.Half)/(a*p*(a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_754():
    f = (-a**2*c*x**2 + c)**p*exp(acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(sympy.S.Half - p)*(1 - 1/(a*x))**(p + sympy.S(-1)/2)*(1 + 1/(a*x))**(p + sympy.S(3)/2)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half - p, -2*p - 1), (-2*p,), 2/(x*(a + 1/x)))/((1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_755():
    f = (-a**2*c*x**2 + c)**p*exp(-acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(-p + sympy.S(-1)/2)*(1 - 1/(a*x))**(p + sympy.S.Half)*(1 + 1/(a*x))**(p + sympy.S.Half)*(-a**2*c*x**2 + c)**p*hyper((-2*p - 1, -p + sympy.S(-1)/2), (-2*p,), 2/(x*(a + 1/x)))/((1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_756():
    f = (-a**2*c*x**2 + c)**p*exp(-2*acoth(a*x))
    F = -2**(p + 1)*(-a**2*c*x**2 + c)**p*hyper((p, -p - 1), (p + 1,), a*x/2 + sympy.S.Half)/(a*p*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_757():
    f = (-a**2*c*x**2 + c)**p*exp(-3*acoth(a*x))
    F = x*((a - 1/x)/(a + 1/x))**(-p + sympy.S(-3)/2)*(1 - 1/(a*x))**(p + sympy.S(3)/2)*(1 + 1/(a*x))**(p + sympy.S(-1)/2)*(-a**2*c*x**2 + c)**p*hyper((-2*p - 1, -p + sympy.S(-3)/2), (-2*p,), 2/(x*(a + 1/x)))/((1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_758():
    f = (c - c/(a**2*x**2))**4*exp(acoth(a*x))
    F = c**4*x*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(9)/2) + 8*c**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/(7*a) + 47*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)/(42*a) + 61*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/(70*a) - 131*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/(280*a) - 91*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(120*a) - 67*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(48*a) - 51*c**4*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(16*a) + 35*c**4*acsc(a*x)/(16*a) + c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_759():
    f = (c - c/(a**2*x**2))**3*exp(acoth(a*x))
    F = c**3*x*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(7)/2) + 6*c**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/(5*a) + 23*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/(20*a) - 43*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(60*a) - 31*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(24*a) - 23*c**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(8*a) + 15*c**3*acsc(a*x)/(8*a) + c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_760():
    f = (c - c/(a**2*x**2))**2*exp(acoth(a*x))
    F = c**2*x*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2) + 4*c**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(3*a) - 7*c**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(6*a) - 5*c**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(2*a) + 3*c**2*acsc(a*x)/(2*a) + c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_761():
    f = (c - c/(a**2*x**2))*exp(acoth(a*x))
    F = c*x*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2) - 2*c*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/a + c*acsc(a*x)/a + c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_762():
    f = exp(acoth(a*x))/(c - c/(a**2*x**2))
    F = x*sqrt(1 + 1/(a*x))/(c*sqrt(1 - 1/(a*x))) + atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c) - 2*sqrt(1 + 1/(a*x))/(a*c*sqrt(1 - 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_763():
    f = exp(acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x/(c**2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))) + 8*sqrt(1 - 1/(a*x))/(3*a*c**2*sqrt(1 + 1/(a*x))) + atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**2) - 11/(3*a*c**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 4/(3*a*c**2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_764():
    f = exp(acoth(a*x))/(c - c/(a**2*x**2))**3
    F = x/(c**3*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)) + 16*sqrt(1 - 1/(a*x))/(5*a*c**3*sqrt(1 + 1/(a*x))) + 21*sqrt(1 - 1/(a*x))/(5*a*c**3*(1 + 1/(a*x))**(sympy.S(3)/2)) + atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**3) - 34/(5*a*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)) - 29/(15*a*c**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)) - 6/(5*a*c**3*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_765():
    f = exp(acoth(a*x))/(c - c/(a**2*x**2))**4
    F = x/(c**4*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)) + 128*sqrt(1 - 1/(a*x))/(35*a*c**4*sqrt(1 + 1/(a*x))) + 163*sqrt(1 - 1/(a*x))/(35*a*c**4*(1 + 1/(a*x))**(sympy.S(3)/2)) + 262*sqrt(1 - 1/(a*x))/(35*a*c**4*(1 + 1/(a*x))**(sympy.S(5)/2)) + atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**4) - 269/(21*a*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)) - 62/(21*a*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)) - 11/(7*a*c**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)) - 8/(7*a*c**4*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_766():
    f = (c - c/(a**2*x**2))**5*exp(2*acoth(a*x))
    F = c**5*x + 2*c**5*log(x)/a + 3*c**5/(a**2*x) + 4*c**5/(a**3*x**2) - 2*c**5/(3*a**4*x**3) - 3*c**5/(a**5*x**4) - 2*c**5/(5*a**6*x**5) + 4*c**5/(3*a**7*x**6) + 3*c**5/(7*a**8*x**7) - c**5/(4*a**9*x**8) - c**5/(9*a**10*x**9)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_767():
    f = (c - c/(a**2*x**2))**4*exp(2*acoth(a*x))
    F = c**4*x + 2*c**4*log(x)/a + 2*c**4/(a**2*x) + 3*c**4/(a**3*x**2) - 3*c**4/(2*a**5*x**4) - 2*c**4/(5*a**6*x**5) + c**4/(3*a**7*x**6) + c**4/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_768():
    f = (c - c/(a**2*x**2))**3*exp(2*acoth(a*x))
    F = c**3*x + 2*c**3*log(x)/a + c**3/(a**2*x) + 2*c**3/(a**3*x**2) + c**3/(3*a**4*x**3) - c**3/(2*a**5*x**4) - c**3/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_769():
    f = (c - c/(a**2*x**2))**2*exp(2*acoth(a*x))
    F = c**2*x + 2*c**2*log(x)/a + c**2/(a**3*x**2) + c**2/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_770():
    f = (c - c/(a**2*x**2))*exp(2*acoth(a*x))
    F = c*x + 2*c*log(x)/a - c/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_771():
    f = exp(2*acoth(a*x))/(c - c/(a**2*x**2))
    F = x/c + 2*log(-a*x + 1)/(a*c) + 1/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_772():
    f = exp(2*acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x/c**2 + 17*log(-a*x + 1)/(8*a*c**2) - log(a*x + 1)/(8*a*c**2) + 7/(4*a*c**2*(-a*x + 1)) - 1/(4*a*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_773():
    f = exp(2*acoth(a*x))/(c - c/(a**2*x**2))**3
    F = x/c**3 + 9*log(-a*x + 1)/(4*a*c**3) - log(a*x + 1)/(4*a*c**3) - 1/(16*a*c**3*(a*x + 1)) + 39/(16*a*c**3*(-a*x + 1)) - 5/(8*a*c**3*(-a*x + 1)**2) + 1/(12*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_774():
    f = exp(2*acoth(a*x))/(c - c/(a**2*x**2))**4
    F = x/c**4 + 303*log(-a*x + 1)/(128*a*c**4) - 47*log(a*x + 1)/(128*a*c**4) - 11/(64*a*c**4*(a*x + 1)) + 1/(64*a*c**4*(a*x + 1)**2) + 99/(32*a*c**4*(-a*x + 1)) - 35/(32*a*c**4*(-a*x + 1)**2) + 13/(48*a*c**4*(-a*x + 1)**3) - 1/(32*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_775():
    f = (c - c/(a**2*x**2))**4*exp(3*acoth(a*x))
    F = c**4*x*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(11)/2) + 8*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(11)/2)/(7*a) + 15*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(11)/2)/(14*a) - 57*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/(70*a) - 303*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/(280*a) - 61*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(40*a) - 37*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(16*a) - 63*c**4*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(16*a) + 15*c**4*acsc(a*x)/(16*a) + 3*c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_776():
    f = (c - c/(a**2*x**2))**3*exp(3*acoth(a*x))
    F = c**3*x*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(9)/2) + 6*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)/(5*a) - 21*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/(20*a) - 29*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(20*a) - 17*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(8*a) - 27*c**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(8*a) + 3*c**3*acsc(a*x)/(8*a) + 3*c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_777():
    f = (c - c/(a**2*x**2))**2*exp(3*acoth(a*x))
    F = c**2*x*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2) - 4*c**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(3*a) - 11*c**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(6*a) - 5*c**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(2*a) - c**2*acsc(a*x)/(2*a) + 3*c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_778():
    f = (c - c/(a**2*x**2))*exp(3*acoth(a*x))
    F = c*x*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2) - 3*c*acsc(a*x)/a + 3*c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_779():
    f = exp(3*acoth(a*x))/(c - c/(a**2*x**2))
    F = x*sqrt(1 + 1/(a*x))/(c*(1 - 1/(a*x))**(sympy.S(3)/2)) + 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c) - 14*sqrt(1 + 1/(a*x))/(3*a*c*sqrt(1 - 1/(a*x))) - 5*sqrt(1 + 1/(a*x))/(3*a*c*(1 - 1/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_780():
    f = exp(3*acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x*sqrt(1 + 1/(a*x))/(c**2*(1 - 1/(a*x))**(sympy.S(5)/2)) + 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**2) - 24*sqrt(1 + 1/(a*x))/(5*a*c**2*sqrt(1 - 1/(a*x))) - 9*sqrt(1 + 1/(a*x))/(5*a*c**2*(1 - 1/(a*x))**(sympy.S(3)/2)) - 6*sqrt(1 + 1/(a*x))/(5*a*c**2*(1 - 1/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_781():
    f = exp(3*acoth(a*x))/(c - c/(a**2*x**2))**3
    F = x/(c**3*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x))) + 176*sqrt(1 - 1/(a*x))/(35*a*c**3*sqrt(1 + 1/(a*x))) + 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**3) - 281/(35*a*c**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))) - 88/(35*a*c**3*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))) - 53/(35*a*c**3*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))) - 8/(7*a*c**3*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_782():
    f = exp(3*acoth(a*x))/(c - c/(a**2*x**2))**4
    F = x/(c**4*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)) + 1664*sqrt(1 - 1/(a*x))/(315*a*c**4*sqrt(1 + 1/(a*x))) + 2609*sqrt(1 - 1/(a*x))/(315*a*c**4*(1 + 1/(a*x))**(sympy.S(3)/2)) + 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**4) - 1462/(105*a*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)) - 1147/(315*a*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)) - 208/(105*a*c**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)) - 29/(21*a*c**4*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)) - 10/(9*a*c**4*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_783():
    f = (c - c/(a**2*x**2))**5*exp(4*acoth(a*x))
    F = c**5*x + 4*c**5*log(x)/a - 3*c**5/(a**2*x) + 4*c**5/(a**3*x**2) + 14*c**5/(3*a**4*x**3) - 14*c**5/(5*a**6*x**5) - 4*c**5/(3*a**7*x**6) + 3*c**5/(7*a**8*x**7) + c**5/(2*a**9*x**8) + c**5/(9*a**10*x**9)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_784():
    f = (c - c/(a**2*x**2))**4*exp(4*acoth(a*x))
    F = c**4*x + 4*c**4*log(x)/a - 4*c**4/(a**2*x) + 2*c**4/(a**3*x**2) + 10*c**4/(3*a**4*x**3) + c**4/(a**5*x**4) - 4*c**4/(5*a**6*x**5) - 2*c**4/(3*a**7*x**6) - c**4/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_785():
    f = (c - c/(a**2*x**2))**3*exp(4*acoth(a*x))
    F = c**3*x + 4*c**3*log(x)/a - 5*c**3/(a**2*x) + 5*c**3/(3*a**4*x**3) + c**3/(a**5*x**4) + c**3/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_786():
    f = (c - c/(a**2*x**2))**2*exp(4*acoth(a*x))
    F = c**2*x + 4*c**2*log(x)/a - 6*c**2/(a**2*x) - 2*c**2/(a**3*x**2) - c**2/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_787():
    f = (c - c/(a**2*x**2))*exp(4*acoth(a*x))
    F = c*x - 4*c*log(x)/a + 8*c*log(-a*x + 1)/a + c/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_788():
    f = exp(4*acoth(a*x))/(c - c/(a**2*x**2))
    F = x/c + 4*log(-a*x + 1)/(a*c) + 5/(a*c*(-a*x + 1)) - 1/(a*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_789():
    f = exp(4*acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x/c**2 + 4*log(-a*x + 1)/(a*c**2) + 6/(a*c**2*(-a*x + 1)) - 2/(a*c**2*(-a*x + 1)**2) + 1/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_790():
    f = exp(4*acoth(a*x))/(c - c/(a**2*x**2))**3
    F = x/c**3 + 129*log(-a*x + 1)/(32*a*c**3) - log(a*x + 1)/(32*a*c**3) + 111/(16*a*c**3*(-a*x + 1)) - 49/(16*a*c**3*(-a*x + 1)**2) + 11/(12*a*c**3*(-a*x + 1)**3) - 1/(8*a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_791():
    f = exp(4*acoth(a*x))/(c - c/(a**2*x**2))**4
    F = x/c**4 + 261*log(-a*x + 1)/(64*a*c**4) - 5*log(a*x + 1)/(64*a*c**4) - 1/(64*a*c**4*(a*x + 1)) + 501/(64*a*c**4*(-a*x + 1)) - 67/(16*a*c**4*(-a*x + 1)**2) + 83/(48*a*c**4*(-a*x + 1)**3) - 7/(16*a*c**4*(-a*x + 1)**4) + 1/(20*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_792():
    f = (c - c/(a**2*x**2))**4*exp(-acoth(a*x))
    F = c**4*x*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(7)/2) + 8*c**4*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/(7*a) + 7*c**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/(6*a) + 29*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)/(30*a) + 19*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)/(40*a) + 7*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(40*a) - c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(16*a) - 19*c**4*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(16*a) + 35*c**4*acsc(a*x)/(16*a) - c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_793():
    f = (c - c/(a**2*x**2))**3*exp(-acoth(a*x))
    F = c**3*x*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(5)/2) + 6*c**3*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(5*a) + 5*c**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(4*a) + 11*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(12*a) + c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(24*a) - 7*c**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(8*a) + 15*c**3*acsc(a*x)/(8*a) - c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_794():
    f = (c - c/(a**2*x**2))**2*exp(-acoth(a*x))
    F = c**2*x*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2) + 4*c**2*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(3*a) + 3*c**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(2*a) - c**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(2*a) + 3*c**2*acsc(a*x)/(2*a) - c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_795():
    f = (c - c/(a**2*x**2))*exp(-acoth(a*x))
    F = c*x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x)) + 2*c*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/a + c*acsc(a*x)/a - c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_796():
    f = exp(-acoth(a*x))/(c - c/(a**2*x**2))
    F = x*sqrt(1 - 1/(a*x))/(c*sqrt(1 + 1/(a*x))) + 2*sqrt(1 - 1/(a*x))/(a*c*sqrt(1 + 1/(a*x))) - atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_797():
    f = exp(-acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x/(c**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)) + 8*sqrt(1 - 1/(a*x))/(3*a*c**2*sqrt(1 + 1/(a*x))) + 5*sqrt(1 - 1/(a*x))/(3*a*c**2*(1 + 1/(a*x))**(sympy.S(3)/2)) - atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**2) - 2/(a*c**2*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_798():
    f = exp(-acoth(a*x))/(c - c/(a**2*x**2))**3
    F = x/(c**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)) + 16*sqrt(1 - 1/(a*x))/(5*a*c**3*sqrt(1 + 1/(a*x))) + 11*sqrt(1 - 1/(a*x))/(5*a*c**3*(1 + 1/(a*x))**(sympy.S(3)/2)) + 14*sqrt(1 - 1/(a*x))/(5*a*c**3*(1 + 1/(a*x))**(sympy.S(5)/2)) - atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**3) - 13/(3*a*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)) - 4/(3*a*c**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_799():
    f = exp(-acoth(a*x))/(c - c/(a**2*x**2))**4
    F = x/(c**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)) + 128*sqrt(1 - 1/(a*x))/(35*a*c**4*sqrt(1 + 1/(a*x))) + 93*sqrt(1 - 1/(a*x))/(35*a*c**4*(1 + 1/(a*x))**(sympy.S(3)/2)) + 122*sqrt(1 - 1/(a*x))/(35*a*c**4*(1 + 1/(a*x))**(sympy.S(5)/2)) + 115*sqrt(1 - 1/(a*x))/(21*a*c**4*(1 + 1/(a*x))**(sympy.S(7)/2)) - atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**4) - 28/(3*a*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)) - 31/(15*a*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(7)/2)) - 6/(5*a*c**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_800():
    f = (c - c/(a**2*x**2))**4*exp(-2*acoth(a*x))
    F = c**4*x - 2*c**4*log(x)/a + 2*c**4/(a**2*x) - 3*c**4/(a**3*x**2) + 3*c**4/(2*a**5*x**4) - 2*c**4/(5*a**6*x**5) - c**4/(3*a**7*x**6) + c**4/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_801():
    f = (c - c/(a**2*x**2))**3*exp(-2*acoth(a*x))
    F = c**3*x - 2*c**3*log(x)/a + c**3/(a**2*x) - 2*c**3/(a**3*x**2) + c**3/(3*a**4*x**3) + c**3/(2*a**5*x**4) - c**3/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_802():
    f = (c - c/(a**2*x**2))**2*exp(-2*acoth(a*x))
    F = c**2*x - 2*c**2*log(x)/a - c**2/(a**3*x**2) + c**2/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_803():
    f = (c - c/(a**2*x**2))*exp(-2*acoth(a*x))
    F = c*x - 2*c*log(x)/a - c/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_804():
    f = exp(-2*acoth(a*x))/(c - c/(a**2*x**2))
    F = x/c - 2*log(a*x + 1)/(a*c) - 1/(a*c*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_805():
    f = exp(-2*acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x/c**2 + log(-a*x + 1)/(8*a*c**2) - 17*log(a*x + 1)/(8*a*c**2) - 7/(4*a*c**2*(a*x + 1)) + 1/(4*a*c**2*(a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_806():
    f = exp(-2*acoth(a*x))/(c - c/(a**2*x**2))**3
    F = x/c**3 + log(-a*x + 1)/(4*a*c**3) - 9*log(a*x + 1)/(4*a*c**3) - 39/(16*a*c**3*(a*x + 1)) + 5/(8*a*c**3*(a*x + 1)**2) - 1/(12*a*c**3*(a*x + 1)**3) + 1/(16*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_807():
    f = exp(-2*acoth(a*x))/(c - c/(a**2*x**2))**4
    F = x/c**4 + 47*log(-a*x + 1)/(128*a*c**4) - 303*log(a*x + 1)/(128*a*c**4) - 99/(32*a*c**4*(a*x + 1)) + 35/(32*a*c**4*(a*x + 1)**2) - 13/(48*a*c**4*(a*x + 1)**3) + 1/(32*a*c**4*(a*x + 1)**4) + 11/(64*a*c**4*(-a*x + 1)) - 1/(64*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_808():
    f = (c - c/(a**2*x**2))**4*exp(-3*acoth(a*x))
    F = c**4*x*(1 - 1/(a*x))**(sympy.S(11)/2)*(1 + 1/(a*x))**(sympy.S(5)/2) + 8*c**4*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(7*a) + 17*c**4*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(14*a) + 11*c**4*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(10*a) + 5*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(5)/2)/(8*a) - 3*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(5)/2)/(8*a) + 27*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(16*a) + 33*c**4*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(16*a) + 15*c**4*acsc(a*x)/(16*a) - 3*c**4*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_809():
    f = (c - c/(a**2*x**2))**3*exp(-3*acoth(a*x))
    F = c**3*x*(1 - 1/(a*x))**(sympy.S(9)/2)*(1 + 1/(a*x))**(sympy.S(3)/2) + 6*c**3*(1 - 1/(a*x))**(sympy.S(7)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(5*a) + 27*c**3*(1 - 1/(a*x))**(sympy.S(5)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(20*a) + 5*c**3*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(3)/2)/(4*a) + 3*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(3)/2)/(8*a) + 21*c**3*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(8*a) + 3*c**3*acsc(a*x)/(8*a) - 3*c**3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_810():
    f = (c - c/(a**2*x**2))**2*exp(-3*acoth(a*x))
    F = c**2*x*(1 - 1/(a*x))**(sympy.S(7)/2)*sqrt(1 + 1/(a*x)) + 4*c**2*(1 - 1/(a*x))**(sympy.S(5)/2)*sqrt(1 + 1/(a*x))/(3*a) + 11*c**2*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x))/(6*a) + 5*c**2*sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x))/(2*a) - c**2*acsc(a*x)/(2*a) - 3*c**2*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_811():
    f = (c - c/(a**2*x**2))*exp(-3*acoth(a*x))
    F = c*x*(1 - 1/(a*x))**(sympy.S(3)/2)*sqrt(1 + 1/(a*x)) - 3*c*acsc(a*x)/a - 3*c*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_812():
    f = exp(-3*acoth(a*x))/(c - c/(a**2*x**2))
    F = x*sqrt(1 - 1/(a*x))/(c*(1 + 1/(a*x))**(sympy.S(3)/2)) + 14*sqrt(1 - 1/(a*x))/(3*a*c*sqrt(1 + 1/(a*x))) + 5*sqrt(1 - 1/(a*x))/(3*a*c*(1 + 1/(a*x))**(sympy.S(3)/2)) - 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_813():
    f = exp(-3*acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x*sqrt(1 - 1/(a*x))/(c**2*(1 + 1/(a*x))**(sympy.S(5)/2)) + 24*sqrt(1 - 1/(a*x))/(5*a*c**2*sqrt(1 + 1/(a*x))) + 9*sqrt(1 - 1/(a*x))/(5*a*c**2*(1 + 1/(a*x))**(sympy.S(3)/2)) + 6*sqrt(1 - 1/(a*x))/(5*a*c**2*(1 + 1/(a*x))**(sympy.S(5)/2)) - 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_814():
    f = exp(-3*acoth(a*x))/(c - c/(a**2*x**2))**3
    F = x/(c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2)) + 176*sqrt(1 - 1/(a*x))/(35*a*c**3*sqrt(1 + 1/(a*x))) + 71*sqrt(1 - 1/(a*x))/(35*a*c**3*(1 + 1/(a*x))**(sympy.S(3)/2)) + 54*sqrt(1 - 1/(a*x))/(35*a*c**3*(1 + 1/(a*x))**(sympy.S(5)/2)) + 11*sqrt(1 - 1/(a*x))/(7*a*c**3*(1 + 1/(a*x))**(sympy.S(7)/2)) - 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**3) - 2/(a*c**3*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_815():
    f = exp(-3*acoth(a*x))/(c - c/(a**2*x**2))**4
    F = x/(c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(9)/2)) + 1664*sqrt(1 - 1/(a*x))/(315*a*c**4*sqrt(1 + 1/(a*x))) + 719*sqrt(1 - 1/(a*x))/(315*a*c**4*(1 + 1/(a*x))**(sympy.S(3)/2)) + 202*sqrt(1 - 1/(a*x))/(105*a*c**4*(1 + 1/(a*x))**(sympy.S(5)/2)) + 139*sqrt(1 - 1/(a*x))/(63*a*c**4*(1 + 1/(a*x))**(sympy.S(7)/2)) + 28*sqrt(1 - 1/(a*x))/(9*a*c**4*(1 + 1/(a*x))**(sympy.S(9)/2)) - 3*atanh(sqrt(1 - 1/(a*x))*sqrt(1 + 1/(a*x)))/(a*c**4) - 5/(a*c**4*sqrt(1 - 1/(a*x))*(1 + 1/(a*x))**(sympy.S(9)/2)) - 4/(3*a*c**4*(1 - 1/(a*x))**(sympy.S(3)/2)*(1 + 1/(a*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_816():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(acoth(a*x))
    F = c**3*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + c**3*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + 3*c**3*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) + 3*c**3*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(c - c/(a**2*x**2))/(a**4*x**3*sqrt(1 - 1/(a**2*x**2))) - 3*c**3*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2))) + c**3*sqrt(c - c/(a**2*x**2))/(5*a**6*x**5*sqrt(1 - 1/(a**2*x**2))) + c**3*sqrt(c - c/(a**2*x**2))/(6*a**7*x**6*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_817():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(acoth(a*x))
    F = c**2*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + c**2*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + 2*c**2*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) + c**2*sqrt(c - c/(a**2*x**2))/(a**3*x**2*sqrt(1 - 1/(a**2*x**2))) - c**2*sqrt(c - c/(a**2*x**2))/(3*a**4*x**3*sqrt(1 - 1/(a**2*x**2))) - c**2*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_818():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(acoth(a*x))
    F = c*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + c*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + c*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) + c*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_819():
    f = sqrt(c - c/(a**2*x**2))*exp(acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_820():
    f = exp(acoth(a*x))/sqrt(c - c/(a**2*x**2))
    F = x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a**2*x**2)) + sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(a*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_821():
    f = exp(acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c*sqrt(c - c/(a**2*x**2))) + 5*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(4*a*c*sqrt(c - c/(a**2*x**2))) - sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(4*a*c*sqrt(c - c/(a**2*x**2))) + sqrt(1 - 1/(a**2*x**2))/(2*a*c*sqrt(c - c/(a**2*x**2))*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_822():
    f = exp(acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**2*sqrt(c - c/(a**2*x**2))) + 23*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) - 7*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) - sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + sqrt(1 - 1/(a**2*x**2))/(a*c**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1)) - sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_823():
    f = exp(acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**3*sqrt(c - c/(a**2*x**2))) + 51*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(32*a*c**3*sqrt(c - c/(a**2*x**2))) - 19*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(32*a*c**3*sqrt(c - c/(a**2*x**2))) - 5*sqrt(1 - 1/(a**2*x**2))/(16*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2) + 3*sqrt(1 - 1/(a**2*x**2))/(2*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)) - 11*sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2) + sqrt(1 - 1/(a**2*x**2))/(24*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_824():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(2*acoth(a*x))
    F = -57*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(16*(-a*x + 1)**3*(a*x + 1)**3) + 2*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) + 25*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(16*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) + 41*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(24*(-a*x + 1)**3*(a*x + 1)**2) + 57*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(80*(-a*x + 1)**3*(a*x + 1)) + 11*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(30*(-a*x + 1)**3) - 13*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)/(40*(-a*x + 1)**3) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)/(15*(-a*x + 1)**2) + x*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)/(-6*a*x + 6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_825():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(2*acoth(a*x))
    F = 25*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(a*x + 1)**2) - 2*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) - 9*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) - 17*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(12*(-a*x + 1)**2*(a*x + 1)) - 5*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)/(6*(-a*x + 1)**2) + x*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)/(-4*a*x + 4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_826():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(2*acoth(a*x))
    F = 2*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) + a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) - 5*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(a*x + 1)) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a*x + 1) + x*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(a*x + 1)/(-2*a*x + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_827():
    f = sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2)) - 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_828():
    f = exp(2*acoth(a*x))/sqrt(c - c/(a**2*x**2))
    F = -(-2*a*x + 2)*(a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))) + 2*sqrt(-a*x + 1)*sqrt(a*x + 1)*asin(a*x)/(a**2*x*sqrt(c - c/(a**2*x**2))) - (a*x + 1)**2/(a**2*x*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_829():
    f = exp(2*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = -(a*x + 1)**2/(3*a**2*x*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + (-4*a*x + 10)*(-a*x + 1)*(a*x + 1)**2/(3*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - 2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)*asin(a*x)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_830():
    f = exp(2*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = -(a*x + 1)**2/(5*a**2*x*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + (-2*a*x + 2)*(a*x + 1)**2/(3*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 58*(-a*x + 1)**2*(a*x + 1)**2/(15*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 2*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)*asin(a*x)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 2*(-a*x + 1)**3*(a*x + 1)**2*(43*a*x + 28)/(15*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_831():
    f = exp(2*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = -(a*x + 1)**2/(7*a**2*x*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + (-2*a*x + 2)*(a*x + 1)**2/(5*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 124*(-a*x + 1)**2*(a*x + 1)**2/(105*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 782*(-a*x + 1)**3*(a*x + 1)**2/(105*a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 142*(-a*x + 1)**4*(a*x + 1)**2/(35*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 2*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)*asin(a*x)/(a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 2*(-a*x + 1)**4*(a*x + 1)**3*(107*a*x + 72)/(35*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_832():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(3*acoth(a*x))
    F = c**4*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + 3*c**4*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + 4*c**4*sqrt(c - c/(a**2*x**2))/(a**3*x**2*sqrt(1 - 1/(a**2*x**2))) + 2*c**4*sqrt(c - c/(a**2*x**2))/(a**4*x**3*sqrt(1 - 1/(a**2*x**2))) - 3*c**4*sqrt(c - c/(a**2*x**2))/(2*a**5*x**4*sqrt(1 - 1/(a**2*x**2))) - 8*c**4*sqrt(c - c/(a**2*x**2))/(5*a**6*x**5*sqrt(1 - 1/(a**2*x**2))) + 3*c**4*sqrt(c - c/(a**2*x**2))/(7*a**8*x**7*sqrt(1 - 1/(a**2*x**2))) + c**4*sqrt(c - c/(a**2*x**2))/(8*a**9*x**8*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_833():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(3*acoth(a*x))
    F = c**3*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + 3*c**3*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) + 5*c**3*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2))) + 5*c**3*sqrt(c - c/(a**2*x**2))/(3*a**4*x**3*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2))) - 3*c**3*sqrt(c - c/(a**2*x**2))/(5*a**6*x**5*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(c - c/(a**2*x**2))/(6*a**7*x**6*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_834():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(3*acoth(a*x))
    F = c**2*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + 3*c**2*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - 2*c**2*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) + c**2*sqrt(c - c/(a**2*x**2))/(a**3*x**2*sqrt(1 - 1/(a**2*x**2))) + c**2*sqrt(c - c/(a**2*x**2))/(a**4*x**3*sqrt(1 - 1/(a**2*x**2))) + c**2*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_835():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(3*acoth(a*x))
    F = c*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + 3*c*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - 3*c*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) - c*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_836():
    f = sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_837():
    f = exp(3*acoth(a*x))/sqrt(c - c/(a**2*x**2))
    F = x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a**2*x**2)) + 3*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(a*sqrt(c - c/(a**2*x**2))) + 2*sqrt(1 - 1/(a**2*x**2))/(a*sqrt(c - c/(a**2*x**2))*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_838():
    f = exp(3*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c*sqrt(c - c/(a**2*x**2))) + 3*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(a*c*sqrt(c - c/(a**2*x**2))) + 3*sqrt(1 - 1/(a**2*x**2))/(a*c*sqrt(c - c/(a**2*x**2))*(-a*x + 1)) - sqrt(1 - 1/(a**2*x**2))/(2*a*c*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_839():
    f = exp(3*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**2*sqrt(c - c/(a**2*x**2))) + 49*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) - sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) + 31*sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1)) - 9*sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2) + sqrt(1 - 1/(a**2*x**2))/(6*a*c**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_840():
    f = exp(3*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**3*sqrt(c - c/(a**2*x**2))) + 201*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(64*a*c**3*sqrt(c - c/(a**2*x**2))) - 9*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(64*a*c**3*sqrt(c - c/(a**2*x**2))) - sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + 75*sqrt(1 - 1/(a**2*x**2))/(16*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)) - 59*sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2) + sqrt(1 - 1/(a**2*x**2))/(2*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**3) - sqrt(1 - 1/(a**2*x**2))/(16*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_841():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(-acoth(a*x))
    F = c**3*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - c**3*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + 3*c**3*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) - 3*c**3*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(c - c/(a**2*x**2))/(a**4*x**3*sqrt(1 - 1/(a**2*x**2))) + 3*c**3*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2))) + c**3*sqrt(c - c/(a**2*x**2))/(5*a**6*x**5*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(c - c/(a**2*x**2))/(6*a**7*x**6*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_842():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(-acoth(a*x))
    F = c**2*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - c**2*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + 2*c**2*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) - c**2*sqrt(c - c/(a**2*x**2))/(a**3*x**2*sqrt(1 - 1/(a**2*x**2))) - c**2*sqrt(c - c/(a**2*x**2))/(3*a**4*x**3*sqrt(1 - 1/(a**2*x**2))) + c**2*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_843():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(-acoth(a*x))
    F = c*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - c*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + c*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) - c*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_844():
    f = sqrt(c - c/(a**2*x**2))*exp(-acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_845():
    f = exp(-acoth(a*x))/sqrt(c - c/(a**2*x**2))
    F = x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a**2*x**2)) - sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(a*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_846():
    f = exp(-acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c*sqrt(c - c/(a**2*x**2))) + sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(4*a*c*sqrt(c - c/(a**2*x**2))) - 5*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(4*a*c*sqrt(c - c/(a**2*x**2))) - sqrt(1 - 1/(a**2*x**2))/(2*a*c*sqrt(c - c/(a**2*x**2))*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_847():
    f = exp(-acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**2*sqrt(c - c/(a**2*x**2))) + 7*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) - 23*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) - sqrt(1 - 1/(a**2*x**2))/(a*c**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2) + sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_848():
    f = exp(-acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**3*sqrt(c - c/(a**2*x**2))) + 19*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(32*a*c**3*sqrt(c - c/(a**2*x**2))) - 51*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(32*a*c**3*sqrt(c - c/(a**2*x**2))) - 3*sqrt(1 - 1/(a**2*x**2))/(2*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + 11*sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2) - sqrt(1 - 1/(a**2*x**2))/(24*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)**3) + 5*sqrt(1 - 1/(a**2*x**2))/(16*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)) - sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_849():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(-2*acoth(a*x))
    F = 7*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(16*(-a*x + 1)**3*(a*x + 1)**3) - 2*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) + 25*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(16*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) + 3*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(8*(-a*x + 1)**3*(a*x + 1)**2) - 19*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(16*(-a*x + 1)**3*(a*x + 1)) + 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(3*(-a*x + 1)**2*(a*x + 1)) - 23*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)/((-120*a*x + 120)*(a*x + 1)) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(15*a*x + 15) + x*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)/(6*a*x + 6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_850():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(-2*acoth(a*x))
    F = -7*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(a*x + 1)**2) + 2*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) - 9*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) + 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/((-a*x + 1)**2*(a*x + 1)) - 7*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/((-24*a*x + 24)*(a*x + 1)) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(6*a*x + 6) + x*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1)/(4*a*x + 4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_851():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(-2*acoth(a*x))
    F = -2*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) + a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) - 5*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(a*x + 1)) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(a*x + 1) + x*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)/(2*a*x + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_852():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2)) + 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_853():
    f = exp(-2*acoth(a*x))/sqrt(c - c/(a**2*x**2))
    F = -(-2*a*x + 2)*(a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))) - 2*sqrt(-a*x + 1)*sqrt(a*x + 1)*asin(a*x)/(a**2*x*sqrt(c - c/(a**2*x**2))) - (-a*x + 1)**2/(a**2*x*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_854():
    f = exp(-2*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = -(-a*x + 1)**2/(3*a**2*x*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + 2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)*asin(a*x)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + 2*(-a*x + 1)**2*(a*x + 1)*(2*a*x + 5)/(3*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_855():
    f = exp(-2*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = -(-a*x + 1)**2/(a**2*x*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 2*(-a*x + 1)**3/(5*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 2*(-a*x + 1)**3*(a*x + 1)/(15*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 2*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)*asin(a*x)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 2*(-a*x + 1)**3*(a*x + 1)**2*(13*a*x + 28)/(15*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_856():
    f = exp(-2*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = -(-a*x + 1)**2/(3*a**2*x*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 10*(-a*x + 1)**3/(3*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 12*(-a*x + 1)**4/(7*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 82*(-a*x + 1)**4*(a*x + 1)/(105*a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 2*(-a*x + 1)**4*(a*x + 1)**2/(35*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 2*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)*asin(a*x)/(a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 2*(-a*x + 1)**4*(a*x + 1)**3*(37*a*x + 72)/(35*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_857():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(-3*acoth(a*x))
    F = c**4*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - 3*c**4*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - 4*c**4*sqrt(c - c/(a**2*x**2))/(a**3*x**2*sqrt(1 - 1/(a**2*x**2))) + 2*c**4*sqrt(c - c/(a**2*x**2))/(a**4*x**3*sqrt(1 - 1/(a**2*x**2))) + 3*c**4*sqrt(c - c/(a**2*x**2))/(2*a**5*x**4*sqrt(1 - 1/(a**2*x**2))) - 8*c**4*sqrt(c - c/(a**2*x**2))/(5*a**6*x**5*sqrt(1 - 1/(a**2*x**2))) + 3*c**4*sqrt(c - c/(a**2*x**2))/(7*a**8*x**7*sqrt(1 - 1/(a**2*x**2))) - c**4*sqrt(c - c/(a**2*x**2))/(8*a**9*x**8*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_858():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(-3*acoth(a*x))
    F = c**3*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - 3*c**3*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - c**3*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) - 5*c**3*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2))) + 5*c**3*sqrt(c - c/(a**2*x**2))/(3*a**4*x**3*sqrt(1 - 1/(a**2*x**2))) + c**3*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2))) - 3*c**3*sqrt(c - c/(a**2*x**2))/(5*a**6*x**5*sqrt(1 - 1/(a**2*x**2))) + c**3*sqrt(c - c/(a**2*x**2))/(6*a**7*x**6*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_859():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(-3*acoth(a*x))
    F = c**2*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - 3*c**2*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - 2*c**2*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) - c**2*sqrt(c - c/(a**2*x**2))/(a**3*x**2*sqrt(1 - 1/(a**2*x**2))) + c**2*sqrt(c - c/(a**2*x**2))/(a**4*x**3*sqrt(1 - 1/(a**2*x**2))) - c**2*sqrt(c - c/(a**2*x**2))/(4*a**5*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_860():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(-3*acoth(a*x))
    F = c*x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - 3*c*sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - 3*c*sqrt(c - c/(a**2*x**2))/(a**2*x*sqrt(1 - 1/(a**2*x**2))) + c*sqrt(c - c/(a**2*x**2))/(2*a**3*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_861():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - 4*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_862():
    f = exp(-3*acoth(a*x))/sqrt(c - c/(a**2*x**2))
    F = x*sqrt(1 - 1/(a**2*x**2))/sqrt(c - c/(a**2*x**2)) - 3*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(a*sqrt(c - c/(a**2*x**2))) - 2*sqrt(1 - 1/(a**2*x**2))/(a*sqrt(c - c/(a**2*x**2))*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_863():
    f = exp(-3*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c*sqrt(c - c/(a**2*x**2))) - 3*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(a*c*sqrt(c - c/(a**2*x**2))) - 3*sqrt(1 - 1/(a**2*x**2))/(a*c*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + sqrt(1 - 1/(a**2*x**2))/(2*a*c*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_864():
    f = exp(-3*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**2*sqrt(c - c/(a**2*x**2))) + sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) - 49*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(16*a*c**2*sqrt(c - c/(a**2*x**2))) - 31*sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + 9*sqrt(1 - 1/(a**2*x**2))/(8*a*c**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2) - sqrt(1 - 1/(a**2*x**2))/(6*a*c**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_865():
    f = exp(-3*acoth(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = x*sqrt(1 - 1/(a**2*x**2))/(c**3*sqrt(c - c/(a**2*x**2))) + 9*sqrt(1 - 1/(a**2*x**2))*log(-a*x + 1)/(64*a*c**3*sqrt(c - c/(a**2*x**2))) - 201*sqrt(1 - 1/(a**2*x**2))*log(a*x + 1)/(64*a*c**3*sqrt(c - c/(a**2*x**2))) - 75*sqrt(1 - 1/(a**2*x**2))/(16*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)) + 59*sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2) - sqrt(1 - 1/(a**2*x**2))/(2*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)**3) + sqrt(1 - 1/(a**2*x**2))/(16*a*c**3*sqrt(c - c/(a**2*x**2))*(a*x + 1)**4) + sqrt(1 - 1/(a**2*x**2))/(32*a*c**3*sqrt(c - c/(a**2*x**2))*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_866():
    f = x**m*sqrt(c - c/(a**2*x**2))*exp(acoth(a*x))
    F = x**(m + 1)*sqrt(c - c/(a**2*x**2))/(sqrt(1 - 1/(a**2*x**2))*(m + 1)) + x**m*sqrt(c - c/(a**2*x**2))/(a*m*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_867():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(acoth(a*x))
    F = x**3*sqrt(c - c/(a**2*x**2))/(3*sqrt(1 - 1/(a**2*x**2))) + x**2*sqrt(c - c/(a**2*x**2))/(2*a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_868():
    f = x*sqrt(c - c/(a**2*x**2))*exp(acoth(a*x))
    F = x**2*sqrt(c - c/(a**2*x**2))/(2*sqrt(1 - 1/(a**2*x**2))) + x*sqrt(c - c/(a**2*x**2))/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_869():
    f = sqrt(c - c/(a**2*x**2))*exp(acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_870():
    f = sqrt(c - c/(a**2*x**2))*exp(acoth(a*x))/x
    F = sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) - sqrt(c - c/(a**2*x**2))/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_871():
    f = sqrt(c - c/(a**2*x**2))*exp(acoth(a*x))/x**2
    F = -sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(2*a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_872():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))
    F = x**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(4*a**2) + x*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(6*a**3) + 7*x*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(24*a**3) + 7*x*sqrt(c - c/(a**2*x**2))/(8*a**3) - 7*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(8*a**3*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_873():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(3*a**2) + x*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(3*a**2) + x*sqrt(c - c/(a**2*x**2))/a**2 - x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(a**2*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_874():
    f = x*sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(2*a) + 3*x*sqrt(c - c/(a**2*x**2))/(2*a) - 3*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(2*a*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_875():
    f = sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2)) - 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_876():
    f = sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))/x
    F = -a*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + 2*a*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + sqrt(c - c/(a**2*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_877():
    f = sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))/x**2
    F = 3*a**2*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 3*a*sqrt(c - c/(a**2*x**2))/2 + sqrt(c - c/(a**2*x**2))*(a*x + 1)/(2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_878():
    f = sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))/x**3
    F = a**3*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + a**2*sqrt(c - c/(a**2*x**2)) + a*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(3*x) + sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_879():
    f = sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))/x**4
    F = 7*a**4*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 4*a**3*sqrt(c - c/(a**2*x**2))/3 + 7*a**2*sqrt(c - c/(a**2*x**2))/(8*x) + 2*a*sqrt(c - c/(a**2*x**2))/(3*x**2) + sqrt(c - c/(a**2*x**2))/(4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_880():
    f = sqrt(c - c/(a**2*x**2))*exp(2*acoth(a*x))/x**5
    F = 3*a**5*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(4*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 6*a**4*sqrt(c - c/(a**2*x**2))/5 + 3*a**3*sqrt(c - c/(a**2*x**2))/(4*x) + 3*a**2*sqrt(c - c/(a**2*x**2))/(5*x**2) + a*sqrt(c - c/(a**2*x**2))/(2*x**3) + sqrt(c - c/(a**2*x**2))/(5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_881():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))
    F = x**4*sqrt(c - c/(a**2*x**2))/(4*sqrt(1 - 1/(a**2*x**2))) + x**3*sqrt(c - c/(a**2*x**2))/(a*sqrt(1 - 1/(a**2*x**2))) + 2*x**2*sqrt(c - c/(a**2*x**2))/(a**2*sqrt(1 - 1/(a**2*x**2))) + 4*x*sqrt(c - c/(a**2*x**2))/(a**3*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_882():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))
    F = x**3*sqrt(c - c/(a**2*x**2))/(3*sqrt(1 - 1/(a**2*x**2))) + 3*x**2*sqrt(c - c/(a**2*x**2))/(2*a*sqrt(1 - 1/(a**2*x**2))) + 4*x*sqrt(c - c/(a**2*x**2))/(a**2*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_883():
    f = x*sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))
    F = x**2*sqrt(c - c/(a**2*x**2))/(2*sqrt(1 - 1/(a**2*x**2))) + 3*x*sqrt(c - c/(a**2*x**2))/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_884():
    f = sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_885():
    f = sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))/x
    F = -3*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(1 - 1/(a**2*x**2)) + sqrt(c - c/(a**2*x**2))/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_886():
    f = sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))/x**2
    F = -4*a*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*a*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(1 - 1/(a**2*x**2)) + 3*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) + sqrt(c - c/(a**2*x**2))/(2*a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_887():
    f = sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))/x**3
    F = -4*a**2*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*a**2*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(1 - 1/(a**2*x**2)) + 4*a*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(c - c/(a**2*x**2))/(2*x**2*sqrt(1 - 1/(a**2*x**2))) + sqrt(c - c/(a**2*x**2))/(3*a*x**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_888():
    f = sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))/x**4
    F = -4*a**3*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*a**3*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(1 - 1/(a**2*x**2)) + 4*a**2*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) + 2*a*sqrt(c - c/(a**2*x**2))/(x**2*sqrt(1 - 1/(a**2*x**2))) + sqrt(c - c/(a**2*x**2))/(x**3*sqrt(1 - 1/(a**2*x**2))) + sqrt(c - c/(a**2*x**2))/(4*a*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_889():
    f = sqrt(c - c/(a**2*x**2))*exp(3*acoth(a*x))/x**5
    F = -4*a**4*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*a**4*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(1 - 1/(a**2*x**2)) + 4*a**3*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) + 2*a**2*sqrt(c - c/(a**2*x**2))/(x**2*sqrt(1 - 1/(a**2*x**2))) + 4*a*sqrt(c - c/(a**2*x**2))/(3*x**3*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(c - c/(a**2*x**2))/(4*x**4*sqrt(1 - 1/(a**2*x**2))) + sqrt(c - c/(a**2*x**2))/(5*a*x**5*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_890():
    f = x**m*sqrt(c - c/(a**2*x**2))*exp(-acoth(a*x))
    F = x**(m + 1)*sqrt(c - c/(a**2*x**2))/(sqrt(1 - 1/(a**2*x**2))*(m + 1)) - x**m*sqrt(c - c/(a**2*x**2))/(a*m*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_891():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(-acoth(a*x))
    F = x**3*sqrt(c - c/(a**2*x**2))/(3*sqrt(1 - 1/(a**2*x**2))) - x**2*sqrt(c - c/(a**2*x**2))/(2*a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_892():
    f = x*sqrt(c - c/(a**2*x**2))*exp(-acoth(a*x))
    F = x**2*sqrt(c - c/(a**2*x**2))/(2*sqrt(1 - 1/(a**2*x**2))) - x*sqrt(c - c/(a**2*x**2))/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_893():
    f = sqrt(c - c/(a**2*x**2))*exp(-acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) - sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_894():
    f = sqrt(c - c/(a**2*x**2))*exp(-acoth(a*x))/x
    F = sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + sqrt(c - c/(a**2*x**2))/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_895():
    f = sqrt(c - c/(a**2*x**2))*exp(-acoth(a*x))/x**2
    F = sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(2*a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_896():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))
    F = x**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(4*a**2) - x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(6*a**3) - 7*x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(24*a**3) - 7*x*sqrt(c - c/(a**2*x**2))/(8*a**3) - 7*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(8*a**3*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_897():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(3*a**2) + x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(3*a**2) + x*sqrt(c - c/(a**2*x**2))/a**2 + x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(a**2*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_898():
    f = x*sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))
    F = -x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(2*a) - 3*x*sqrt(c - c/(a**2*x**2))/(2*a) - 3*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(2*a*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_899():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2)) + 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_900():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))/x
    F = -a*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*a*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + sqrt(c - c/(a**2*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_901():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))/x**2
    F = 3*a**2*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 3*a*sqrt(c - c/(a**2*x**2))/2 + sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_902():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))/x**3
    F = -a**3*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + a**2*sqrt(c - c/(a**2*x**2)) - a*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(3*x) + sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_903():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))/x**4
    F = 7*a**4*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 4*a**3*sqrt(c - c/(a**2*x**2))/3 + 7*a**2*sqrt(c - c/(a**2*x**2))/(8*x) - 2*a*sqrt(c - c/(a**2*x**2))/(3*x**2) + sqrt(c - c/(a**2*x**2))/(4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_904():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*acoth(a*x))/x**5
    F = -3*a**5*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(4*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 6*a**4*sqrt(c - c/(a**2*x**2))/5 - 3*a**3*sqrt(c - c/(a**2*x**2))/(4*x) + 3*a**2*sqrt(c - c/(a**2*x**2))/(5*x**2) - a*sqrt(c - c/(a**2*x**2))/(2*x**3) + sqrt(c - c/(a**2*x**2))/(5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_905():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))
    F = x**4*sqrt(c - c/(a**2*x**2))/(4*sqrt(1 - 1/(a**2*x**2))) - x**3*sqrt(c - c/(a**2*x**2))/(a*sqrt(1 - 1/(a**2*x**2))) + 2*x**2*sqrt(c - c/(a**2*x**2))/(a**2*sqrt(1 - 1/(a**2*x**2))) - 4*x*sqrt(c - c/(a**2*x**2))/(a**3*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_906():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))
    F = x**3*sqrt(c - c/(a**2*x**2))/(3*sqrt(1 - 1/(a**2*x**2))) - 3*x**2*sqrt(c - c/(a**2*x**2))/(2*a*sqrt(1 - 1/(a**2*x**2))) + 4*x*sqrt(c - c/(a**2*x**2))/(a**2*sqrt(1 - 1/(a**2*x**2))) - 4*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_907():
    f = x*sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))
    F = x**2*sqrt(c - c/(a**2*x**2))/(2*sqrt(1 - 1/(a**2*x**2))) - 3*x*sqrt(c - c/(a**2*x**2))/(a*sqrt(1 - 1/(a**2*x**2))) + 4*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_908():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))
    F = x*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + sqrt(c - c/(a**2*x**2))*log(x)/(a*sqrt(1 - 1/(a**2*x**2))) - 4*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_909():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))/x
    F = -3*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(1 - 1/(a**2*x**2)) - sqrt(c - c/(a**2*x**2))/(a*x*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_910():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))/x**2
    F = 4*a*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) - 4*a*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(1 - 1/(a**2*x**2)) + 3*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) - sqrt(c - c/(a**2*x**2))/(2*a*x**2*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_911():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))/x**3
    F = -4*a**2*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*a**2*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(1 - 1/(a**2*x**2)) - 4*a*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(c - c/(a**2*x**2))/(2*x**2*sqrt(1 - 1/(a**2*x**2))) - sqrt(c - c/(a**2*x**2))/(3*a*x**3*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_912():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))/x**4
    F = 4*a**3*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) - 4*a**3*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(1 - 1/(a**2*x**2)) + 4*a**2*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) - 2*a*sqrt(c - c/(a**2*x**2))/(x**2*sqrt(1 - 1/(a**2*x**2))) + sqrt(c - c/(a**2*x**2))/(x**3*sqrt(1 - 1/(a**2*x**2))) - sqrt(c - c/(a**2*x**2))/(4*a*x**4*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_913():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*acoth(a*x))/x**5
    F = -4*a**4*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(1 - 1/(a**2*x**2)) + 4*a**4*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(1 - 1/(a**2*x**2)) - 4*a**3*sqrt(c - c/(a**2*x**2))/(x*sqrt(1 - 1/(a**2*x**2))) + 2*a**2*sqrt(c - c/(a**2*x**2))/(x**2*sqrt(1 - 1/(a**2*x**2))) - 4*a*sqrt(c - c/(a**2*x**2))/(3*x**3*sqrt(1 - 1/(a**2*x**2))) + 3*sqrt(c - c/(a**2*x**2))/(4*x**4*sqrt(1 - 1/(a**2*x**2))) - sqrt(c - c/(a**2*x**2))/(5*a*x**5*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_914():
    f = (c - c/(a**2*x**2))*exp(n*acoth(a*x))
    F = -2**(n/2 + 1)*c*(1 - 1/(a*x))**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), (a - 1/x)/(2*a))/(a*(2 - n)) + 4*c*(1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (a - 1/x)/(a + 1/x))/(a*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_915():
    f = exp(n*acoth(a*x))/(c - c/(a**2*x**2))
    F = x*(1 + 1/(a*x))**(n/2)/(c*(1 - 1/(a*x))**(n/2)) + 2*(1 + 1/(a*x))**(n/2)*hyper((1, n/2), (n/2 + 1,), (a + 1/x)/(a - 1/x))/(a*c*(1 - 1/(a*x))**(n/2)) - (1 + 1/(a*x))**(n/2)*(n + 1)/(a*c*n*(1 - 1/(a*x))**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_916():
    f = exp(n*acoth(a*x))/(c - c/(a**2*x**2))**2
    F = x*(1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 - 1)/c**2 - (1 - 1/(a*x))**(-n/2 - 1)*(1 + 1/(a*x))**(n/2 - 1)*(n + 3)/(a*c**2*(n + 2)) + 2*(1 + 1/(a*x))**(n/2)*hyper((1, n/2), (n/2 + 1,), (a + 1/x)/(a - 1/x))/(a*c**2*(1 - 1/(a*x))**(n/2)) + (1 - 1/(a*x))**(1 - n/2)*(1 + 1/(a*x))**(n/2 - 1)*(-n**3 - n**2 + 4*n + 6)/(a*c**2*n*(2 - n)*(n + 2)) - (1 + 1/(a*x))**(n/2 - 1)*(n**2 + 4*n + 6)/(a*c**2*n*(1 - 1/(a*x))**(n/2)*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_917():
    f = sqrt(c - c/(a**2*x**2))*exp(n*acoth(a*x))
    F = -2**(n/2 + sympy.S.Half)*(1 - 1/(a*x))**(sympy.S.Half - n/2)*sqrt(c - c/(a**2*x**2))*hyper((sympy.S.Half - n/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (a - 1/x)/(2*a))/(a*(1 - n)*sqrt(1 - 1/(a**2*x**2))) + x*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S.Half)*sqrt(c - c/(a**2*x**2))/sqrt(1 - 1/(a**2*x**2)) + 2*n*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*sqrt(c - c/(a**2*x**2))*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (a - 1/x)/(a + 1/x))/(a*(1 - n)*sqrt(1 - 1/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_918():
    f = exp(n*acoth(a*x))/sqrt(c - c/(a**2*x**2))
    F = x*sqrt(1 - 1/(a**2*x**2))*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S.Half)/sqrt(c - c/(a**2*x**2)) + 2*n*sqrt(1 - 1/(a**2*x**2))*(1 - 1/(a*x))**(sympy.S.Half - n/2)*(1 + 1/(a*x))**(n/2 + sympy.S(-1)/2)*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (a - 1/x)/(a + 1/x))/(a*(1 - n)*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_919():
    f = (c - c/(a**2*x**2))**p*exp(n*acoth(a*x))
    F = -2**(-n/2 + p + 1)*(1 + 1/(a*x))**(n/2 + p + 1)*(c - c/(a**2*x**2))**p*appellf1(n/2 + p + 1, 2, n/2 - p, n/2 + p + 2, 1 + 1/(a*x), (a + 1/x)/(2*a))/(a*(1 - 1/(a**2*x**2))**p*(n + 2*p + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_920():
    f = (c - c/(a**2*x**2))**p*exp(-2*p*acoth(a*x))
    F = (1 - 1/(a*x))**(2*p + 1)*(c - c/(a**2*x**2))**p*hyper((2, 2*p + 1), (2*p + 2,), 1 - 1/(a*x))/(a*(1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_4_Inverse_hyperbolic_cotangent_7_4_2_Exponentials_of_inverse_hyperbolic_cotangent_functions_921():
    f = (c - c/(a**2*x**2))**p*exp(2*p*acoth(a*x))
    F = -(1 + 1/(a*x))**(2*p + 1)*(c - c/(a**2*x**2))**p*hyper((2, 2*p + 1), (2*p + 2,), 1 + 1/(a*x))/(a*(1 - 1/(a**2*x**2))**p*(2*p + 1))
    assert integrate(f, x) == F

