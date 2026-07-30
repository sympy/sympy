"""Generated from MathematicaSyntaxTestSuite.

Source: 7 Inverse hyperbolic functions/7.3 Inverse hyperbolic tangent/7.3.6 Exponentials of inverse hyperbolic tangent functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, m, n, p = symbols('a b c m n p')

def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1():
    f = x**4*exp(atanh(a*x))
    F = -x**4*sqrt(-a**2*x**2 + 1)/(5*a) - x**3*sqrt(-a**2*x**2 + 1)/(4*a**2) - 4*x**2*sqrt(-a**2*x**2 + 1)/(15*a**3) - (45*a*x + 64)*sqrt(-a**2*x**2 + 1)/(120*a**5) + 3*asin(a*x)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_2():
    f = x**3*exp(atanh(a*x))
    F = -x**3*sqrt(-a**2*x**2 + 1)/(4*a) - x**2*sqrt(-a**2*x**2 + 1)/(3*a**2) - (9*a*x + 16)*sqrt(-a**2*x**2 + 1)/(24*a**4) + 3*asin(a*x)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_3():
    f = x**2*exp(atanh(a*x))
    F = -x*sqrt(-a**2*x**2 + 1)/(2*a**2) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**3) - sqrt(-a**2*x**2 + 1)/a**3 + asin(a*x)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_4():
    f = x*exp(atanh(a*x))
    F = -(a*x + 2)*sqrt(-a**2*x**2 + 1)/(2*a**2) + asin(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_5():
    f = exp(atanh(a*x))
    F = -sqrt(-a**2*x**2 + 1)/a + asin(a*x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_6():
    f = exp(atanh(a*x))/x
    F = asin(a*x) - atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_7():
    f = exp(atanh(a*x))/x**2
    F = -a*atanh(sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_8():
    f = exp(atanh(a*x))/x**3
    F = -a**2*atanh(sqrt(-a**2*x**2 + 1))/2 - a*sqrt(-a**2*x**2 + 1)/x - sqrt(-a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_9():
    f = exp(atanh(a*x))/x**4
    F = -a**3*atanh(sqrt(-a**2*x**2 + 1))/2 - 2*a**2*sqrt(-a**2*x**2 + 1)/(3*x) - a*sqrt(-a**2*x**2 + 1)/(2*x**2) - sqrt(-a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_10():
    f = exp(atanh(a*x))/x**5
    F = -3*a**4*atanh(sqrt(-a**2*x**2 + 1))/8 - 2*a**3*sqrt(-a**2*x**2 + 1)/(3*x) - 3*a**2*sqrt(-a**2*x**2 + 1)/(8*x**2) - a*sqrt(-a**2*x**2 + 1)/(3*x**3) - sqrt(-a**2*x**2 + 1)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_11():
    f = x**3*exp(2*atanh(a*x))
    F = -x**4/4 - 2*x**3/(3*a) - x**2/a**2 - 2*x/a**3 - 2*log(-a*x + 1)/a**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_12():
    f = x**2*exp(2*atanh(a*x))
    F = -x**3/3 - x**2/a - 2*x/a**2 - 2*log(-a*x + 1)/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_13():
    f = x*exp(2*atanh(a*x))
    F = -x**2/2 - 2*x/a - 2*log(-a*x + 1)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_14():
    f = exp(2*atanh(a*x))
    F = -x - 2*log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_15():
    f = exp(2*atanh(a*x))/x
    F = log(x) - 2*log(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_16():
    f = exp(2*atanh(a*x))/x**2
    F = 2*a*log(x) - 2*a*log(-a*x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_17():
    f = exp(2*atanh(a*x))/x**3
    F = 2*a**2*log(x) - 2*a**2*log(-a*x + 1) - 2*a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_18():
    f = exp(2*atanh(a*x))/x**4
    F = 2*a**3*log(x) - 2*a**3*log(-a*x + 1) - 2*a**2/x - a/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_19():
    f = x**2*exp(3*atanh(a*x))
    F = (a*x + 1)**3/(a**3*sqrt(-a**2*x**2 + 1)) + (a*x + 3)**2*sqrt(-a**2*x**2 + 1)/(3*a**3) + (3*a*x + 28)*sqrt(-a**2*x**2 + 1)/(6*a**3) - 11*asin(a*x)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_20():
    f = x*exp(3*atanh(a*x))
    F = 9*sqrt(-a**2*x**2 + 1)/(2*a**2) - 9*asin(a*x)/(2*a**2) + 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**2*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_21():
    f = exp(3*atanh(a*x))
    F = 2*(a*x + 1)**2/(a*sqrt(-a**2*x**2 + 1)) + 3*sqrt(-a**2*x**2 + 1)/a - 3*asin(a*x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_22():
    f = exp(3*atanh(a*x))/x
    F = -asin(a*x) - atanh(sqrt(-a**2*x**2 + 1)) + 4*sqrt(-a**2*x**2 + 1)/(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_23():
    f = exp(3*atanh(a*x))/x**2
    F = -3*a*atanh(sqrt(-a**2*x**2 + 1)) + 4*a*sqrt(-a**2*x**2 + 1)/(-a*x + 1) - sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_24():
    f = exp(3*atanh(a*x))/x**3
    F = -9*a**2*atanh(sqrt(-a**2*x**2 + 1))/2 + 4*a**2*sqrt(-a**2*x**2 + 1)/(-a*x + 1) - 3*a*sqrt(-a**2*x**2 + 1)/x - sqrt(-a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_25():
    f = exp(3*atanh(a*x))/x**4
    F = -11*a**3*atanh(sqrt(-a**2*x**2 + 1))/2 + 4*a**3*sqrt(-a**2*x**2 + 1)/(-a*x + 1) - 14*a**2*sqrt(-a**2*x**2 + 1)/(3*x) - 3*a*sqrt(-a**2*x**2 + 1)/(2*x**2) - sqrt(-a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_26():
    f = x**3*exp(4*atanh(a*x))
    F = x**4/4 + 4*x**3/(3*a) + 4*x**2/a**2 + 12*x/a**3 + 16*log(-a*x + 1)/a**4 + 4/(a**4*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_27():
    f = x**2*exp(4*atanh(a*x))
    F = x**3/3 + 2*x**2/a + 8*x/a**2 + 12*log(-a*x + 1)/a**3 + 4/(a**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_28():
    f = x*exp(4*atanh(a*x))
    F = x**2/2 + 4*x/a + 8*log(-a*x + 1)/a**2 + 4/(a**2*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_29():
    f = exp(4*atanh(a*x))
    F = x + 4*log(-a*x + 1)/a + 4/(a*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_30():
    f = exp(4*atanh(a*x))/x
    F = log(x) + 4/(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_31():
    f = exp(4*atanh(a*x))/x**2
    F = 4*a*log(x) - 4*a*log(-a*x + 1) + 4*a/(-a*x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_32():
    f = exp(4*atanh(a*x))/x**3
    F = 8*a**2*log(x) - 8*a**2*log(-a*x + 1) + 4*a**2/(-a*x + 1) - 4*a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_33():
    f = exp(4*atanh(a*x))/x**4
    F = 12*a**3*log(x) - 12*a**3*log(-a*x + 1) + 4*a**3/(-a*x + 1) - 8*a**2/x - 2*a/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_34():
    f = x**3*exp(-atanh(a*x))
    F = x**3*sqrt(-a**2*x**2 + 1)/(4*a) - x**2*sqrt(-a**2*x**2 + 1)/(3*a**2) - (-9*a*x + 16)*sqrt(-a**2*x**2 + 1)/(24*a**4) - 3*asin(a*x)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_35():
    f = x**2*exp(-atanh(a*x))
    F = -x*sqrt(-a**2*x**2 + 1)/(2*a**2) - (-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**3) + sqrt(-a**2*x**2 + 1)/a**3 + asin(a*x)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_36():
    f = x*exp(-atanh(a*x))
    F = -(-a*x + 2)*sqrt(-a**2*x**2 + 1)/(2*a**2) - asin(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_37():
    f = exp(-atanh(a*x))
    F = sqrt(-a**2*x**2 + 1)/a + asin(a*x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_38():
    f = exp(-atanh(a*x))/x
    F = -asin(a*x) - atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_39():
    f = exp(-atanh(a*x))/x**2
    F = a*atanh(sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_40():
    f = exp(-atanh(a*x))/x**3
    F = -a**2*atanh(sqrt(-a**2*x**2 + 1))/2 + a*sqrt(-a**2*x**2 + 1)/x - sqrt(-a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_41():
    f = exp(-atanh(a*x))/x**4
    F = a**3*atanh(sqrt(-a**2*x**2 + 1))/2 - 2*a**2*sqrt(-a**2*x**2 + 1)/(3*x) + a*sqrt(-a**2*x**2 + 1)/(2*x**2) - sqrt(-a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_42():
    f = exp(-atanh(a*x))/x**5
    F = -3*a**4*atanh(sqrt(-a**2*x**2 + 1))/8 + 2*a**3*sqrt(-a**2*x**2 + 1)/(3*x) - 3*a**2*sqrt(-a**2*x**2 + 1)/(8*x**2) + a*sqrt(-a**2*x**2 + 1)/(3*x**3) - sqrt(-a**2*x**2 + 1)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_43():
    f = x**3*exp(-2*atanh(a*x))
    F = -x**4/4 + 2*x**3/(3*a) - x**2/a**2 + 2*x/a**3 - 2*log(a*x + 1)/a**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_44():
    f = x**2*exp(-2*atanh(a*x))
    F = -x**3/3 + x**2/a - 2*x/a**2 + 2*log(a*x + 1)/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_45():
    f = x*exp(-2*atanh(a*x))
    F = -x**2/2 + 2*x/a - 2*log(a*x + 1)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_46():
    f = exp(-2*atanh(a*x))
    F = -x + 2*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_47():
    f = exp(-2*atanh(a*x))/x
    F = log(x) - 2*log(a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_48():
    f = exp(-2*atanh(a*x))/x**2
    F = -2*a*log(x) + 2*a*log(a*x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_49():
    f = exp(-2*atanh(a*x))/x**3
    F = 2*a**2*log(x) - 2*a**2*log(a*x + 1) + 2*a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_50():
    f = exp(-2*atanh(a*x))/x**4
    F = -2*a**3*log(x) + 2*a**3*log(a*x + 1) - 2*a**2/x + a/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_51():
    f = x**3*exp(-3*atanh(a*x))
    F = -x**3*sqrt(-a**2*x**2 + 1)/(4*a) + x**2*sqrt(-a**2*x**2 + 1)/a**2 + (-27*a*x + 18)*sqrt(-a**2*x**2 + 1)/(8*a**4) + (-a*x + 1)**3/(a**4*sqrt(-a**2*x**2 + 1)) + 27*sqrt(-a**2*x**2 + 1)/(4*a**4) + 51*asin(a*x)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_52():
    f = x**2*exp(-3*atanh(a*x))
    F = -(-3*a*x + 28)*sqrt(-a**2*x**2 + 1)/(6*a**3) - (-a*x + 1)**3/(a**3*sqrt(-a**2*x**2 + 1)) - (-a*x + 3)**2*sqrt(-a**2*x**2 + 1)/(3*a**3) - 11*asin(a*x)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_53():
    f = x*exp(-3*atanh(a*x))
    F = 9*sqrt(-a**2*x**2 + 1)/(2*a**2) + 9*asin(a*x)/(2*a**2) + 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**2*(a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**2*(a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_54():
    f = exp(-3*atanh(a*x))
    F = -2*(-a*x + 1)**2/(a*sqrt(-a**2*x**2 + 1)) - 3*sqrt(-a**2*x**2 + 1)/a - 3*asin(a*x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_55():
    f = exp(-3*atanh(a*x))/x
    F = asin(a*x) - atanh(sqrt(-a**2*x**2 + 1)) + 4*sqrt(-a**2*x**2 + 1)/(a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_56():
    f = exp(-3*atanh(a*x))/x**2
    F = 3*a*atanh(sqrt(-a**2*x**2 + 1)) - 4*a*sqrt(-a**2*x**2 + 1)/(a*x + 1) - sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_57():
    f = exp(-3*atanh(a*x))/x**3
    F = -9*a**2*atanh(sqrt(-a**2*x**2 + 1))/2 + 4*a**2*sqrt(-a**2*x**2 + 1)/(a*x + 1) + 3*a*sqrt(-a**2*x**2 + 1)/x - sqrt(-a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_58():
    f = exp(-3*atanh(a*x))/x**4
    F = 11*a**3*atanh(sqrt(-a**2*x**2 + 1))/2 - 4*a**3*sqrt(-a**2*x**2 + 1)/(a*x + 1) - 14*a**2*sqrt(-a**2*x**2 + 1)/(3*x) + 3*a*sqrt(-a**2*x**2 + 1)/(2*x**2) - sqrt(-a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_59():
    f = exp(-3*atanh(a*x))/x**5
    F = -51*a**4*atanh(sqrt(-a**2*x**2 + 1))/8 + 4*a**4*sqrt(-a**2*x**2 + 1)/(a*x + 1) + 6*a**3*sqrt(-a**2*x**2 + 1)/x - 19*a**2*sqrt(-a**2*x**2 + 1)/(8*x**2) + a*sqrt(-a**2*x**2 + 1)/x**3 - sqrt(-a**2*x**2 + 1)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_60():
    f = x**m*exp(atanh(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/4, sympy.S(1)/4, m + 2, -a*x, a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_61():
    f = x**2*exp(atanh(a*x)/2)
    F = -x*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)/(3*a**2) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)/(12*a**3) - 3*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(8*a**3) - 3*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) + 3*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_62():
    f = x*exp(atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)/(2*a**2) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*a**2) - sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) + sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_63():
    f = exp(atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/a - sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) + sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_64():
    f = exp(atanh(a*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 + sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 - 2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_65():
    f = exp(atanh(a*x)/2)/x**2
    F = -a*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - a*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_66():
    f = exp(atanh(a*x)/2)/x**3
    F = -a**2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - a**2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - a*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*x) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_67():
    f = exp(atanh(a*x)/2)/x**4
    F = -3*a**3*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 3*a**3*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 11*a**2*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(24*x) - 5*a*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(12*x**2) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_68():
    f = exp(atanh(a*x)/2)/x**5
    F = -11*a**4*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 11*a**4*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 83*a**3*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(192*x) - 29*a**2*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(96*x**2) - 7*a*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(24*x**3) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_69():
    f = exp(atanh(a*x)/2)/x**6
    F = -31*a**5*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/128 - 31*a**5*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/128 - 611*a**4*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(1920*x) - 269*a**3*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(960*x**2) - 11*a**2*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(48*x**3) - 9*a*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(40*x**4) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_70():
    f = x**m*exp(3*atanh(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-3)/4, sympy.S(3)/4, m + 2, -a*x, a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_71():
    f = x**3*exp(3*atanh(a*x)/2)
    F = -x**2*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(7)/4)/(4*a**2) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(7)/4)*(4*a*x + 11)/(32*a**4) - 41*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(64*a**4) + 123*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) - 123*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) - 123*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) - 123*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_72():
    f = x**2*exp(3*atanh(a*x)/2)
    F = -x*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(7)/4)/(3*a**2) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(7)/4)/(4*a**3) - 17*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(24*a**3) + 17*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) - 17*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) - 17*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 17*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_73():
    f = x*exp(3*atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(7)/4)/(2*a**2) - 3*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*a**2) + 9*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) - 9*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) - 9*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) - 9*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_74():
    f = exp(3*atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/a + 3*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) - 3*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_75():
    f = exp(3*atanh(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 + 2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_76():
    f = exp(3*atanh(a*x)/2)/x**2
    F = 3*a*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - 3*a*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_77():
    f = exp(3*atanh(a*x)/2)/x**3
    F = 9*a**2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - 9*a**2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - 3*a*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*x) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(7)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_78():
    f = exp(3*atanh(a*x)/2)/x**4
    F = 17*a**3*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 17*a**3*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 23*a**2*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(24*x) - 7*a*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(12*x**2) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_79():
    f = exp(3*atanh(a*x)/2)/x**5
    F = 123*a**4*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 123*a**4*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 63*a**3*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(64*x) - 15*a**2*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(32*x**2) - 3*a*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(8*x**3) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_80():
    f = x**m*exp(5*atanh(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-5)/4, sympy.S(5)/4, m + 2, -a*x, a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_81():
    f = x**3*exp(5*atanh(a*x)/2)
    F = 4*x**3*(a*x + 1)**(sympy.S(5)/4)/(a*(-a*x + 1)**(sympy.S(1)/4)) + 17*x**2*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)/(4*a**2) + (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)*(452*a*x + 521)/(96*a**4) + 475*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(64*a**4) + 475*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) - 475*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) + 475*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) + 475*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_82():
    f = x**2*exp(5*atanh(a*x)/2)
    F = (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(9)/4)/(3*a**3) + 11*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)/(4*a**3) + 55*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(8*a**3) + 55*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) - 55*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) + 55*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) + 55*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3) + 2*(a*x + 1)**(sympy.S(9)/4)/(a**3*(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_83():
    f = x*exp(5*atanh(a*x)/2)
    F = 5*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(5)/4)/(2*a**2) + 25*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*a**2) + 25*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) - 25*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) + 25*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) + 25*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2) + 2*(a*x + 1)**(sympy.S(9)/4)/(a**2*(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_84():
    f = exp(5*atanh(a*x)/2)
    F = 5*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/a + 5*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) - 5*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) + 5*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) + 5*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(2*a) + 4*(a*x + 1)**(sympy.S(5)/4)/(a*(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_85():
    f = exp(5*atanh(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 - 2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + 8*(a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_86():
    f = exp(5*atanh(a*x)/2)/x**2
    F = -5*a*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - 5*a*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + 10*a*(a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4) - (a*x + 1)**(sympy.S(5)/4)/(x*(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_87():
    f = exp(5*atanh(a*x)/2)/x**3
    F = -25*a**2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - 25*a**2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 + 25*a**2*(a*x + 1)**(sympy.S(1)/4)/(2*(-a*x + 1)**(sympy.S(1)/4)) - 5*a*(a*x + 1)**(sympy.S(5)/4)/(4*x*(-a*x + 1)**(sympy.S(1)/4)) - (a*x + 1)**(sympy.S(9)/4)/(2*x**2*(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_88():
    f = exp(5*atanh(a*x)/2)/x**4
    F = -55*a**3*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 55*a**3*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 + 287*a**3*(a*x + 1)**(sympy.S(1)/4)/(24*(-a*x + 1)**(sympy.S(1)/4)) - 61*a**2*(a*x + 1)**(sympy.S(1)/4)/(24*x*(-a*x + 1)**(sympy.S(1)/4)) - 13*a*(a*x + 1)**(sympy.S(1)/4)/(12*x**2*(-a*x + 1)**(sympy.S(1)/4)) - (a*x + 1)**(sympy.S(1)/4)/(3*x**3*(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_89():
    f = exp(5*atanh(a*x)/2)/x**5
    F = -475*a**4*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 475*a**4*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 + 2467*a**4*(a*x + 1)**(sympy.S(1)/4)/(192*(-a*x + 1)**(sympy.S(1)/4)) - 521*a**3*(a*x + 1)**(sympy.S(1)/4)/(192*x*(-a*x + 1)**(sympy.S(1)/4)) - 113*a**2*(a*x + 1)**(sympy.S(1)/4)/(96*x**2*(-a*x + 1)**(sympy.S(1)/4)) - 17*a*(a*x + 1)**(sympy.S(1)/4)/(24*x**3*(-a*x + 1)**(sympy.S(1)/4)) - (a*x + 1)**(sympy.S(1)/4)/(4*x**4*(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_90():
    f = x**m*exp(-atanh(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/4, sympy.S(1)/4, m + 2, a*x, -a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_91():
    f = x**3*exp(-atanh(a*x)/2)
    F = -x**2*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*a**2) - (-4*a*x + 25)*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(96*a**4) - 11*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(64*a**4) - 11*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) + 11*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) + 11*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) + 11*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_92():
    f = x**2*exp(-atanh(a*x)/2)
    F = -x*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(3*a**2) + (-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(12*a**3) + 3*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(8*a**3) + 3*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) - 3*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_93():
    f = x*exp(-atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(2*a**2) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*a**2) - sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) + sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_94():
    f = exp(-atanh(a*x)/2)
    F = (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/a + sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) - sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_95():
    f = exp(-atanh(a*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 + sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 + 2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_96():
    f = exp(-atanh(a*x)/2)/x**2
    F = -a*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + a*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_97():
    f = exp(-atanh(a*x)/2)/x**3
    F = a**2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - a**2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 + a*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*x) - (-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_98():
    f = exp(-atanh(a*x)/2)/x**4
    F = -3*a**3*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 + 3*a**3*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 11*a**2*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(24*x) + 5*a*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(12*x**2) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_99():
    f = exp(-atanh(a*x)/2)/x**5
    F = 11*a**4*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 11*a**4*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 + 83*a**3*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(192*x) - 29*a**2*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(96*x**2) + 7*a*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(24*x**3) - (-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_100():
    f = x**m*exp(-3*atanh(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-3)/4, sympy.S(3)/4, m + 2, a*x, -a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_101():
    f = x**3*exp(-3*atanh(a*x)/2)
    F = -x**2*(-a*x + 1)**(sympy.S(7)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*a**2) - (-4*a*x + 11)*(-a*x + 1)**(sympy.S(7)/4)*(a*x + 1)**(sympy.S(1)/4)/(32*a**4) - 41*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(64*a**4) + 123*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) - 123*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) + 123*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) + 123*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_102():
    f = x**2*exp(-3*atanh(a*x)/2)
    F = -x*(-a*x + 1)**(sympy.S(7)/4)*(a*x + 1)**(sympy.S(1)/4)/(3*a**2) + (-a*x + 1)**(sympy.S(7)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*a**3) + 17*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(24*a**3) - 17*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) + 17*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) - 17*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 17*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_103():
    f = x*exp(-3*atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(7)/4)*(a*x + 1)**(sympy.S(1)/4)/(2*a**2) - 3*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*a**2) + 9*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) - 9*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) + 9*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) + 9*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_104():
    f = exp(-3*atanh(a*x)/2)
    F = (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/a - 3*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) + 3*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_105():
    f = exp(-3*atanh(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 - 2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_106():
    f = exp(-3*atanh(a*x)/2)/x**2
    F = 3*a*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + 3*a*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_107():
    f = exp(-3*atanh(a*x)/2)/x**3
    F = -9*a**2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - 9*a**2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 + 3*a*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*x) - (-a*x + 1)**(sympy.S(7)/4)*(a*x + 1)**(sympy.S(1)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_108():
    f = exp(-3*atanh(a*x)/2)/x**4
    F = 17*a**3*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 + 17*a**3*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 23*a**2*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(24*x) + 7*a*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(12*x**2) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_109():
    f = exp(-3*atanh(a*x)/2)/x**5
    F = -123*a**4*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 123*a**4*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 + 63*a**3*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(64*x) - 15*a**2*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(32*x**2) + 3*a*(-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(8*x**3) - (-a*x + 1)**(sympy.S(3)/4)*(a*x + 1)**(sympy.S(1)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_110():
    f = x**m*exp(-5*atanh(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-5)/4, sympy.S(5)/4, m + 2, a*x, -a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_111():
    f = x**3*exp(-5*atanh(a*x)/2)
    F = -4*x**3*(-a*x + 1)**(sympy.S(5)/4)/(a*(a*x + 1)**(sympy.S(1)/4)) + 17*x**2*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*a**2) + (-452*a*x + 521)*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(96*a**4) + 475*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(64*a**4) + 475*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) - 475*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a**4) - 475*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) - 475*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_112():
    f = x**2*exp(-5*atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(9)/4)*(a*x + 1)**(sympy.S(3)/4)/(3*a**3) - 2*(-a*x + 1)**(sympy.S(9)/4)/(a**3*(a*x + 1)**(sympy.S(1)/4)) - 11*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*a**3) - 55*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(8*a**3) - 55*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) + 55*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(32*a**3) + 55*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) + 55*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_113():
    f = x*exp(-5*atanh(a*x)/2)
    F = 2*(-a*x + 1)**(sympy.S(9)/4)/(a**2*(a*x + 1)**(sympy.S(1)/4)) + 5*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(2*a**2) + 25*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*a**2) + 25*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) - 25*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a**2) - 25*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) - 25*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_114():
    f = exp(-5*atanh(a*x)/2)
    F = -4*(-a*x + 1)**(sympy.S(5)/4)/(a*(a*x + 1)**(sympy.S(1)/4)) - 5*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/a - 5*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) + 5*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(4*a) + 5*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) + 5*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_115():
    f = exp(-5*atanh(a*x)/2)/x
    F = 8*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/2 + 2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_116():
    f = exp(-5*atanh(a*x)/2)/x**2
    F = -10*a*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 5*a*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) + 5*a*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(5)/4)/(x*(a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_117():
    f = exp(-5*atanh(a*x)/2)/x**3
    F = 25*a**2*(-a*x + 1)**(sympy.S(1)/4)/(2*(a*x + 1)**(sympy.S(1)/4)) + 25*a**2*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 - 25*a**2*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/4 + 5*a*(-a*x + 1)**(sympy.S(5)/4)/(4*x*(a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(9)/4)/(2*x**2*(a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_118():
    f = exp(-5*atanh(a*x)/2)/x**4
    F = -287*a**3*(-a*x + 1)**(sympy.S(1)/4)/(24*(a*x + 1)**(sympy.S(1)/4)) - 55*a**3*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 + 55*a**3*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/8 - 61*a**2*(-a*x + 1)**(sympy.S(1)/4)/(24*x*(a*x + 1)**(sympy.S(1)/4)) + 13*a*(-a*x + 1)**(sympy.S(1)/4)/(12*x**2*(a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(1)/4)/(3*x**3*(a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_119():
    f = exp(-5*atanh(a*x)/2)/x**5
    F = 2467*a**4*(-a*x + 1)**(sympy.S(1)/4)/(192*(a*x + 1)**(sympy.S(1)/4)) + 475*a**4*atan((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 - 475*a**4*atanh((a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4))/64 + 521*a**3*(-a*x + 1)**(sympy.S(1)/4)/(192*x*(a*x + 1)**(sympy.S(1)/4)) - 113*a**2*(-a*x + 1)**(sympy.S(1)/4)/(96*x**2*(a*x + 1)**(sympy.S(1)/4)) + 17*a*(-a*x + 1)**(sympy.S(1)/4)/(24*x**3*(a*x + 1)**(sympy.S(1)/4)) - (-a*x + 1)**(sympy.S(1)/4)/(4*x**4*(a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_120():
    f = x**m*exp(atanh(x)/3)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/6, sympy.S(1)/6, m + 2, -x, x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_121():
    f = x**2*exp(atanh(x)/3)
    F = -x*(1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(7)/6)/3 - (1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(7)/6)/18 - 19*(1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(1)/6)/54 - 19*sqrt(3)*log(-sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/324 + 19*sqrt(3)*log(sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/324 - 19*atan((1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6))/81 - 19*atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) - sqrt(3))/162 - 19*atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + sqrt(3))/162
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_122():
    f = x*exp(atanh(x)/3)
    F = -(1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(7)/6)/2 - (1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(1)/6)/6 - sqrt(3)*log(-sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/36 + sqrt(3)*log(sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/36 - atan((1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6))/9 - atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) - sqrt(3))/18 - atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + sqrt(3))/18
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_123():
    f = exp(atanh(x)/3)
    F = -(1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(1)/6) - sqrt(3)*log(-sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/6 + sqrt(3)*log(sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/6 - 2*atan((1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6))/3 - atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) - sqrt(3))/3 - atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + sqrt(3))/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_124():
    f = exp(atanh(x)/3)/x
    F = log(1 + (x + 1)**(sympy.S(1)/3)/(1 - x)**(sympy.S(1)/3) - (x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/2 - log(1 + (x + 1)**(sympy.S(1)/3)/(1 - x)**(sympy.S(1)/3) + (x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/2 - sqrt(3)*log(-sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/2 + sqrt(3)*log(sqrt(3)*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + (1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*(x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/3) - sqrt(3)*atan(sqrt(3)*(1 + 2*(x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/3) - 2*atan((1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6)) - atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) - sqrt(3)) - atan(2*(1 - x)**(sympy.S(1)/6)/(x + 1)**(sympy.S(1)/6) + sqrt(3)) - 2*atanh((x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_125():
    f = exp(atanh(x)/3)/x**2
    F = log(1 + (x + 1)**(sympy.S(1)/3)/(1 - x)**(sympy.S(1)/3) - (x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/6 - log(1 + (x + 1)**(sympy.S(1)/3)/(1 - x)**(sympy.S(1)/3) + (x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/6 + sqrt(3)*atan(sqrt(3)*(1 - 2*(x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/3)/3 - sqrt(3)*atan(sqrt(3)*(1 + 2*(x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/3)/3 - 2*atanh((x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/3 - (1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(1)/6)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_126():
    f = exp(atanh(x)/3)/x**3
    F = log(1 + (x + 1)**(sympy.S(1)/3)/(1 - x)**(sympy.S(1)/3) - (x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/36 - log(1 + (x + 1)**(sympy.S(1)/3)/(1 - x)**(sympy.S(1)/3) + (x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/36 + sqrt(3)*atan(sqrt(3)*(1 - 2*(x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/3)/18 - sqrt(3)*atan(sqrt(3)*(1 + 2*(x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/3)/18 - atanh((x + 1)**(sympy.S(1)/6)/(1 - x)**(sympy.S(1)/6))/9 - (1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(1)/6)/(6*x) - (1 - x)**(sympy.S(5)/6)*(x + 1)**(sympy.S(7)/6)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_127():
    f = x**m*exp(2*atanh(x)/3)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/3, sympy.S(1)/3, m + 2, -x, x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_128():
    f = x**2*exp(2*atanh(x)/3)
    F = -x*(1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(4)/3)/3 - (1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(4)/3)/9 - 11*(1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(1)/3)/27 + 11*log(x + 1)/81 + 11*log((1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/27 - 22*sqrt(3)*atan(2*sqrt(3)*(1 - x)**(sympy.S(1)/3)/(3*(x + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/81
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_129():
    f = x*exp(2*atanh(x)/3)
    F = -(1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(4)/3)/2 - (1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(1)/3)/3 + log(x + 1)/9 + log((1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/3 - 2*sqrt(3)*atan(2*sqrt(3)*(1 - x)**(sympy.S(1)/3)/(3*(x + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/9
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_130():
    f = exp(2*atanh(x)/3)
    F = -(1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(1)/3) + log(x + 1)/3 + log((1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1) - 2*sqrt(3)*atan(2*sqrt(3)*(1 - x)**(sympy.S(1)/3)/(3*(x + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_131():
    f = exp(2*atanh(x)/3)/x
    F = -log(x)/2 + log(x + 1)/2 + 3*log((1 - x)**(sympy.S(1)/3)/(x + 1)**(sympy.S(1)/3) + 1)/2 + 3*log((1 - x)**(sympy.S(1)/3) - (x + 1)**(sympy.S(1)/3))/2 - sqrt(3)*atan(2*sqrt(3)*(1 - x)**(sympy.S(1)/3)/(3*(x + 1)**(sympy.S(1)/3)) - sqrt(3)/3) + sqrt(3)*atan(2*sqrt(3)*(1 - x)**(sympy.S(1)/3)/(3*(x + 1)**(sympy.S(1)/3)) + sqrt(3)/3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_132():
    f = exp(2*atanh(x)/3)/x**2
    F = -log(x)/3 + log((1 - x)**(sympy.S(1)/3) - (x + 1)**(sympy.S(1)/3)) + 2*sqrt(3)*atan(2*sqrt(3)*(1 - x)**(sympy.S(1)/3)/(3*(x + 1)**(sympy.S(1)/3)) + sqrt(3)/3)/3 - (1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(1)/3)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_133():
    f = exp(2*atanh(x)/3)/x**3
    F = -log(x)/9 + log((1 - x)**(sympy.S(1)/3) - (x + 1)**(sympy.S(1)/3))/3 + 2*sqrt(3)*atan(2*sqrt(3)*(1 - x)**(sympy.S(1)/3)/(3*(x + 1)**(sympy.S(1)/3)) + sqrt(3)/3)/9 - (1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(1)/3)/(3*x) - (1 - x)**(sympy.S(2)/3)*(x + 1)**(sympy.S(4)/3)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_134():
    f = x**m*exp(atanh(a*x)/4)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/8, sympy.S(1)/8, m + 2, -a*x, a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_135():
    f = x**2*exp(atanh(a*x)/4)
    F = -x*(-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(9)/8)/(3*a**2) - (-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(9)/8)/(24*a**3) - 11*(-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(1)/8)/(32*a**3) - 11*sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) + 11*sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) - 11*sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) + 11*sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) + 11*sqrt(2 - sqrt(2))*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(128*a**3) - 11*sqrt(2 - sqrt(2))*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(128*a**3) + 11*sqrt(sqrt(2) + 2)*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(128*a**3) - 11*sqrt(sqrt(2) + 2)*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(128*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_136():
    f = x*exp(atanh(a*x)/4)
    F = -(-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(9)/8)/(2*a**2) - (-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(1)/8)/(8*a**2) - sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) + sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) - sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) + sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) + sqrt(2 - sqrt(2))*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(32*a**2) - sqrt(2 - sqrt(2))*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(32*a**2) + sqrt(sqrt(2) + 2)*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(32*a**2) - sqrt(sqrt(2) + 2)*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(32*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_137():
    f = exp(atanh(a*x)/4)
    F = -(-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(1)/8)/a - sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) + sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) - sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) + sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) + sqrt(2 - sqrt(2))*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*a) - sqrt(2 - sqrt(2))*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*a) + sqrt(sqrt(2) + 2)*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*a) - sqrt(sqrt(2) + 2)*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_138():
    f = exp(atanh(a*x)/4)/x
    F = sqrt(2)*log(1 + (a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4) - sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/2 - sqrt(2)*log(1 + (a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4) + sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/2 - sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/2 + sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/2 - sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/2 + sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + (-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/2 + sqrt(2 - sqrt(2))*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) - sqrt(2 - sqrt(2))*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) + sqrt(sqrt(2) + 2)*atan((-2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2)) - sqrt(sqrt(2) + 2)*atan((2*(-a*x + 1)**(sympy.S(1)/8)/(a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2)) - 2*atan((a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8)) + sqrt(2)*atan(1 - sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8)) - sqrt(2)*atan(1 + sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8)) - 2*atanh((a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_139():
    f = exp(atanh(a*x)/4)/x**2
    F = sqrt(2)*a*log(1 + (a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4) - sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/8 - sqrt(2)*a*log(1 + (a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4) + sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/8 - a*atan((a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/2 + sqrt(2)*a*atan(1 - sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/4 - sqrt(2)*a*atan(1 + sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/4 - a*atanh((a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/2 - (-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(1)/8)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_140():
    f = exp(atanh(a*x)/4)/x**3
    F = sqrt(2)*a**2*log(1 + (a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4) - sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/64 - sqrt(2)*a**2*log(1 + (a*x + 1)**(sympy.S(1)/4)/(-a*x + 1)**(sympy.S(1)/4) + sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/64 - a**2*atan((a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/16 + sqrt(2)*a**2*atan(1 - sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/32 - sqrt(2)*a**2*atan(1 + sqrt(2)*(a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/32 - a**2*atanh((a*x + 1)**(sympy.S(1)/8)/(-a*x + 1)**(sympy.S(1)/8))/16 - a*(-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(1)/8)/(8*x) - (-a*x + 1)**(sympy.S(7)/8)*(a*x + 1)**(sympy.S(9)/8)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_141():
    f = x**m*exp(4*atanh(a*x))
    F = -4*x**(m + 1)*hyper((1, m + 1), (m + 2,), a*x) + 4*x**(m + 1)/(-a*x + 1) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_142():
    f = x**m*exp(3*atanh(a*x))
    F = -a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) + 4*a*x**(m + 2)*hyper((sympy.S(3)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) - 3*x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1) + 4*x**(m + 1)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_143():
    f = x**m*exp(2*atanh(a*x))
    F = 2*x**(m + 1)*hyper((1, m + 1), (m + 2,), a*x)/(m + 1) - x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_144():
    f = x**m*exp(atanh(a*x))
    F = a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) + x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_145():
    f = x**m*exp(-atanh(a*x))
    F = -a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) + x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_146():
    f = x**m*exp(-2*atanh(a*x))
    F = 2*x**(m + 1)*hyper((1, m + 1), (m + 2,), -a*x)/(m + 1) - x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_147():
    f = x**m*exp(-3*atanh(a*x))
    F = a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) - 4*a*x**(m + 2)*hyper((sympy.S(3)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) - 3*x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1) + 4*x**(m + 1)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_148():
    f = x**m*exp(n*atanh(a*x))
    F = x**(m + 1)*appellf1(m + 1, -n/2, n/2, m + 2, -a*x, a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_149():
    f = x**3*exp(n*atanh(a*x))
    F = -2**(n/2 - 2)*n*(n**2 + 8)*(-a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(3*a**4*(2 - n)) - x**2*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)/(4*a**2) - (-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)*(2*a*n*x + n**2 + 6)/(24*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_150():
    f = x**2*exp(n*atanh(a*x))
    F = -2**(n/2)*(n**2 + 2)*(-a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(3*a**3*(2 - n)) - x*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)/(3*a**2) - n*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_151():
    f = x*exp(n*atanh(a*x))
    F = -2**(n/2)*n*(-a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(a**2*(2 - n)) - (-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_152():
    f = exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*(-a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(a*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_153():
    f = exp(n*atanh(a*x))/x
    F = -2**(n/2 + 1)*hyper((-n/2, -n/2), (1 - n/2,), -a*x/2 + sympy.S.Half)/(n*(-a*x + 1)**(n/2)) + 2*(a*x + 1)**(n/2)*hyper((1, -n/2), (1 - n/2,), (-a*x + 1)/(a*x + 1))/(n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_154():
    f = exp(n*atanh(a*x))/x**2
    F = -4*a*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (-a*x + 1)/(a*x + 1))/(2 - n)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_155():
    f = exp(n*atanh(a*x))/x**3
    F = -2*a**2*n*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (-a*x + 1)/(a*x + 1))/(2 - n) - (-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_156():
    f = exp(n*atanh(a*x))/x**4
    F = -2*a**3*(n**2 + 2)*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (-a*x + 1)/(a*x + 1))/(6 - 3*n) - a*n*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)/(6*x**2) - (-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_157():
    f = (-a*c*x + c)**p*exp(atanh(a*x))
    F = -2*sqrt(2)*(-a*c*x + c)**(p + 1)*hyper((sympy.S(-1)/2, p + sympy.S.Half), (p + sympy.S(3)/2,), -a*x/2 + sympy.S.Half)/(a*c*(2*p + 1)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_158():
    f = (-a*c*x + c)**4*exp(atanh(a*x))
    F = 7*c**4*x*sqrt(-a**2*x**2 + 1)/8 + c**4*(-a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a) + 7*c**4*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(20*a) + 7*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(12*a) + 7*c**4*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_159():
    f = (-a*c*x + c)**3*exp(atanh(a*x))
    F = 5*c**3*x*sqrt(-a**2*x**2 + 1)/8 + c**3*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*a) + 5*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(12*a) + 5*c**3*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_160():
    f = (-a*c*x + c)**2*exp(atanh(a*x))
    F = c**2*x*sqrt(-a**2*x**2 + 1)/2 + c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a) + c**2*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_161():
    f = (-a*c*x + c)*exp(atanh(a*x))
    F = c*x*sqrt(-a**2*x**2 + 1)/2 + c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_162():
    f = exp(atanh(a*x))/(-a*c*x + c)
    F = -asin(a*x)/(a*c) + 2*sqrt(-a**2*x**2 + 1)/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_163():
    f = exp(atanh(a*x))/(-a*c*x + c)**2
    F = (-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_164():
    f = exp(atanh(a*x))/(-a*c*x + c)**3
    F = (-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a*c**3*(-a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_165():
    f = exp(atanh(a*x))/(-a*c*x + c)**4
    F = 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a*c**4*(-a*x + 1)**3) + 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a*c**4*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_166():
    f = exp(atanh(a*x))/(-a*c*x + c)**5
    F = 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(315*a*c**5*(-a*x + 1)**3) + 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a*c**5*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(21*a*c**5*(-a*x + 1)**5) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(9*a*c**5*(-a*x + 1)**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_167():
    f = (-a*c*x + c)**p*exp(2*atanh(a*x))
    F = -2*(-a*c*x + c)**p/(a*p) + (-a*c*x + c)**(p + 1)/(a*c*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_168():
    f = (-a*c*x + c)**5*exp(2*atanh(a*x))
    F = c**5*(-a*x + 1)**6/(6*a) - 2*c**5*(-a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_169():
    f = (-a*c*x + c)**4*exp(2*atanh(a*x))
    F = c**4*(-a*x + 1)**5/(5*a) - c**4*(-a*x + 1)**4/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_170():
    f = (-a*c*x + c)**3*exp(2*atanh(a*x))
    F = c**3*(-a*x + 1)**4/(4*a) - 2*c**3*(-a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_171():
    f = (-a*c*x + c)**2*exp(2*atanh(a*x))
    F = -a**2*c**2*x**3/3 + c**2*x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_172():
    f = exp(2*atanh(a*x))/(-a*c*x + c)
    F = log(-a*x + 1)/(a*c) + 2/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_173():
    f = exp(2*atanh(a*x))/(-a*c*x + c)**2
    F = x/(c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_174():
    f = exp(2*atanh(a*x))/(-a*c*x + c)**3
    F = -1/(2*a*c**3*(-a*x + 1)**2) + 2/(3*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_175():
    f = exp(2*atanh(a*x))/(-a*c*x + c)**4
    F = -1/(3*a*c**4*(-a*x + 1)**3) + 1/(2*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_176():
    f = (-a*c*x + c)**p*exp(3*atanh(a*x))
    F = 4*sqrt(2)*(-a*c*x + c)**(p + 1)*hyper((sympy.S(-3)/2, p + sympy.S(-1)/2), (p + sympy.S.Half,), -a*x/2 + sympy.S.Half)/(a*c*(1 - 2*p)*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_177():
    f = (-a*c*x + c)**4*exp(3*atanh(a*x))
    F = c**4*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/4 + 3*c**4*x*sqrt(-a**2*x**2 + 1)/8 + c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a) + 3*c**4*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_178():
    f = (-a*c*x + c)**3*exp(3*atanh(a*x))
    F = c**3*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/4 + 3*c**3*x*sqrt(-a**2*x**2 + 1)/8 + 3*c**3*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_179():
    f = (-a*c*x + c)**2*exp(3*atanh(a*x))
    F = c**2*x*sqrt(-a**2*x**2 + 1)/2 - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a) + c**2*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_180():
    f = (-a*c*x + c)*exp(3*atanh(a*x))
    F = -3*c*sqrt(-a**2*x**2 + 1)/(2*a) + 3*c*asin(a*x)/(2*a) - c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_181():
    f = exp(3*atanh(a*x))/(-a*c*x + c)
    F = asin(a*x)/(a*c) - 2*sqrt(-a**2*x**2 + 1)/(a*c*(-a*x + 1)) + 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*c*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_182():
    f = exp(3*atanh(a*x))/(-a*c*x + c)**2
    F = (-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a*c**2*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_183():
    f = exp(3*atanh(a*x))/(-a*c*x + c)**3
    F = (-a**2*x**2 + 1)**(sympy.S(5)/2)/(35*a*c**3*(-a*x + 1)**5) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(7*a*c**3*(-a*x + 1)**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_184():
    f = exp(3*atanh(a*x))/(-a*c*x + c)**4
    F = 2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(315*a*c**4*(-a*x + 1)**5) + 2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(63*a*c**4*(-a*x + 1)**6) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(9*a*c**4*(-a*x + 1)**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_185():
    f = exp(3*atanh(a*x))/(-a*c*x + c)**5
    F = 2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(1155*a*c**5*(-a*x + 1)**5) + 2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(231*a*c**5*(-a*x + 1)**6) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(33*a*c**5*(-a*x + 1)**7) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(11*a*c**5*(-a*x + 1)**8)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_186():
    f = (-a*c*x + c)**p*exp(4*atanh(a*x))
    F = 4*c*(-a*c*x + c)**(p - 1)/(a*(1 - p)) + 4*(-a*c*x + c)**p/(a*p) - (-a*c*x + c)**(p + 1)/(a*c*(p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_187():
    f = (-a*c*x + c)**5*exp(4*atanh(a*x))
    F = -c**5*(-a*x + 1)**6/(6*a) + 4*c**5*(-a*x + 1)**5/(5*a) - c**5*(-a*x + 1)**4/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_188():
    f = (-a*c*x + c)**4*exp(4*atanh(a*x))
    F = a**4*c**4*x**5/5 - 2*a**2*c**4*x**3/3 + c**4*x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_189():
    f = (-a*c*x + c)**3*exp(4*atanh(a*x))
    F = -c**3*(a*x + 1)**4/(4*a) + 2*c**3*(a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_190():
    f = (-a*c*x + c)**2*exp(4*atanh(a*x))
    F = c**2*(a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_191():
    f = (-a*c*x + c)*exp(4*atanh(a*x))
    F = -a*c*x**2/2 - 3*c*x - 4*c*log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_192():
    f = exp(4*atanh(a*x))/(-a*c*x + c)
    F = -log(-a*x + 1)/(a*c) - 4/(a*c*(-a*x + 1)) + 2/(a*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_193():
    f = exp(4*atanh(a*x))/(-a*c*x + c)**2
    F = (a*x + 1)**3/(6*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_194():
    f = exp(4*atanh(a*x))/(-a*c*x + c)**3
    F = 1/(2*a*c**3*(-a*x + 1)**2) - 4/(3*a*c**3*(-a*x + 1)**3) + 1/(a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_195():
    f = exp(4*atanh(a*x))/(-a*c*x + c)**4
    F = 1/(3*a*c**4*(-a*x + 1)**3) - 1/(a*c**4*(-a*x + 1)**4) + 4/(5*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_196():
    f = (-a*c*x + c)**p*exp(-atanh(a*x))
    F = -sqrt(2)*sqrt(-a*x + 1)*(-a*c*x + c)**(p + 1)*hyper((sympy.S.Half, p + sympy.S(3)/2), (p + sympy.S(5)/2,), -a*x/2 + sympy.S.Half)/(a*c*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_197():
    f = (-a*c*x + c)**3*exp(-atanh(a*x))
    F = c**3*(-a*x + 1)**3*sqrt(-a**2*x**2 + 1)/(4*a) + 7*c**3*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1)/(12*a) + 35*c**3*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(24*a) + 35*c**3*sqrt(-a**2*x**2 + 1)/(8*a) + 35*c**3*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_198():
    f = (-a*c*x + c)**2*exp(-atanh(a*x))
    F = c**2*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1)/(3*a) + 5*c**2*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(6*a) + 5*c**2*sqrt(-a**2*x**2 + 1)/(2*a) + 5*c**2*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_199():
    f = (-a*c*x + c)*exp(-atanh(a*x))
    F = c*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(2*a) + 3*c*sqrt(-a**2*x**2 + 1)/(2*a) + 3*c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_200():
    f = exp(-atanh(a*x))/(-a*c*x + c)
    F = asin(a*x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_201():
    f = exp(-atanh(a*x))/(-a*c*x + c)**2
    F = sqrt(-a**2*x**2 + 1)/(a*c**2*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_202():
    f = exp(-atanh(a*x))/(-a*c*x + c)**3
    F = sqrt(-a**2*x**2 + 1)/(3*a*c**3*(-a*x + 1)) + sqrt(-a**2*x**2 + 1)/(3*a*c**3*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_203():
    f = exp(-atanh(a*x))/(-a*c*x + c)**4
    F = 2*sqrt(-a**2*x**2 + 1)/(15*a*c**4*(-a*x + 1)) + 2*sqrt(-a**2*x**2 + 1)/(15*a*c**4*(-a*x + 1)**2) + sqrt(-a**2*x**2 + 1)/(5*a*c**4*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_204():
    f = exp(-atanh(a*x))/(-a*c*x + c)**5
    F = 2*sqrt(-a**2*x**2 + 1)/(35*a*c**5*(-a*x + 1)) + 2*sqrt(-a**2*x**2 + 1)/(35*a*c**5*(-a*x + 1)**2) + 3*sqrt(-a**2*x**2 + 1)/(35*a*c**5*(-a*x + 1)**3) + sqrt(-a**2*x**2 + 1)/(7*a*c**5*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_205():
    f = (-a*c*x + c)**p*exp(-2*atanh(a*x))
    F = -(-a*c*x + c)**(p + 2)*hyper((1, p + 2), (p + 3,), -a*x/2 + sympy.S.Half)/(2*a*c**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_206():
    f = (-a*c*x + c)**4*exp(-2*atanh(a*x))
    F = -16*c**4*x + c**4*(-a*x + 1)**5/(5*a) + c**4*(-a*x + 1)**4/(2*a) + 4*c**4*(-a*x + 1)**3/(3*a) + 4*c**4*(-a*x + 1)**2/a + 32*c**4*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_207():
    f = (-a*c*x + c)**3*exp(-2*atanh(a*x))
    F = -8*c**3*x + c**3*(-a*x + 1)**4/(4*a) + 2*c**3*(-a*x + 1)**3/(3*a) + 2*c**3*(-a*x + 1)**2/a + 16*c**3*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_208():
    f = (-a*c*x + c)**2*exp(-2*atanh(a*x))
    F = -4*c**2*x + c**2*(-a*x + 1)**3/(3*a) + c**2*(-a*x + 1)**2/a + 8*c**2*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_209():
    f = (-a*c*x + c)*exp(-2*atanh(a*x))
    F = a*c*x**2/2 - 3*c*x + 4*c*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_210():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)
    F = log(a*x + 1)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_211():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**2
    F = atanh(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_212():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**3
    F = atanh(a*x)/(2*a*c**3) + 1/(2*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_213():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**4
    F = atanh(a*x)/(4*a*c**4) + 1/(4*a*c**4*(-a*x + 1)) + 1/(4*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_214():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**5
    F = atanh(a*x)/(8*a*c**5) + 1/(8*a*c**5*(-a*x + 1)) + 1/(8*a*c**5*(-a*x + 1)**2) + 1/(6*a*c**5*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_215():
    f = (-a*c*x + c)**p*exp(-3*atanh(a*x))
    F = -sqrt(2)*(-a*x + 1)**(sympy.S(3)/2)*(-a*c*x + c)**(p + 1)*hyper((sympy.S(3)/2, p + sympy.S(5)/2), (p + sympy.S(7)/2,), -a*x/2 + sympy.S.Half)/(2*a*c*(2*p + 5))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_216():
    f = (-a*c*x + c)**3*exp(-3*atanh(a*x))
    F = -2*c**3*(-a*x + 1)**5/(a*sqrt(-a**2*x**2 + 1)) - 9*c**3*(-a*x + 1)**3*sqrt(-a**2*x**2 + 1)/(4*a) - 21*c**3*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1)/(4*a) - 105*c**3*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(8*a) - 315*c**3*sqrt(-a**2*x**2 + 1)/(8*a) - 315*c**3*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_217():
    f = (-a*c*x + c)**2*exp(-3*atanh(a*x))
    F = -2*c**2*(-a*x + 1)**4/(a*sqrt(-a**2*x**2 + 1)) - 7*c**2*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1)/(3*a) - 35*c**2*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(6*a) - 35*c**2*sqrt(-a**2*x**2 + 1)/(2*a) - 35*c**2*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_218():
    f = (-a*c*x + c)*exp(-3*atanh(a*x))
    F = -2*c*(-a*x + 1)**3/(a*sqrt(-a**2*x**2 + 1)) - 5*c*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(2*a) - 15*c*sqrt(-a**2*x**2 + 1)/(2*a) - 15*c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_219():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)
    F = -(-2*a*x + 2)/(a*c*sqrt(-a**2*x**2 + 1)) - asin(a*x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_220():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**2
    F = -(-a*x + 1)/(a*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_221():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**3
    F = x/(c**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_222():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**4
    F = 2*x/(3*c**4*sqrt(-a**2*x**2 + 1)) + 1/(3*a*c**4*(-a*x + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_223():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**5
    F = 2*x/(5*c**5*sqrt(-a**2*x**2 + 1)) + 1/(5*a*c**5*(-a*x + 1)*sqrt(-a**2*x**2 + 1)) + 1/(5*a*c**5*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_224():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**6
    F = 8*x/(35*c**6*sqrt(-a**2*x**2 + 1)) + 4/(35*a*c**6*(-a*x + 1)*sqrt(-a**2*x**2 + 1)) + 4/(35*a*c**6*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1)) + 1/(7*a*c**6*(-a*x + 1)**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_225():
    f = (-a*c*x + c)**(sympy.S(9)/2)*exp(atanh(a*x))
    F = 4096*c**6*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3465*a*(-a*c*x + c)**(sympy.S(3)/2)) + 1024*c**5*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(1155*a*sqrt(-a*c*x + c)) + 128*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(231*a) + 32*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(99*a) + 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(5)/2)/(11*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_226():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(atanh(a*x))
    F = 256*c**5*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(315*a*(-a*c*x + c)**(sympy.S(3)/2)) + 64*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a*sqrt(-a*c*x + c)) + 8*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(21*a) + 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(9*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_227():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(atanh(a*x))
    F = 64*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a*(-a*c*x + c)**(sympy.S(3)/2)) + 16*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a*sqrt(-a*c*x + c)) + 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)*sqrt(-a*c*x + c)/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_228():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(atanh(a*x))
    F = 8*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a*(-a*c*x + c)**(sympy.S(3)/2)) + 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_229():
    f = sqrt(-a*c*x + c)*exp(atanh(a*x))
    F = 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_230():
    f = exp(atanh(a*x))/sqrt(-a*c*x + c)
    F = -2*sqrt(-a**2*x**2 + 1)/(a*sqrt(-a*c*x + c)) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_231():
    f = exp(atanh(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = sqrt(-a**2*x**2 + 1)/(a*(-a*c*x + c)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(2*a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_232():
    f = exp(atanh(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = sqrt(-a**2*x**2 + 1)/(2*a*(-a*c*x + c)**(sympy.S(5)/2)) - sqrt(-a**2*x**2 + 1)/(8*a*c*(-a*c*x + c)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(16*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_233():
    f = exp(atanh(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = sqrt(-a**2*x**2 + 1)/(3*a*(-a*c*x + c)**(sympy.S(7)/2)) - sqrt(-a**2*x**2 + 1)/(24*a*c*(-a*c*x + c)**(sympy.S(5)/2)) - sqrt(-a**2*x**2 + 1)/(32*a*c**2*(-a*c*x + c)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(64*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_234():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(2*atanh(a*x))
    F = -4*(-a*c*x + c)**(sympy.S(7)/2)/(7*a) + 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_235():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = -4*(-a*c*x + c)**(sympy.S(5)/2)/(5*a) + 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_236():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = -4*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) + 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_237():
    f = sqrt(-a*c*x + c)*exp(2*atanh(a*x))
    F = -4*sqrt(-a*c*x + c)/a + 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_238():
    f = exp(2*atanh(a*x))/sqrt(-a*c*x + c)
    F = 4/(a*sqrt(-a*c*x + c)) + 2*sqrt(-a*c*x + c)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_239():
    f = exp(2*atanh(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = 4/(3*a*(-a*c*x + c)**(sympy.S(3)/2)) - 2/(a*c*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_240():
    f = exp(2*atanh(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = 4/(5*a*(-a*c*x + c)**(sympy.S(5)/2)) - 2/(3*a*c*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_241():
    f = exp(2*atanh(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = 4/(7*a*(-a*c*x + c)**(sympy.S(7)/2)) - 2/(5*a*c*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_242():
    f = (-a*c*x + c)**(sympy.S(9)/2)*exp(3*atanh(a*x))
    F = 256*c**7*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(1155*a*(-a*c*x + c)**(sympy.S(5)/2)) + 64*c**6*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(231*a*(-a*c*x + c)**(sympy.S(3)/2)) + 8*c**5*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(33*a*sqrt(-a*c*x + c)) + 2*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)*sqrt(-a*c*x + c)/(11*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_243():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(3*atanh(a*x))
    F = 64*c**6*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(315*a*(-a*c*x + c)**(sympy.S(5)/2)) + 16*c**5*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(63*a*(-a*c*x + c)**(sympy.S(3)/2)) + 2*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(9*a*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_244():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(3*atanh(a*x))
    F = 8*c**5*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(35*a*(-a*c*x + c)**(sympy.S(5)/2)) + 2*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(7*a*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_245():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(3*atanh(a*x))
    F = 2*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_246():
    f = sqrt(-a*c*x + c)*exp(3*atanh(a*x))
    F = 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/a - 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*(-a*c*x + c)**(sympy.S(3)/2)) - 4*c*sqrt(-a**2*x**2 + 1)/(a*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_247():
    f = exp(3*atanh(a*x))/sqrt(-a*c*x + c)
    F = c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(a*(-a*c*x + c)**(sympy.S(5)/2)) + 3*sqrt(-a**2*x**2 + 1)/(a*sqrt(-a*c*x + c)) - 3*sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_248():
    f = exp(3*atanh(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a*(-a*c*x + c)**(sympy.S(7)/2)) - 3*sqrt(-a**2*x**2 + 1)/(4*a*(-a*c*x + c)**(sympy.S(3)/2)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(8*a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_249():
    f = exp(3*atanh(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*(-a*c*x + c)**(sympy.S(9)/2)) - sqrt(-a**2*x**2 + 1)/(4*a*(-a*c*x + c)**(sympy.S(5)/2)) + sqrt(-a**2*x**2 + 1)/(16*a*c*(-a*c*x + c)**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(32*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_250():
    f = exp(3*atanh(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*a*(-a*c*x + c)**(sympy.S(11)/2)) - sqrt(-a**2*x**2 + 1)/(8*a*(-a*c*x + c)**(sympy.S(7)/2)) + sqrt(-a**2*x**2 + 1)/(64*a*c*(-a*c*x + c)**(sympy.S(5)/2)) + 3*sqrt(-a**2*x**2 + 1)/(256*a*c**2*(-a*c*x + c)**(sympy.S(3)/2)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(512*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_251():
    f = (-a*c*x + c)**(sympy.S(9)/2)*exp(-atanh(a*x))
    F = 16384*c**5*sqrt(-a**2*x**2 + 1)/(693*a*sqrt(-a*c*x + c)) + 4096*c**4*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(693*a) + 512*c**3*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(231*a) + 640*c**2*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(5)/2)/(693*a) + 40*c*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(7)/2)/(99*a) + 2*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(9)/2)/(11*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_252():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(-atanh(a*x))
    F = 4096*c**4*sqrt(-a**2*x**2 + 1)/(315*a*sqrt(-a*c*x + c)) + 1024*c**3*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(315*a) + 128*c**2*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(105*a) + 32*c*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(5)/2)/(63*a) + 2*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(7)/2)/(9*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_253():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(-atanh(a*x))
    F = 256*c**3*sqrt(-a**2*x**2 + 1)/(35*a*sqrt(-a*c*x + c)) + 64*c**2*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(35*a) + 24*c*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(35*a) + 2*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(5)/2)/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_254():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(-atanh(a*x))
    F = 64*c**2*sqrt(-a**2*x**2 + 1)/(15*a*sqrt(-a*c*x + c)) + 16*c*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(15*a) + 2*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_255():
    f = sqrt(-a*c*x + c)*exp(-atanh(a*x))
    F = 8*c*sqrt(-a**2*x**2 + 1)/(3*a*sqrt(-a*c*x + c)) + 2*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_256():
    f = exp(-atanh(a*x))/sqrt(-a*c*x + c)
    F = 2*sqrt(-a**2*x**2 + 1)/(a*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_257():
    f = exp(-atanh(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_258():
    f = exp(-atanh(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = sqrt(-a**2*x**2 + 1)/(2*a*c*(-a*c*x + c)**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(4*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_259():
    f = exp(-atanh(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = sqrt(-a**2*x**2 + 1)/(4*a*c*(-a*c*x + c)**(sympy.S(5)/2)) + 3*sqrt(-a**2*x**2 + 1)/(16*a*c**2*(-a*c*x + c)**(sympy.S(3)/2)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(32*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_260():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(-2*atanh(a*x))
    F = -32*sqrt(2)*c**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a + 32*c**3*sqrt(-a*c*x + c)/a + 16*c**2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) + 8*c*(-a*c*x + c)**(sympy.S(5)/2)/(5*a) + 4*(-a*c*x + c)**(sympy.S(7)/2)/(7*a) + 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_261():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(-2*atanh(a*x))
    F = -16*sqrt(2)*c**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a + 16*c**2*sqrt(-a*c*x + c)/a + 8*c*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) + 4*(-a*c*x + c)**(sympy.S(5)/2)/(5*a) + 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_262():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(-2*atanh(a*x))
    F = -8*sqrt(2)*c**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a + 8*c*sqrt(-a*c*x + c)/a + 4*(-a*c*x + c)**(sympy.S(3)/2)/(3*a) + 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_263():
    f = sqrt(-a*c*x + c)*exp(-2*atanh(a*x))
    F = -4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a + 4*sqrt(-a*c*x + c)/a + 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_264():
    f = exp(-2*atanh(a*x))/sqrt(-a*c*x + c)
    F = 2*sqrt(-a*c*x + c)/(a*c) - 2*sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_265():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_266():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = 1/(a*c**2*sqrt(-a*c*x + c)) - sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(2*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_267():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = 1/(3*a*c**2*(-a*c*x + c)**(sympy.S(3)/2)) + 1/(2*a*c**3*sqrt(-a*c*x + c)) - sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(4*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_268():
    f = exp(-2*atanh(a*x))/(-a*c*x + c)**(sympy.S(9)/2)
    F = 1/(5*a*c**2*(-a*c*x + c)**(sympy.S(5)/2)) + 1/(6*a*c**3*(-a*c*x + c)**(sympy.S(3)/2)) + 1/(4*a*c**4*sqrt(-a*c*x + c)) - sqrt(2)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/(8*a*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_269():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(-3*atanh(a*x))
    F = -4096*c**2*sqrt(-a*c*x + c)/(35*a*sqrt(-a**2*x**2 + 1)) + 1024*c*(-a*c*x + c)**(sympy.S(3)/2)/(35*a*sqrt(-a**2*x**2 + 1)) + 128*(-a*c*x + c)**(sympy.S(5)/2)/(35*a*sqrt(-a**2*x**2 + 1)) + 32*(-a*c*x + c)**(sympy.S(7)/2)/(35*a*c*sqrt(-a**2*x**2 + 1)) + 2*(-a*c*x + c)**(sympy.S(9)/2)/(7*a*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_270():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(-3*atanh(a*x))
    F = -256*c*sqrt(-a*c*x + c)/(5*a*sqrt(-a**2*x**2 + 1)) + 64*(-a*c*x + c)**(sympy.S(3)/2)/(5*a*sqrt(-a**2*x**2 + 1)) + 8*(-a*c*x + c)**(sympy.S(5)/2)/(5*a*c*sqrt(-a**2*x**2 + 1)) + 2*(-a*c*x + c)**(sympy.S(7)/2)/(5*a*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_271():
    f = sqrt(-a*c*x + c)*exp(-3*atanh(a*x))
    F = -64*sqrt(-a*c*x + c)/(3*a*sqrt(-a**2*x**2 + 1)) + 16*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c*sqrt(-a**2*x**2 + 1)) + 2*(-a*c*x + c)**(sympy.S(5)/2)/(3*a*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_272():
    f = exp(-3*atanh(a*x))/sqrt(-a*c*x + c)
    F = -8*sqrt(-a*c*x + c)/(a*c*sqrt(-a**2*x**2 + 1)) + 2*(-a*c*x + c)**(sympy.S(3)/2)/(a*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_273():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = -2*sqrt(-a*c*x + c)/(a*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_274():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = -sqrt(-a*c*x + c)/(a*c**3*sqrt(-a**2*x**2 + 1)) + sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(2*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_275():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = 1/(2*a*c**3*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)) - 3*sqrt(-a*c*x + c)/(4*a*c**4*sqrt(-a**2*x**2 + 1)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(8*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_276():
    f = exp(-3*atanh(a*x))/(-a*c*x + c)**(sympy.S(9)/2)
    F = 1/(4*a*c**3*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)) + 5/(16*a*c**4*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)) - 15*sqrt(-a*c*x + c)/(32*a*c**5*sqrt(-a**2*x**2 + 1)) + 15*sqrt(2)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/(64*a*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_277():
    f = (-a*c*x + c)**(sympy.S(7)/2)*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*(-a*c*x + c)**(sympy.S(9)/2)*hyper((-n/2, sympy.S(9)/2 - n/2), (sympy.S(11)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*(9 - n)*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_278():
    f = (-a*c*x + c)**(sympy.S(5)/2)*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*(-a*c*x + c)**(sympy.S(7)/2)*hyper((-n/2, sympy.S(7)/2 - n/2), (sympy.S(9)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*(7 - n)*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_279():
    f = (-a*c*x + c)**(sympy.S(3)/2)*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*(-a*c*x + c)**(sympy.S(5)/2)*hyper((-n/2, sympy.S(5)/2 - n/2), (sympy.S(7)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*(5 - n)*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_280():
    f = sqrt(-a*c*x + c)*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)*hyper((-n/2, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*(3 - n)*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_281():
    f = exp(n*atanh(a*x))/sqrt(-a*c*x + c)
    F = -2**(n/2 + 1)*sqrt(-a*c*x + c)*hyper((-n/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*(1 - n)*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_282():
    f = exp(n*atanh(a*x))/(-a*c*x + c)**(sympy.S(3)/2)
    F = 2**(n/2 + 1)*hyper((-n/2, -n/2 + sympy.S(-1)/2), (sympy.S.Half - n/2,), -a*x/2 + sympy.S.Half)/(a*c*(n + 1)*(-a*x + 1)**(n/2)*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_283():
    f = exp(n*atanh(a*x))/(-a*c*x + c)**(sympy.S(5)/2)
    F = 2**(n/2 + 1)*hyper((-n/2, -n/2 + sympy.S(-3)/2), (-n/2 + sympy.S(-1)/2,), -a*x/2 + sympy.S.Half)/(a*c*(n + 3)*(-a*x + 1)**(n/2)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_284():
    f = exp(n*atanh(a*x))/(-a*c*x + c)**(sympy.S(7)/2)
    F = 2**(n/2 + 1)*hyper((-n/2, -n/2 + sympy.S(-5)/2), (-n/2 + sympy.S(-3)/2,), -a*x/2 + sympy.S.Half)/(a*c*(n + 5)*(-a*x + 1)**(n/2)*(-a*c*x + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_285():
    f = x**4*(-a*c*x + c)*exp(atanh(a*x))
    F = c*x**5*sqrt(-a**2*x**2 + 1)/6 - c*x**3*sqrt(-a**2*x**2 + 1)/(24*a**2) - c*x*sqrt(-a**2*x**2 + 1)/(16*a**4) + c*asin(a*x)/(16*a**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_286():
    f = x**3*(-a*c*x + c)*exp(atanh(a*x))
    F = c*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a**4) - c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_287():
    f = x**2*(-a*c*x + c)*exp(atanh(a*x))
    F = c*x**3*sqrt(-a**2*x**2 + 1)/4 - c*x*sqrt(-a**2*x**2 + 1)/(8*a**2) + c*asin(a*x)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_288():
    f = x*(-a*c*x + c)*exp(atanh(a*x))
    F = -c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_289():
    f = (-a*c*x + c)*exp(atanh(a*x))
    F = c*x*sqrt(-a**2*x**2 + 1)/2 + c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_290():
    f = (-a*c*x + c)*exp(atanh(a*x))/x
    F = c*sqrt(-a**2*x**2 + 1) - c*atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_291():
    f = (-a*c*x + c)*exp(atanh(a*x))/x**2
    F = -a*c*asin(a*x) - c*sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_292():
    f = (-a*c*x + c)*exp(atanh(a*x))/x**3
    F = a**2*c*atanh(sqrt(-a**2*x**2 + 1))/2 - c*sqrt(-a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_293():
    f = (-a*c*x + c)*exp(atanh(a*x))/x**4
    F = -c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_294():
    f = x**3*(-a*c*x + c)**2*exp(atanh(a*x))
    F = c**2*x**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6*a) - c**2*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**2) - c**2*x*sqrt(-a**2*x**2 + 1)/(16*a**3) - c**2*(-15*a*x + 16)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(120*a**4) - c**2*asin(a*x)/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_295():
    f = x**2*(-a*c*x + c)**2*exp(atanh(a*x))
    F = -c**2*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*a**2) + c**2*x*sqrt(-a**2*x**2 + 1)/(8*a**2) - c**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a**3) + c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**3) + c**2*asin(a*x)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_296():
    f = x*(-a*c*x + c)**2*exp(atanh(a*x))
    F = -c**2*x*sqrt(-a**2*x**2 + 1)/(8*a) - c**2*(-3*a*x + 4)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(12*a**2) - c**2*asin(a*x)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_297():
    f = (-a*c*x + c)**2*exp(atanh(a*x))
    F = c**2*x*sqrt(-a**2*x**2 + 1)/2 + c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a) + c**2*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_298():
    f = (-a*c*x + c)**2*exp(atanh(a*x))/x
    F = c**2*(-a*x + 2)*sqrt(-a**2*x**2 + 1)/2 - c**2*asin(a*x)/2 - c**2*atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_299():
    f = (-a*c*x + c)**2*exp(atanh(a*x))/x**2
    F = -a*c**2*asin(a*x) + a*c**2*atanh(sqrt(-a**2*x**2 + 1)) - c**2*(a*x + 1)*sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_300():
    f = (-a*c*x + c)**2*exp(atanh(a*x))/x**3
    F = a**2*c**2*asin(a*x) + a**2*c**2*atanh(sqrt(-a**2*x**2 + 1))/2 - c**2*(-2*a*x + 1)*sqrt(-a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_301():
    f = (-a*c*x + c)**2*exp(atanh(a*x))/x**4
    F = -a**3*c**2*atanh(sqrt(-a**2*x**2 + 1))/2 + a*c**2*sqrt(-a**2*x**2 + 1)/(2*x**2) - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_302():
    f = (-a*c*x + c)**2*exp(atanh(a*x))/x**5
    F = a**4*c**2*atanh(sqrt(-a**2*x**2 + 1))/8 - a**2*c**2*sqrt(-a**2*x**2 + 1)/(8*x**2) + a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*x**3) - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_303():
    f = (-a*c*x + c)**2*exp(atanh(a*x))/x**6
    F = -a**5*c**2*atanh(sqrt(-a**2*x**2 + 1))/8 + a**3*c**2*sqrt(-a**2*x**2 + 1)/(8*x**2) - 2*a**2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*x**3) + a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*x**4) - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_304():
    f = (-a*c*x + c)**2*exp(atanh(a*x))/x**7
    F = a**6*c**2*atanh(sqrt(-a**2*x**2 + 1))/16 - a**4*c**2*sqrt(-a**2*x**2 + 1)/(16*x**2) + 2*a**3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*x**3) - a**2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(8*x**4) + a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*x**5) - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_305():
    f = x**3*(-a*c*x + c)**3*exp(atanh(a*x))
    F = -c**3*x**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/7 + c**3*x**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a) - 11*c**3*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a**2) - c**3*x*sqrt(-a**2*x**2 + 1)/(8*a**3) - c**3*(-105*a*x + 88)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(420*a**4) - c**3*asin(a*x)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_306():
    f = x**2*(-a*c*x + c)**3*exp(atanh(a*x))
    F = -c**3*x**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/6 + 2*c**3*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a) + 3*c**3*x*sqrt(-a**2*x**2 + 1)/(16*a**2) + c**3*(-45*a*x + 32)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(120*a**3) + 3*c**3*asin(a*x)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_307():
    f = x*(-a*c*x + c)**3*exp(atanh(a*x))
    F = -c**3*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/5 - c**3*x*sqrt(-a**2*x**2 + 1)/(4*a) - c**3*(-15*a*x + 14)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(30*a**2) - c**3*asin(a*x)/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_308():
    f = (-a*c*x + c)**3*exp(atanh(a*x))
    F = 5*c**3*x*sqrt(-a**2*x**2 + 1)/8 + c**3*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*a) + 5*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(12*a) + 5*c**3*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_309():
    f = (-a*c*x + c)**3*exp(atanh(a*x))/x
    F = c**3*(-a*x + 1)*sqrt(-a**2*x**2 + 1) - c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/3 - c**3*asin(a*x) - c**3*atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_310():
    f = (-a*c*x + c)**3*exp(atanh(a*x))/x**2
    F = -a*c**3*(a*x + 4)*sqrt(-a**2*x**2 + 1)/2 - a*c**3*asin(a*x)/2 + 2*a*c**3*atanh(sqrt(-a**2*x**2 + 1)) - c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_311():
    f = (-a*c*x + c)**3*exp(atanh(a*x))/x**3
    F = 2*a**2*c**3*asin(a*x) - a**2*c**3*atanh(sqrt(-a**2*x**2 + 1))/2 + a*c**3*(a*x + 4)*sqrt(-a**2*x**2 + 1)/(2*x) - c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_312():
    f = (-a*c*x + c)**3*exp(atanh(a*x))/x**4
    F = -a**3*c**3*asin(a*x) - a**3*c**3*atanh(sqrt(-a**2*x**2 + 1)) + a*c**3*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/x**2 - c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_313():
    f = (-a*c*x + c)**3*exp(atanh(a*x))/x**5
    F = 5*a**4*c**3*atanh(sqrt(-a**2*x**2 + 1))/8 - 5*a**2*c**3*sqrt(-a**2*x**2 + 1)/(8*x**2) + 2*a*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*x**3) - c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_314():
    f = (-a*c*x + c)**3*exp(atanh(a*x))/x**6
    F = -a**5*c**3*atanh(sqrt(-a**2*x**2 + 1))/4 + a**3*c**3*sqrt(-a**2*x**2 + 1)/(4*x**2) - 7*a**2*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*x**3) + a*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*x**4) - c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_315():
    f = x**3*(-a*c*x + c)**4*exp(atanh(a*x))
    F = a*c**4*x**5*(-a**2*x**2 + 1)**(sympy.S(3)/2)/8 - 3*c**4*x**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/7 + 29*c**4*x**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(48*a) - 19*c**4*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a**2) - 29*c**4*x*sqrt(-a**2*x**2 + 1)/(128*a**3) - c**4*(-3045*a*x + 2432)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6720*a**4) - 29*c**4*asin(a*x)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_316():
    f = x**2*(-a*c*x + c)**4*exp(atanh(a*x))
    F = a*c**4*x**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/7 - c**4*x**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/2 + 5*c**4*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a) + 5*c**4*x*sqrt(-a**2*x**2 + 1)/(16*a**2) + 5*c**4*(-21*a*x + 16)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(168*a**3) + 5*c**4*asin(a*x)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_317():
    f = x*(-a*c*x + c)**4*exp(atanh(a*x))
    F = -7*c**4*x*sqrt(-a**2*x**2 + 1)/(16*a) - c**4*(-a*x + 1)**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6*a**2) - c**4*(-a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(10*a**2) - 7*c**4*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(40*a**2) - 7*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(24*a**2) - 7*c**4*asin(a*x)/(16*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_318():
    f = (-a*c*x + c)**4*exp(atanh(a*x))
    F = 7*c**4*x*sqrt(-a**2*x**2 + 1)/8 + c**4*(-a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a) + 7*c**4*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(20*a) + 7*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(12*a) + 7*c**4*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_319():
    f = (-a*c*x + c)**4*exp(atanh(a*x))/x
    F = a*c**4*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/4 + c**4*(-13*a*x + 8)*sqrt(-a**2*x**2 + 1)/8 - c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2) - 13*c**4*asin(a*x)/8 - c**4*atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_320():
    f = (-a*c*x + c)**4*exp(atanh(a*x))/x**2
    F = -a*c**4*(-a*x + 6)*sqrt(-a**2*x**2 + 1)/2 + a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/3 + a*c**4*asin(a*x)/2 + 3*a*c**4*atanh(sqrt(-a**2*x**2 + 1)) - c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_321():
    f = (-a*c*x + c)**4*exp(atanh(a*x))/x**3
    F = 5*a**2*c**4*(a*x + 1)*sqrt(-a**2*x**2 + 1)/2 + 5*a**2*c**4*asin(a*x)/2 - 5*a**2*c**4*atanh(sqrt(-a**2*x**2 + 1))/2 + 3*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/x - c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_322():
    f = (-a*c*x + c)**4*exp(atanh(a*x))/x**4
    F = -3*a**3*c**4*asin(a*x) - a**3*c**4*atanh(sqrt(-a**2*x**2 + 1))/2 - a**2*c**4*(-a*x + 6)*sqrt(-a**2*x**2 + 1)/(2*x) + 3*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*x**2) - c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_323():
    f = (-a*c*x + c)**4*exp(atanh(a*x))/x**6
    F = -7*a**5*c**4*atanh(sqrt(-a**2*x**2 + 1))/8 + 7*a**3*c**4*sqrt(-a**2*x**2 + 1)/(8*x**2) - 17*a**2*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*x**3) + 3*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(4*x**4) - c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_324():
    f = (-a*c*x + c)**4*exp(atanh(a*x))/x**7
    F = 7*a**6*c**4*atanh(sqrt(-a**2*x**2 + 1))/16 - 7*a**4*c**4*sqrt(-a**2*x**2 + 1)/(16*x**2) + 11*a**3*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*x**3) - 7*a**2*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(8*x**4) + 3*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*x**5) - c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_325():
    f = x**4*exp(atanh(a*x))/(-a*c*x + c)
    F = x**3*sqrt(-a**2*x**2 + 1)/(4*a**2*c) + 2*x**2*sqrt(-a**2*x**2 + 1)/(3*a**3*c) + 11*x*sqrt(-a**2*x**2 + 1)/(8*a**4*c) + (a*x + 1)**2/(a**5*c*sqrt(-a**2*x**2 + 1)) + 13*sqrt(-a**2*x**2 + 1)/(3*a**5*c) - 27*asin(a*x)/(8*a**5*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_326():
    f = x**3*exp(atanh(a*x))/(-a*c*x + c)
    F = x**2*sqrt(-a**2*x**2 + 1)/(3*a**2*c) + x*sqrt(-a**2*x**2 + 1)/(a**3*c) + (a*x + 1)**2/(a**4*c*sqrt(-a**2*x**2 + 1)) + 11*sqrt(-a**2*x**2 + 1)/(3*a**4*c) - 3*asin(a*x)/(a**4*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_327():
    f = x**2*exp(atanh(a*x))/(-a*c*x + c)
    F = (a*x + 1)**2/(a**3*c*sqrt(-a**2*x**2 + 1)) + (a*x + 6)*sqrt(-a**2*x**2 + 1)/(2*a**3*c) - 5*asin(a*x)/(2*a**3*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_328():
    f = x*exp(atanh(a*x))/(-a*c*x + c)
    F = 2*sqrt(-a**2*x**2 + 1)/(a**2*c) - 2*asin(a*x)/(a**2*c) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**2*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_329():
    f = exp(atanh(a*x))/(-a*c*x + c)
    F = -asin(a*x)/(a*c) + 2*sqrt(-a**2*x**2 + 1)/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_330():
    f = exp(atanh(a*x))/(x*(-a*c*x + c))
    F = (2*a*x + 2)/(c*sqrt(-a**2*x**2 + 1)) - atanh(sqrt(-a**2*x**2 + 1))/c
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_331():
    f = exp(atanh(a*x))/(x**2*(-a*c*x + c))
    F = 2*a*(a*x + 1)/(c*sqrt(-a**2*x**2 + 1)) - 2*a*atanh(sqrt(-a**2*x**2 + 1))/c - sqrt(-a**2*x**2 + 1)/(c*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_332():
    f = exp(atanh(a*x))/(x**3*(-a*c*x + c))
    F = 2*a**2*(a*x + 1)/(c*sqrt(-a**2*x**2 + 1)) - 5*a**2*atanh(sqrt(-a**2*x**2 + 1))/(2*c) - 2*a*sqrt(-a**2*x**2 + 1)/(c*x) - sqrt(-a**2*x**2 + 1)/(2*c*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_333():
    f = exp(atanh(a*x))/(x**4*(-a*c*x + c))
    F = 2*a**3*(a*x + 1)/(c*sqrt(-a**2*x**2 + 1)) - 3*a**3*atanh(sqrt(-a**2*x**2 + 1))/c - 8*a**2*sqrt(-a**2*x**2 + 1)/(3*c*x) - a*sqrt(-a**2*x**2 + 1)/(c*x**2) - sqrt(-a**2*x**2 + 1)/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_334():
    f = x**4*exp(atanh(a*x))/(-a*c*x + c)**2
    F = -2*(a*x + 1)**3/(a**5*c**2*sqrt(-a**2*x**2 + 1)) + (a*x + 1)**3/(3*a**5*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - (a*x + 5)**2*sqrt(-a**2*x**2 + 1)/(3*a**5*c**2) - (a*x + 5)*sqrt(-a**2*x**2 + 1)/(6*a**5*c**2) - 5*sqrt(-a**2*x**2 + 1)/(2*a**5*c**2) + 17*asin(a*x)/(2*a**5*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_335():
    f = x**3*exp(atanh(a*x))/(-a*c*x + c)**2
    F = (a*x + 1)**3/(3*a**4*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 3*(a*x + 1)**2/(a**4*c**2*sqrt(-a**2*x**2 + 1)) - (a*x + 12)*sqrt(-a**2*x**2 + 1)/(2*a**4*c**2) + 11*asin(a*x)/(2*a**4*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_336():
    f = x**2*exp(atanh(a*x))/(-a*c*x + c)**2
    F = 3*asin(a*x)/(a**3*c**2) - 6*sqrt(-a**2*x**2 + 1)/(a**3*c**2*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**3*c**2*(-a*x + 1)**2) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**3*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_337():
    f = x*exp(atanh(a*x))/(-a*c*x + c)**2
    F = asin(a*x)/(a**2*c**2) - 2*sqrt(-a**2*x**2 + 1)/(a**2*c**2*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**2*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_338():
    f = exp(atanh(a*x))/(-a*c*x + c)**2
    F = (-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_339():
    f = exp(atanh(a*x))/(x*(-a*c*x + c)**2)
    F = (4*a*x + 4)/(3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (5*a*x + 3)/(3*c**2*sqrt(-a**2*x**2 + 1)) - atanh(sqrt(-a**2*x**2 + 1))/c**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_340():
    f = exp(atanh(a*x))/(x**2*(-a*c*x + c)**2)
    F = 4*a*(a*x + 1)/(3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a*(11*a*x + 9)/(3*c**2*sqrt(-a**2*x**2 + 1)) - 3*a*atanh(sqrt(-a**2*x**2 + 1))/c**2 - sqrt(-a**2*x**2 + 1)/(c**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_341():
    f = exp(atanh(a*x))/(x**3*(-a*c*x + c)**2)
    F = 4*a**2*(a*x + 1)/(3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a**2*(17*a*x + 15)/(3*c**2*sqrt(-a**2*x**2 + 1)) - 11*a**2*atanh(sqrt(-a**2*x**2 + 1))/(2*c**2) - 3*a*sqrt(-a**2*x**2 + 1)/(c**2*x) - sqrt(-a**2*x**2 + 1)/(2*c**2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_342():
    f = exp(atanh(a*x))/(x**4*(-a*c*x + c)**2)
    F = 4*a**3*(a*x + 1)/(3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a**3*(23*a*x + 21)/(3*c**2*sqrt(-a**2*x**2 + 1)) - 17*a**3*atanh(sqrt(-a**2*x**2 + 1))/(2*c**2) - 17*a**2*sqrt(-a**2*x**2 + 1)/(3*c**2*x) - 3*a*sqrt(-a**2*x**2 + 1)/(2*c**2*x**2) - sqrt(-a**2*x**2 + 1)/(3*c**2*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_343():
    f = x**4*exp(atanh(a*x))/(-a*c*x + c)**3
    F = (a*x + 1)**4/(5*a**5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 19*(a*x + 1)**3/(15*a**5*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 6*(a*x + 1)**2/(a**5*c**3*sqrt(-a**2*x**2 + 1)) + (a*x + 20)*sqrt(-a**2*x**2 + 1)/(2*a**5*c**3) - 19*asin(a*x)/(2*a**5*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_344():
    f = x**3*exp(atanh(a*x))/(-a*c*x + c)**3
    F = -4*asin(a*x)/(a**4*c**3) + 8*sqrt(-a**2*x**2 + 1)/(a**4*c**3*(-a*x + 1)) - (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**4*c**3*(-a*x + 1)**2) - 14*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a**4*c**3*(-a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**4*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_345():
    f = x**2*exp(atanh(a*x))/(-a*c*x + c)**3
    F = -asin(a*x)/(a**3*c**3) + 2*sqrt(-a**2*x**2 + 1)/(a**3*c**3*(-a*x + 1)) - 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**3*c**3*(-a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**3*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_346():
    f = x*exp(atanh(a*x))/(-a*c*x + c)**3
    F = -4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a**2*c**3*(-a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**2*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_347():
    f = exp(atanh(a*x))/(-a*c*x + c)**3
    F = (-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a*c**3*(-a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_348():
    f = exp(atanh(a*x))/(x*(-a*c*x + c)**3)
    F = 4*a*x/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (8*a*x + 5)/(5*c**3*sqrt(-a**2*x**2 + 1)) + (8*a*x + 8)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - atanh(sqrt(-a**2*x**2 + 1))/c**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_349():
    f = exp(atanh(a*x))/(x**2*(-a*c*x + c)**3)
    F = 8*a*(a*x + 1)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 4*a*(8*a*x + 5)/(15*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a*(79*a*x + 60)/(15*c**3*sqrt(-a**2*x**2 + 1)) - 4*a*atanh(sqrt(-a**2*x**2 + 1))/c**3 - sqrt(-a**2*x**2 + 1)/(c**3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_350():
    f = exp(atanh(a*x))/(x**3*(-a*c*x + c)**3)
    F = 8*a**2*(a*x + 1)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 4*a**2*(13*a*x + 10)/(15*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a**2*(164*a*x + 135)/(15*c**3*sqrt(-a**2*x**2 + 1)) - 19*a**2*atanh(sqrt(-a**2*x**2 + 1))/(2*c**3) - 4*a*sqrt(-a**2*x**2 + 1)/(c**3*x) - sqrt(-a**2*x**2 + 1)/(2*c**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_351():
    f = exp(atanh(a*x))/(x**4*(-a*c*x + c)**3)
    F = 8*a**3*(a*x + 1)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 4*a**3*(6*a*x + 5)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a**3*(93*a*x + 80)/(5*c**3*sqrt(-a**2*x**2 + 1)) - 18*a**3*atanh(sqrt(-a**2*x**2 + 1))/c**3 - 29*a**2*sqrt(-a**2*x**2 + 1)/(3*c**3*x) - 2*a*sqrt(-a**2*x**2 + 1)/(c**3*x**2) - sqrt(-a**2*x**2 + 1)/(3*c**3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_352():
    f = x**5*exp(atanh(a*x))/(-a*c*x + c)**4
    F = (a*x + 1)**5/(7*a**6*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - 33*(a*x + 1)**4/(35*a**6*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 317*(a*x + 1)**3/(105*a**6*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 10*(a*x + 1)**2/(a**6*c**4*sqrt(-a**2*x**2 + 1)) - (a*x + 30)*sqrt(-a**2*x**2 + 1)/(2*a**6*c**4) + 29*asin(a*x)/(2*a**6*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_353():
    f = x**4*exp(atanh(a*x))/(-a*c*x + c)**4
    F = 5*asin(a*x)/(a**5*c**4) - 10*sqrt(-a**2*x**2 + 1)/(a**5*c**4*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**5*c**4*(-a*x + 1)**2) + 184*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a**5*c**4*(-a*x + 1)**3) - 26*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a**5*c**4*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a**5*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_354():
    f = x**3*exp(atanh(a*x))/(-a*c*x + c)**4
    F = asin(a*x)/(a**4*c**4) - 2*sqrt(-a**2*x**2 + 1)/(a**4*c**4*(-a*x + 1)) + 86*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a**4*c**4*(-a*x + 1)**3) - 19*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a**4*c**4*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a**4*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_355():
    f = x**2*exp(atanh(a*x))/(-a*c*x + c)**4
    F = 23*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a**3*c**4*(-a*x + 1)**3) - 12*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a**3*c**4*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a**3*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_356():
    f = x*exp(atanh(a*x))/(-a*c*x + c)**4
    F = -(-a**2*x**2 + 1)**(sympy.S(3)/2)/(21*a**2*c**4*(-a*x + 1)**3) - (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a**2*c**4*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a**2*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_357():
    f = exp(atanh(a*x))/(-a*c*x + c)**4
    F = 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a*c**4*(-a*x + 1)**3) + 2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a*c**4*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_358():
    f = exp(atanh(a*x))/(x*(-a*c*x + c)**4)
    F = -(-12*a*x + 28)/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + (16*a*x + 16)/(7*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + (83*a*x + 35)/(105*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (166*a*x + 105)/(105*c**4*sqrt(-a**2*x**2 + 1)) - atanh(sqrt(-a**2*x**2 + 1))/c**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_359():
    f = exp(atanh(a*x))/(x**2*(-a*c*x + c)**4)
    F = 16*a*(a*x + 1)/(7*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + 4*a*(17*a*x + 7)/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + a*(307*a*x + 175)/(105*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a*(719*a*x + 525)/(105*c**4*sqrt(-a**2*x**2 + 1)) - 5*a*atanh(sqrt(-a**2*x**2 + 1))/c**4 - sqrt(-a**2*x**2 + 1)/(c**4*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_360():
    f = exp(atanh(a*x))/(x**3*(-a*c*x + c)**4)
    F = 16*a**2*(a*x + 1)/(7*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + 4*a**2*(31*a*x + 21)/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + a**2*(671*a*x + 455)/(105*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + a**2*(1867*a*x + 1470)/(105*c**4*sqrt(-a**2*x**2 + 1)) - 29*a**2*atanh(sqrt(-a**2*x**2 + 1))/(2*c**4) - 5*a*sqrt(-a**2*x**2 + 1)/(c**4*x) - sqrt(-a**2*x**2 + 1)/(2*c**4*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_361():
    f = x*(x + 1)*exp(atanh(x))
    F = -sqrt(1 - x)*(x + 1)**(sympy.S(5)/2)/3 - sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)/3 - sqrt(1 - x)*sqrt(x + 1) + asin(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_362():
    f = (x + 1)*exp(atanh(x))
    F = -sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)/2 - 3*sqrt(1 - x)*sqrt(x + 1)/2 + 3*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_363():
    f = x*(x + 1)**2*exp(atanh(x))
    F = -sqrt(1 - x)*(x + 1)**(sympy.S(7)/2)/4 - sqrt(1 - x)*(x + 1)**(sympy.S(5)/2)/4 - 5*sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)/8 - 15*sqrt(1 - x)*sqrt(x + 1)/8 + 15*asin(x)/8
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_364():
    f = (x + 1)**2*exp(atanh(x))
    F = -sqrt(1 - x)*(x + 1)**(sympy.S(5)/2)/3 - 5*sqrt(1 - x)*(x + 1)**(sympy.S(3)/2)/6 - 5*sqrt(1 - x)*sqrt(x + 1)/2 + 5*asin(x)/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_365():
    f = x*exp(atanh(x))/(x + 1)
    F = -sqrt(1 - x)*sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_366():
    f = exp(atanh(x))/(x + 1)
    F = asin(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_367():
    f = x*exp(atanh(x))/(x + 1)**2
    F = sqrt(1 - x)/sqrt(x + 1) + asin(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_368():
    f = exp(atanh(x))/(x + 1)**2
    F = -sqrt(1 - x)/sqrt(x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_369():
    f = x*(x + 1)**(sympy.S(3)/2)*exp(atanh(x))
    F = 2*(1 - x)**(sympy.S(7)/2)/7 - 2*(1 - x)**(sympy.S(5)/2) + 16*(1 - x)**(sympy.S(3)/2)/3 - 8*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_370():
    f = (x + 1)**(sympy.S(3)/2)*exp(atanh(x))
    F = -2*(1 - x)**(sympy.S(5)/2)/5 + 8*(1 - x)**(sympy.S(3)/2)/3 - 8*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_371():
    f = (1 - x)**(sympy.S(3)/2)*exp(atanh(x))
    F = -2*(x + 1)**(sympy.S(5)/2)/5 + 4*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_372():
    f = x*sqrt(x + 1)*exp(atanh(x))
    F = -2*(1 - x)**(sympy.S(5)/2)/5 + 2*(1 - x)**(sympy.S(3)/2) - 4*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_373():
    f = sqrt(x + 1)*exp(atanh(x))
    F = 2*(1 - x)**(sympy.S(3)/2)/3 - 4*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_374():
    f = x*sqrt(1 - x)*exp(atanh(x))
    F = 2*(x + 1)**(sympy.S(5)/2)/5 - 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_375():
    f = sqrt(1 - x)*exp(atanh(x))
    F = 2*(x + 1)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_376():
    f = x*exp(atanh(x))/sqrt(x + 1)
    F = 2*(1 - x)**(sympy.S(3)/2)/3 - 2*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_377():
    f = exp(atanh(x))/sqrt(x + 1)
    F = -2*sqrt(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_378():
    f = x*exp(atanh(x))/sqrt(1 - x)
    F = -2*(x + 1)**(sympy.S(3)/2)/3 - 2*sqrt(x + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(x + 1)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_379():
    f = exp(atanh(x))/sqrt(1 - x)
    F = -2*sqrt(x + 1) + 2*sqrt(2)*atanh(sqrt(2)*sqrt(x + 1)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_380():
    f = x*exp(atanh(x))/(x + 1)**(sympy.S(3)/2)
    F = -2*sqrt(1 - x) + sqrt(2)*atanh(sqrt(2)*sqrt(1 - x)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_381():
    f = exp(atanh(x))/(x + 1)**(sympy.S(3)/2)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(1 - x)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_382():
    f = x*exp(atanh(x))/(1 - x)**(sympy.S(3)/2)
    F = 5*sqrt(x + 1)/2 - 5*sqrt(2)*atanh(sqrt(2)*sqrt(x + 1)/2)/2 + (x + 1)**(sympy.S(3)/2)/(2 - 2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_383():
    f = exp(atanh(x))/(1 - x)**(sympy.S(3)/2)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(x + 1)/2)/2 + sqrt(x + 1)/(1 - x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_384():
    f = x**m*sqrt(-a*c*x + c)*exp(atanh(a*x))
    F = 2*c*x**m*(a*x + 1)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S(3)/2, -m), (sympy.S(5)/2,), a*x + 1)/(3*a*(-a*x)**m*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_385():
    f = x**2*sqrt(-a*c*x + c)*exp(atanh(a*x))
    F = 2*c**2*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a*(-a*c*x + c)**(sympy.S(3)/2)) - 8*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a**3*(-a*c*x + c)**(sympy.S(3)/2)) + 8*c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a**3*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_386():
    f = x*sqrt(-a*c*x + c)*exp(atanh(a*x))
    F = 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a**2*(-a*c*x + c)**(sympy.S(3)/2)) - 2*c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a**2*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_387():
    f = sqrt(-a*c*x + c)*exp(atanh(a*x))
    F = 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_388():
    f = sqrt(-a*c*x + c)*exp(atanh(a*x))/x
    F = -2*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c)) + 2*c*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_389():
    f = sqrt(-a*c*x + c)*exp(atanh(a*x))/x**2
    F = -a*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c)) - c*sqrt(-a**2*x**2 + 1)/(x*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_390():
    f = x**3*sqrt(-a*c*x + c)*exp(2*atanh(a*x))
    F = -4*sqrt(-a*c*x + c)/a**4 + 14*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**4*c) - 18*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**4*c**2) + 10*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**4*c**3) - 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a**4*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_391():
    f = x**2*sqrt(-a*c*x + c)*exp(2*atanh(a*x))
    F = -4*sqrt(-a*c*x + c)/a**3 + 10*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**3*c) - 8*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**3*c**2) + 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**3*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_392():
    f = x*sqrt(-a*c*x + c)*exp(2*atanh(a*x))
    F = -4*sqrt(-a*c*x + c)/a**2 + 2*(-a*c*x + c)**(sympy.S(3)/2)/(a**2*c) - 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_393():
    f = sqrt(-a*c*x + c)*exp(2*atanh(a*x))
    F = -4*sqrt(-a*c*x + c)/a + 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_394():
    f = sqrt(-a*c*x + c)*exp(2*atanh(a*x))/x
    F = -2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) - 2*sqrt(-a*c*x + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_395():
    f = sqrt(-a*c*x + c)*exp(2*atanh(a*x))/x**2
    F = -3*a*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) - sqrt(-a*c*x + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_396():
    f = sqrt(-a*c*x + c)*exp(2*atanh(a*x))/x**3
    F = -7*a**2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/4 - 7*a*sqrt(-a*c*x + c)/(4*x) - sqrt(-a*c*x + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_397():
    f = sqrt(-a*c*x + c)*exp(2*atanh(a*x))/x**4
    F = -11*a**3*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/8 - 11*a**2*sqrt(-a*c*x + c)/(8*x) - 11*a*sqrt(-a*c*x + c)/(12*x**2) - sqrt(-a*c*x + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_398():
    f = sqrt(-a*c*x + c)*exp(2*atanh(a*x))/x**5
    F = -75*a**4*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/64 - 75*a**3*sqrt(-a*c*x + c)/(64*x) - 25*a**2*sqrt(-a*c*x + c)/(32*x**2) - 5*a*sqrt(-a*c*x + c)/(8*x**3) - sqrt(-a*c*x + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_399():
    f = x**3*sqrt(-a*c*x + c)*exp(3*atanh(a*x))
    F = -2*(a*x + 1)**(sympy.S(9)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(9*a**4*c*(-a*x + 1)**(sympy.S(3)/2)) + 2*(a*x + 1)**(sympy.S(7)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(7*a**4*c*(-a*x + 1)**(sympy.S(3)/2)) - 2*(a*x + 1)**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(5*a**4*c*(-a*x + 1)**(sympy.S(3)/2)) - 2*(a*x + 1)**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**4*c*(-a*x + 1)**(sympy.S(3)/2)) - 4*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(a**4*c*(-a*x + 1)**(sympy.S(3)/2)) + 4*sqrt(2)*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(a**4*c*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_400():
    f = x**2*sqrt(-a*c*x + c)*exp(3*atanh(a*x))
    F = -2*(a*x + 1)**(sympy.S(7)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(7*a**3*c*(-a*x + 1)**(sympy.S(3)/2)) - 2*(a*x + 1)**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**3*c*(-a*x + 1)**(sympy.S(3)/2)) - 4*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(a**3*c*(-a*x + 1)**(sympy.S(3)/2)) + 4*sqrt(2)*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(a**3*c*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_401():
    f = x*sqrt(-a*c*x + c)*exp(3*atanh(a*x))
    F = -2*(a*x + 1)**(sympy.S(5)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(5*a**2*c*(-a*x + 1)**(sympy.S(3)/2)) - 2*(a*x + 1)**(sympy.S(3)/2)*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**2*c*(-a*x + 1)**(sympy.S(3)/2)) - 4*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(a**2*c*(-a*x + 1)**(sympy.S(3)/2)) + 4*sqrt(2)*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(a**2*c*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_402():
    f = sqrt(-a*c*x + c)*exp(3*atanh(a*x))
    F = 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c)*sqrt(-a**2*x**2 + 1)/(2*sqrt(-a*c*x + c)))/a - 2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*(-a*c*x + c)**(sympy.S(3)/2)) - 4*c*sqrt(-a**2*x**2 + 1)/(a*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_403():
    f = sqrt(-a*c*x + c)*exp(3*atanh(a*x))/x
    F = -2*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(c*(-a*x + 1)**(sympy.S(3)/2)) + 4*sqrt(2)*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(c*(-a*x + 1)**(sympy.S(3)/2)) - 2*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(c*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_404():
    f = sqrt(-a*c*x + c)*exp(3*atanh(a*x))/x**2
    F = 4*sqrt(2)*a*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(c*(-a*x + 1)**(sympy.S(3)/2)) - 5*a*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(c*(-a*x + 1)**(sympy.S(3)/2)) - sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(c*x*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_405():
    f = sqrt(-a*c*x + c)*exp(3*atanh(a*x))/x**3
    F = 4*sqrt(2)*a**2*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(c*(-a*x + 1)**(sympy.S(3)/2)) - 23*a**2*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(4*c*(-a*x + 1)**(sympy.S(3)/2)) - 9*a*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(4*c*x*(-a*x + 1)**(sympy.S(3)/2)) - sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(2*c*x**2*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_406():
    f = sqrt(-a*c*x + c)*exp(3*atanh(a*x))/x**4
    F = 4*sqrt(2)*a**3*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(c*(-a*x + 1)**(sympy.S(3)/2)) - 45*a**3*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(8*c*(-a*x + 1)**(sympy.S(3)/2)) - 19*a**2*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(8*c*x*(-a*x + 1)**(sympy.S(3)/2)) - 13*a*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(12*c*x**2*(-a*x + 1)**(sympy.S(3)/2)) - sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(3*c*x**3*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_407():
    f = sqrt(-a*c*x + c)*exp(3*atanh(a*x))/x**5
    F = 4*sqrt(2)*a**4*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*x + 1)/2)/(c*(-a*x + 1)**(sympy.S(3)/2)) - 363*a**4*(-a*c*x + c)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(64*c*(-a*x + 1)**(sympy.S(3)/2)) - 149*a**3*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(64*c*x*(-a*x + 1)**(sympy.S(3)/2)) - 107*a**2*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(96*c*x**2*(-a*x + 1)**(sympy.S(3)/2)) - 17*a*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(24*c*x**3*(-a*x + 1)**(sympy.S(3)/2)) - sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(4*c*x**4*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_408():
    f = x**m*sqrt(-a*c*x + c)*exp(-atanh(a*x))
    F = -2*c*x**(m + 1)*sqrt(-a**2*x**2 + 1)/((2*m + 3)*sqrt(-a*c*x + c)) + x**m*(8*m + 10)*(a*x + 1)*sqrt(-a*c*x + c)*hyper((sympy.S.Half, -m), (sympy.S(3)/2,), a*x + 1)/(a*(-a*x)**m*(2*m + 3)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_409():
    f = x**2*sqrt(-a*c*x + c)*exp(-atanh(a*x))
    F = -2*c*x**3*sqrt(-a**2*x**2 + 1)/(7*sqrt(-a*c*x + c)) + 26*c*x**2*sqrt(-a**2*x**2 + 1)/(35*a*sqrt(-a*c*x + c)) + 104*c*sqrt(-a**2*x**2 + 1)/(105*a**3*sqrt(-a*c*x + c)) + 104*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(105*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_410():
    f = x*sqrt(-a*c*x + c)*exp(-atanh(a*x))
    F = -8*c*sqrt(-a**2*x**2 + 1)/(5*a**2*sqrt(-a*c*x + c)) - 2*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(5*a**2) - 2*sqrt(-a**2*x**2 + 1)*(-a*c*x + c)**(sympy.S(3)/2)/(5*a**2*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_411():
    f = sqrt(-a*c*x + c)*exp(-atanh(a*x))
    F = 8*c*sqrt(-a**2*x**2 + 1)/(3*a*sqrt(-a*c*x + c)) + 2*sqrt(-a**2*x**2 + 1)*sqrt(-a*c*x + c)/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_412():
    f = sqrt(-a*c*x + c)*exp(-atanh(a*x))/x
    F = -2*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c)) - 2*c*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_413():
    f = sqrt(-a*c*x + c)*exp(-atanh(a*x))/x**2
    F = 3*a*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c)) - c*sqrt(-a**2*x**2 + 1)/(x*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_414():
    f = sqrt(-a*c*x + c)*exp(-atanh(a*x))/x**3
    F = -7*a**2*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c))/4 + 7*a*c*sqrt(-a**2*x**2 + 1)/(4*x*sqrt(-a*c*x + c)) - c*sqrt(-a**2*x**2 + 1)/(2*x**2*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_415():
    f = sqrt(-a*c*x + c)*exp(-atanh(a*x))/x**4
    F = 11*a**3*sqrt(c)*atanh(sqrt(c)*sqrt(-a**2*x**2 + 1)/sqrt(-a*c*x + c))/8 - 11*a**2*c*sqrt(-a**2*x**2 + 1)/(8*x*sqrt(-a*c*x + c)) + 11*a*c*sqrt(-a**2*x**2 + 1)/(12*x**2*sqrt(-a*c*x + c)) - c*sqrt(-a**2*x**2 + 1)/(3*x**3*sqrt(-a*c*x + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_416():
    f = x**3*sqrt(-a*c*x + c)*exp(-2*atanh(a*x))
    F = 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a**4 - 4*sqrt(-a*c*x + c)/a**4 - 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**4*c) - 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**4*c**2) + 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**4*c**3) - 2*(-a*c*x + c)**(sympy.S(9)/2)/(9*a**4*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_417():
    f = x**2*sqrt(-a*c*x + c)*exp(-2*atanh(a*x))
    F = -4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a**3 + 4*sqrt(-a*c*x + c)/a**3 + 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**3*c) + 2*(-a*c*x + c)**(sympy.S(7)/2)/(7*a**3*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_418():
    f = x*sqrt(-a*c*x + c)*exp(-2*atanh(a*x))
    F = 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a**2 - 4*sqrt(-a*c*x + c)/a**2 - 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a**2*c) - 2*(-a*c*x + c)**(sympy.S(5)/2)/(5*a**2*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_419():
    f = sqrt(-a*c*x + c)*exp(-2*atanh(a*x))
    F = -4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c)))/a + 4*sqrt(-a*c*x + c)/a + 2*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_420():
    f = sqrt(-a*c*x + c)*exp(-2*atanh(a*x))/x
    F = -2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) + 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) - 2*sqrt(-a*c*x + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_421():
    f = sqrt(-a*c*x + c)*exp(-2*atanh(a*x))/x**2
    F = 5*a*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c)) - 4*sqrt(2)*a*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) - sqrt(-a*c*x + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_422():
    f = sqrt(-a*c*x + c)*exp(-2*atanh(a*x))/x**3
    F = -23*a**2*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/4 + 4*sqrt(2)*a**2*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) + 9*a*sqrt(-a*c*x + c)/(4*x) - sqrt(-a*c*x + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_423():
    f = sqrt(-a*c*x + c)*exp(-2*atanh(a*x))/x**4
    F = 45*a**3*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/8 - 4*sqrt(2)*a**3*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) - 19*a**2*sqrt(-a*c*x + c)/(8*x) + 13*a*sqrt(-a*c*x + c)/(12*x**2) - sqrt(-a*c*x + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_424():
    f = sqrt(-a*c*x + c)*exp(-2*atanh(a*x))/x**5
    F = -363*a**4*sqrt(c)*atanh(sqrt(-a*c*x + c)/sqrt(c))/64 + 4*sqrt(2)*a**4*sqrt(c)*atanh(sqrt(2)*sqrt(-a*c*x + c)/(2*sqrt(c))) + 149*a**3*sqrt(-a*c*x + c)/(64*x) - 107*a**2*sqrt(-a*c*x + c)/(96*x**2) + 17*a*sqrt(-a*c*x + c)/(24*x**3) - sqrt(-a*c*x + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_425():
    f = x**3*sqrt(-a*c*x + c)*exp(-3*atanh(a*x))
    F = 2*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(9)/2)/(9*a**4*(-a*c*x + c)**(sympy.S(3)/2)) - 2*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(7)/2)/(a**4*(-a*c*x + c)**(sympy.S(3)/2)) + 38*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(5)/2)/(5*a**4*(-a*c*x + c)**(sympy.S(3)/2)) - 50*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)/(3*a**4*(-a*c*x + c)**(sympy.S(3)/2)) + 32*c**2*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(a**4*(-a*c*x + c)**(sympy.S(3)/2)) + 8*c**2*(-a*x + 1)**(sympy.S(3)/2)/(a**4*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_426():
    f = x**2*sqrt(-a*c*x + c)*exp(-3*atanh(a*x))
    F = 2*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(7)/2)/(7*a**3*(-a*c*x + c)**(sympy.S(3)/2)) - 12*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(5)/2)/(5*a**3*(-a*c*x + c)**(sympy.S(3)/2)) + 26*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)/(3*a**3*(-a*c*x + c)**(sympy.S(3)/2)) - 24*c**2*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(a**3*(-a*c*x + c)**(sympy.S(3)/2)) - 8*c**2*(-a*x + 1)**(sympy.S(3)/2)/(a**3*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_427():
    f = x*sqrt(-a*c*x + c)*exp(-3*atanh(a*x))
    F = 2*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(5)/2)/(5*a**2*(-a*c*x + c)**(sympy.S(3)/2)) - 10*c**2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)/(3*a**2*(-a*c*x + c)**(sympy.S(3)/2)) + 16*c**2*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(a**2*(-a*c*x + c)**(sympy.S(3)/2)) + 8*c**2*(-a*x + 1)**(sympy.S(3)/2)/(a**2*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_428():
    f = sqrt(-a*c*x + c)*exp(-3*atanh(a*x))
    F = -64*sqrt(-a*c*x + c)/(3*a*sqrt(-a**2*x**2 + 1)) + 16*(-a*c*x + c)**(sympy.S(3)/2)/(3*a*c*sqrt(-a**2*x**2 + 1)) + 2*(-a*c*x + c)**(sympy.S(5)/2)/(3*a*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_429():
    f = sqrt(-a*c*x + c)*exp(-3*atanh(a*x))/x
    F = 2*c**2*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(-a*c*x + c)**(sympy.S(3)/2) - 2*c**2*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(-a*c*x + c)**(sympy.S(3)/2) + 8*c**2*(-a*x + 1)**(sympy.S(3)/2)/(sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_430():
    f = sqrt(-a*c*x + c)*exp(-3*atanh(a*x))/x**2
    F = 7*a*c**2*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(-a*c*x + c)**(sympy.S(3)/2) - 9*a*c**2*(-a*x + 1)**(sympy.S(3)/2)/(sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) - c**2*(-a*x + 1)**(sympy.S(3)/2)/(x*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_431():
    f = sqrt(-a*c*x + c)*exp(-3*atanh(a*x))/x**3
    F = -47*a**2*c**2*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(4*(-a*c*x + c)**(sympy.S(3)/2)) + 47*a**2*c**2*(-a*x + 1)**(sympy.S(3)/2)/(4*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) + 13*a*c**2*(-a*x + 1)**(sympy.S(3)/2)/(4*x*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) - c**2*(-a*x + 1)**(sympy.S(3)/2)/(2*x**2*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_432():
    f = sqrt(-a*c*x + c)*exp(-3*atanh(a*x))/x**4
    F = 119*a**3*c**2*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(8*(-a*c*x + c)**(sympy.S(3)/2)) - 119*a**3*c**2*(-a*x + 1)**(sympy.S(3)/2)/(8*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) - 119*a**2*c**2*(-a*x + 1)**(sympy.S(3)/2)/(24*x*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) + 19*a*c**2*(-a*x + 1)**(sympy.S(3)/2)/(12*x**2*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) - c**2*(-a*x + 1)**(sympy.S(3)/2)/(3*x**3*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_433():
    f = sqrt(-a*c*x + c)*exp(-3*atanh(a*x))/x**5
    F = -1115*a**4*c**2*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(a*x + 1))/(64*(-a*c*x + c)**(sympy.S(3)/2)) + 1115*a**4*c**2*(-a*x + 1)**(sympy.S(3)/2)/(64*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) + 1115*a**3*c**2*(-a*x + 1)**(sympy.S(3)/2)/(192*x*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) - 223*a**2*c**2*(-a*x + 1)**(sympy.S(3)/2)/(96*x**2*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) + 25*a*c**2*(-a*x + 1)**(sympy.S(3)/2)/(24*x**3*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2)) - c**2*(-a*x + 1)**(sympy.S(3)/2)/(4*x**4*sqrt(a*x + 1)*(-a*c*x + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_434():
    f = (-a*c*x + c)**p*exp(-2*p*atanh(a*x))
    F = -(-a*x + 1)**p*(-a*c*x + c)**(p + 1)*hyper((p, 2*p + 1), (2*p + 2,), -a*x/2 + sympy.S.Half)/(2**p*a*c*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_435():
    f = (-a*c*x + c)**p*exp(2*p*atanh(a*x))
    F = (a*x + 1)**(p + 1)*(-a*c*x + c)**p/(a*(p + 1)*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_436():
    f = (-a*c*x + c)**p*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*(-a*c*x + c)**(p + 1)*hyper((-n/2, -n/2 + p + 1), (-n/2 + p + 2,), -a*x/2 + sympy.S.Half)/(a*c*(-a*x + 1)**(n/2)*(-n + 2*p + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_437():
    f = (-a*c*x + c)**3*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*c**3*(-a*x + 1)**(4 - n/2)*hyper((-n/2, 4 - n/2), (5 - n/2,), -a*x/2 + sympy.S.Half)/(a*(8 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_438():
    f = (-a*c*x + c)**2*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*c**2*(-a*x + 1)**(3 - n/2)*hyper((-n/2, 3 - n/2), (4 - n/2,), -a*x/2 + sympy.S.Half)/(a*(6 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_439():
    f = (-a*c*x + c)*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*c*(-a*x + 1)**(2 - n/2)*hyper((-n/2, 2 - n/2), (3 - n/2,), -a*x/2 + sympy.S.Half)/(a*(4 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_440():
    f = exp(n*atanh(a*x))/(-a*c*x + c)
    F = 2**(n/2 + 1)*hyper((-n/2, -n/2), (1 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_441():
    f = exp(n*atanh(a*x))/(-a*c*x + c)**2
    F = (-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 + 1)/(a*c**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_442():
    f = exp(n*atanh(a*x))/(-a*c*x + c)**3
    F = (-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 + 1)/(a*c**3*(n**2 + 6*n + 8)) + (-a*x + 1)**(-n/2 - 2)*(a*x + 1)**(n/2 + 1)/(a*c**3*(n + 4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_443():
    f = exp(n*atanh(a*x))/(-a*c*x + c)**4
    F = (-a*x + 1)**(-n/2 - 3)*(a*x + 1)**(n/2 + 1)/(a*c**4*(n + 6)) + 2*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 + 1)/(a*c**4*(n + 6)*(n**2 + 6*n + 8)) + 2*(-a*x + 1)**(-n/2 - 2)*(a*x + 1)**(n/2 + 1)/(a*c**4*(n + 4)*(n + 6))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_444():
    f = (c - c/(a*x))**p*exp(atanh(a*x))
    F = x*(c - c/(a*x))**p*appellf1(1 - p, sympy.S(-1)/2, sympy.S.Half - p, 2 - p, -a*x, a*x)/((1 - p)*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_445():
    f = (c - c/(a*x))**4*exp(atanh(a*x))
    F = -3*c**4*asin(a*x)/a - c**4*atanh(sqrt(-a**2*x**2 + 1))/(2*a) - c**4*(-a*x + 6)*sqrt(-a**2*x**2 + 1)/(2*a**2*x) + 3*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**3*x**2) - c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_446():
    f = (c - c/(a*x))**3*exp(atanh(a*x))
    F = -2*c**3*asin(a*x)/a + c**3*atanh(sqrt(-a**2*x**2 + 1))/(2*a) - c**3*(a*x + 4)*sqrt(-a**2*x**2 + 1)/(2*a**2*x) + c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_447():
    f = (c - c/(a*x))**2*exp(atanh(a*x))
    F = -c**2*asin(a*x)/a + c**2*atanh(sqrt(-a**2*x**2 + 1))/a - c**2*(a*x + 1)*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_448():
    f = (c - c/(a*x))*exp(atanh(a*x))
    F = -c*sqrt(-a**2*x**2 + 1)/a + c*atanh(sqrt(-a**2*x**2 + 1))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_449():
    f = exp(atanh(a*x))/(c - c/(a*x))
    F = -2*sqrt(-a**2*x**2 + 1)/(a*c) + 2*asin(a*x)/(a*c) - (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_450():
    f = exp(atanh(a*x))/(c - c/(a*x))**2
    F = 3*asin(a*x)/(a*c**2) - 6*sqrt(-a**2*x**2 + 1)/(a*c**2*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a*c**2*(-a*x + 1)**2) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_451():
    f = exp(atanh(a*x))/(c - c/(a*x))**3
    F = 4*asin(a*x)/(a*c**3) - 8*sqrt(-a**2*x**2 + 1)/(a*c**3*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a*c**3*(-a*x + 1)**2) + 14*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(15*a*c**3*(-a*x + 1)**3) - (-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_452():
    f = exp(atanh(a*x))/(c - c/(a*x))**4
    F = 5*asin(a*x)/(a*c**4) - 10*sqrt(-a**2*x**2 + 1)/(a*c**4*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a*c**4*(-a*x + 1)**2) + 184*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(105*a*c**4*(-a*x + 1)**3) - 26*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(35*a*c**4*(-a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(7*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_453():
    f = (c - c/(a*x))**p*exp(2*atanh(a*x))
    F = -x*(c - c/(a*x))**p - (2 - p)*(c - c/(a*x))**p*hyper((1, p), (p + 1,), 1 - 1/(a*x))/(a*p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_454():
    f = (c - c/(a*x))**5*exp(2*atanh(a*x))
    F = -c**5*x + 3*c**5*log(x)/a + 2*c**5/(a**2*x) + c**5/(a**3*x**2) - c**5/(a**4*x**3) + c**5/(4*a**5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_455():
    f = (c - c/(a*x))**4*exp(2*atanh(a*x))
    F = -c**4*x + 2*c**4*log(x)/a + c**4/(a**3*x**2) - c**4/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_456():
    f = (c - c/(a*x))**3*exp(2*atanh(a*x))
    F = -c**3*x + c**3*log(x)/a - c**3/(a**2*x) + c**3/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_457():
    f = (c - c/(a*x))**2*exp(2*atanh(a*x))
    F = -c**2*x - c**2/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_458():
    f = (c - c/(a*x))*exp(2*atanh(a*x))
    F = -c*x - c*log(x)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_459():
    f = exp(2*atanh(a*x))/(c - c/(a*x))
    F = -x/c - 3*log(-a*x + 1)/(a*c) - 2/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_460():
    f = exp(2*atanh(a*x))/(c - c/(a*x))**2
    F = -x/c**2 - 4*log(-a*x + 1)/(a*c**2) - 5/(a*c**2*(-a*x + 1)) + 1/(a*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_461():
    f = exp(2*atanh(a*x))/(c - c/(a*x))**3
    F = -x/c**3 - 5*log(-a*x + 1)/(a*c**3) - 9/(a*c**3*(-a*x + 1)) + 7/(2*a*c**3*(-a*x + 1)**2) - 2/(3*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_462():
    f = exp(2*atanh(a*x))/(c - c/(a*x))**4
    F = -x/c**4 - 6*log(-a*x + 1)/(a*c**4) - 14/(a*c**4*(-a*x + 1)) + 8/(a*c**4*(-a*x + 1)**2) - 3/(a*c**4*(-a*x + 1)**3) + 1/(2*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_463():
    f = (c - c/(a*x))**4*exp(3*atanh(a*x))
    F = c**4*asin(a*x)/a - 3*c**4*atanh(sqrt(-a**2*x**2 + 1))/(2*a) + c**4*(3*a*x + 2)*sqrt(-a**2*x**2 + 1)/(2*a**2*x) - c**4*(-3*a*x + 2)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_464():
    f = (c - c/(a*x))**3*exp(3*atanh(a*x))
    F = 3*c**3*sqrt(-a**2*x**2 + 1)/(2*a) - 3*c**3*atanh(sqrt(-a**2*x**2 + 1))/(2*a) + c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_465():
    f = (c - c/(a*x))**2*exp(3*atanh(a*x))
    F = -c**2*asin(a*x)/a - c**2*atanh(sqrt(-a**2*x**2 + 1))/a - c**2*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_466():
    f = (c - c/(a*x))*exp(3*atanh(a*x))
    F = c*sqrt(-a**2*x**2 + 1)/a - 2*c*asin(a*x)/a + c*atanh(sqrt(-a**2*x**2 + 1))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_467():
    f = exp(3*atanh(a*x))/(c - c/(a*x))
    F = 4*sqrt(-a**2*x**2 + 1)/(a*c) - 4*asin(a*x)/(a*c) + 8*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a*c*(-a*x + 1)**2) - (-a**2*x**2 + 1)**(sympy.S(5)/2)/(3*a*c*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_468():
    f = exp(3*atanh(a*x))/(c - c/(a*x))**2
    F = (a*x + 1)**5/(5*a*c**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 2*(a*x + 1)**4/(3*a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 10*(a*x + 1)**2/(3*a*c**2*sqrt(-a**2*x**2 + 1)) + 5*sqrt(-a**2*x**2 + 1)/(a*c**2) - 5*asin(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_469():
    f = exp(3*atanh(a*x))/(c - c/(a*x))**3
    F = -(a*x + 1)**6/(7*a*c**3*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + 4*(a*x + 1)**5/(7*a*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - (a*x + 1)**4/(a*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 4*(a*x + 1)**2/(a*c**3*sqrt(-a**2*x**2 + 1)) + 6*sqrt(-a**2*x**2 + 1)/(a*c**3) - 6*asin(a*x)/(a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_470():
    f = exp(3*atanh(a*x))/(c - c/(a*x))**4
    F = (a*x + 1)**7/(9*a*c**4*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - 34*(a*x + 1)**6/(63*a*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + 344*(a*x + 1)**5/(315*a*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 4*(a*x + 1)**4/(3*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 14*(a*x + 1)**2/(3*a*c**4*sqrt(-a**2*x**2 + 1)) + 7*sqrt(-a**2*x**2 + 1)/(a*c**4) - 7*asin(a*x)/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_471():
    f = (c - c/(a*x))**p*exp(4*atanh(a*x))
    F = c*x*(c - c/(a*x))**(p - 1) - c*(5 - p)*(c - c/(a*x))**(p - 1)/(a*(1 - p)) + (4 - p)*(c - c/(a*x))**p*hyper((1, p), (p + 1,), 1 - 1/(a*x))/(a*p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_472():
    f = (c - c/(a*x))**5*exp(4*atanh(a*x))
    F = c**5*x - c**5*log(x)/a + 2*c**5/(a**2*x) - c**5/(a**3*x**2) - c**5/(3*a**4*x**3) + c**5/(4*a**5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_473():
    f = (c - c/(a*x))**4*exp(4*atanh(a*x))
    F = c**4*x + 2*c**4/(a**2*x) - c**4/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_474():
    f = (c - c/(a*x))**3*exp(4*atanh(a*x))
    F = c**3*x + c**3*log(x)/a + c**3/(a**2*x) + c**3/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_475():
    f = (c - c/(a*x))**2*exp(4*atanh(a*x))
    F = c**2*x + 2*c**2*log(x)/a - c**2/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_476():
    f = (c - c/(a*x))*exp(4*atanh(a*x))
    F = c*x - c*log(x)/a + 4*c*log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_477():
    f = exp(4*atanh(a*x))/(c - c/(a*x))
    F = x/c + 5*log(-a*x + 1)/(a*c) + 8/(a*c*(-a*x + 1)) - 2/(a*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_478():
    f = exp(4*atanh(a*x))/(c - c/(a*x))**2
    F = x/c**2 + 6*log(-a*x + 1)/(a*c**2) + 13/(a*c**2*(-a*x + 1)) - 6/(a*c**2*(-a*x + 1)**2) + 4/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_479():
    f = exp(4*atanh(a*x))/(c - c/(a*x))**3
    F = x/c**3 + 7*log(-a*x + 1)/(a*c**3) + 19/(a*c**3*(-a*x + 1)) - 25/(2*a*c**3*(-a*x + 1)**2) + 16/(3*a*c**3*(-a*x + 1)**3) - 1/(a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_480():
    f = exp(4*atanh(a*x))/(c - c/(a*x))**4
    F = x/c**4 + 8*log(-a*x + 1)/(a*c**4) + 26/(a*c**4*(-a*x + 1)) - 22/(a*c**4*(-a*x + 1)**2) + 41/(3*a*c**4*(-a*x + 1)**3) - 5/(a*c**4*(-a*x + 1)**4) + 4/(5*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_481():
    f = (c - c/(a*x))**p*exp(-atanh(a*x))
    F = x*(c - c/(a*x))**p*appellf1(1 - p, sympy.S.Half, -p + sympy.S(-1)/2, 2 - p, -a*x, a*x)/((1 - p)*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_482():
    f = (c - c/(a*x))**4*exp(-atanh(a*x))
    F = c**4*sqrt(-a**2*x**2 + 1)/a + 5*c**4*asin(a*x)/a + 25*c**4*atanh(sqrt(-a**2*x**2 + 1))/(2*a) - 32*c**4*sqrt(-a**2*x**2 + 1)/(3*a**2*x) + 5*c**4*sqrt(-a**2*x**2 + 1)/(2*a**3*x**2) - c**4*sqrt(-a**2*x**2 + 1)/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_483():
    f = (c - c/(a*x))**3*exp(-atanh(a*x))
    F = c**3*sqrt(-a**2*x**2 + 1)/a + 4*c**3*asin(a*x)/a + 13*c**3*atanh(sqrt(-a**2*x**2 + 1))/(2*a) - 4*c**3*sqrt(-a**2*x**2 + 1)/(a**2*x) + c**3*sqrt(-a**2*x**2 + 1)/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_484():
    f = (c - c/(a*x))**2*exp(-atanh(a*x))
    F = c**2*sqrt(-a**2*x**2 + 1)/a + 3*c**2*asin(a*x)/a + 3*c**2*atanh(sqrt(-a**2*x**2 + 1))/a - c**2*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_485():
    f = (c - c/(a*x))*exp(-atanh(a*x))
    F = c*sqrt(-a**2*x**2 + 1)/a + 2*c*asin(a*x)/a + c*atanh(sqrt(-a**2*x**2 + 1))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_486():
    f = exp(-atanh(a*x))/(c - c/(a*x))
    F = sqrt(-a**2*x**2 + 1)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_487():
    f = exp(-atanh(a*x))/(c - c/(a*x))**2
    F = sqrt(-a**2*x**2 + 1)/(a*c**2) - asin(a*x)/(a*c**2) + sqrt(-a**2*x**2 + 1)/(a*c**2*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_488():
    f = exp(-atanh(a*x))/(c - c/(a*x))**3
    F = -(a*x + 1)**2/(3*a*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (8*a*x + 8)/(3*a*c**3*sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(a*c**3) - 2*asin(a*x)/(a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_489():
    f = exp(-atanh(a*x))/(c - c/(a*x))**4
    F = (a*x + 1)**3/(5*a*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 6*(a*x + 1)**2/(5*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (24*a*x + 24)/(5*a*c**4*sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(a*c**4) - 3*asin(a*x)/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_490():
    f = (c - c/(a*x))**p*exp(-2*atanh(a*x))
    F = -x*(c - c/(a*x))**(p + 2)/c**2 + (c - c/(a*x))**(p + 2)*hyper((1, p + 2), (p + 3,), 1 - 1/(a*x))/(a*c**2) - (c - c/(a*x))**(p + 2)*hyper((1, p + 2), (p + 3,), (a - 1/x)/(2*a))/(2*a*c**2*(p + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_491():
    f = (c - c/(a*x))**4*exp(-2*atanh(a*x))
    F = -c**4*x - 26*c**4*log(x)/a + 32*c**4*log(a*x + 1)/a - 16*c**4/(a**2*x) + 3*c**4/(a**3*x**2) - c**4/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_492():
    f = (c - c/(a*x))**3*exp(-2*atanh(a*x))
    F = -c**3*x - 11*c**3*log(x)/a + 16*c**3*log(a*x + 1)/a - 5*c**3/(a**2*x) + c**3/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_493():
    f = (c - c/(a*x))**2*exp(-2*atanh(a*x))
    F = -c**2*x - 4*c**2*log(x)/a + 8*c**2*log(a*x + 1)/a - c**2/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_494():
    f = (c - c/(a*x))*exp(-2*atanh(a*x))
    F = -c*x - c*log(x)/a + 4*c*log(a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_495():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))
    F = -x/c + log(a*x + 1)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_496():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))**2
    F = -x/c**2 + atanh(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_497():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))**3
    F = -x/c**3 - 5*log(-a*x + 1)/(4*a*c**3) + log(a*x + 1)/(4*a*c**3) - 1/(2*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_498():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))**4
    F = -x/c**4 - 17*log(-a*x + 1)/(8*a*c**4) + log(a*x + 1)/(8*a*c**4) - 7/(4*a*c**4*(-a*x + 1)) + 1/(4*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_499():
    f = (c - c/(a*x))**3*exp(-3*atanh(a*x))
    F = -32*c**3*(-a*x + 1)/(a*sqrt(-a**2*x**2 + 1)) - c**3*sqrt(-a**2*x**2 + 1)/a - 6*c**3*asin(a*x)/a + 33*c**3*atanh(sqrt(-a**2*x**2 + 1))/(2*a) - 6*c**3*sqrt(-a**2*x**2 + 1)/(a**2*x) + c**3*sqrt(-a**2*x**2 + 1)/(2*a**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_500():
    f = (c - c/(a*x))**2*exp(-3*atanh(a*x))
    F = -16*c**2*(-a*x + 1)/(a*sqrt(-a**2*x**2 + 1)) - c**2*sqrt(-a**2*x**2 + 1)/a - 5*c**2*asin(a*x)/a + 5*c**2*atanh(sqrt(-a**2*x**2 + 1))/a - c**2*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_501():
    f = (c - c/(a*x))*exp(-3*atanh(a*x))
    F = -8*c*(-a*x + 1)/(a*sqrt(-a**2*x**2 + 1)) - c*sqrt(-a**2*x**2 + 1)/a - 4*c*asin(a*x)/a + c*atanh(sqrt(-a**2*x**2 + 1))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_502():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))
    F = -(-a*x + 1)**2/(a*c*sqrt(-a**2*x**2 + 1)) - 2*sqrt(-a**2*x**2 + 1)/(a*c) - 2*asin(a*x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_503():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))**2
    F = -(-a*x + 1)/(a*c**2*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*x**2 + 1)/(a*c**2) - asin(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_504():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))**3
    F = -sqrt(-a**2*x**2 + 1)/(a*c**3) - 1/(a*c**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_505():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))**4
    F = a**2*x**3*(a*x + 1)/(3*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(4*a*x + 3)/(3*c**4*sqrt(-a**2*x**2 + 1)) - 8*sqrt(-a**2*x**2 + 1)/(3*a*c**4) + asin(a*x)/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_506():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))**5
    F = -(a*x + 1)**2/(5*a*c**5*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + (22*a*x + 22)/(15*a*c**5*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - (46*a*x + 60)/(15*a*c**5*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*x**2 + 1)/(a*c**5) + 2*asin(a*x)/(a*c**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_507():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(atanh(a*x))
    F = -7*a**(sympy.S(7)/2)*x**(sympy.S(9)/2)*(c - c/(a*x))**(sympy.S(9)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(9)/2) - a**3*x**4*(c - c/(a*x))**(sympy.S(9)/2)*(-227*a*x + 54)*sqrt(a*x + 1)/(105*(-a*x + 1)**(sympy.S(9)/2)) - 10*a**2*x**3*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)/(21*(-a*x + 1)**(sympy.S(5)/2)) + 2*a*x**2*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)/(5*(-a*x + 1)**(sympy.S(3)/2)) - 2*x*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)/(7*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_508():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(atanh(a*x))
    F = 5*a**(sympy.S(5)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(7)/2) - a**2*x**3*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)*(31*a*x + 18)/(15*(-a*x + 1)**(sympy.S(7)/2)) + 2*a*x**2*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)/(3*(-a*x + 1)**(sympy.S(3)/2)) - 2*x*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)/(5*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_509():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(atanh(a*x))
    F = -3*a**(sympy.S(3)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(5)/2) - 3*a**2*x**3*(c - c/(a*x))**(sympy.S(5)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(5)/2) + 4*a*x**2*(c - c/(a*x))**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(3)/2)/(-a*x + 1)**(sympy.S(5)/2) - 2*x*(c - c/(a*x))**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(3)/2)/(3*(-a*x + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_510():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(atanh(a*x))
    F = sqrt(a)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(3)/2) + a*x**2*(c - c/(a*x))**(sympy.S(3)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(3)/2) - 2*x*(c - c/(a*x))**(sympy.S(3)/2)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(-a*x + 1)**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_511():
    f = exp(atanh(a*x))/sqrt(c - c/(a*x))
    F = -sqrt(-a*x + 1)*sqrt(a*x + 1)/(a*sqrt(c - c/(a*x))) - 3*sqrt(-a*x + 1)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(3)/2)*sqrt(x)*sqrt(c - c/(a*x))) + 2*sqrt(2)*sqrt(-a*x + 1)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(a**(sympy.S(3)/2)*sqrt(x)*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_512():
    f = exp(atanh(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = sqrt(-a*x + 1)*sqrt(a*x + 1)/(a*(c - c/(a*x))**(sympy.S(3)/2)) + 2*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(a**2*x*(c - c/(a*x))**(sympy.S(3)/2)) + 5*(-a*x + 1)**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2)) - 7*sqrt(2)*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(2*a**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_513():
    f = exp(atanh(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = sqrt(-a*x + 1)*sqrt(a*x + 1)/(2*a*(c - c/(a*x))**(sympy.S(5)/2)) - 11*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(8*a**2*x*(c - c/(a*x))**(sympy.S(5)/2)) - 23*(-a*x + 1)**(sympy.S(5)/2)*sqrt(a*x + 1)/(8*a**3*x**2*(c - c/(a*x))**(sympy.S(5)/2)) - 7*(-a*x + 1)**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)) + 79*sqrt(2)*(-a*x + 1)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(16*a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_514():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(9)/2) + 5*c**(sympy.S(9)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 5*c**4*sqrt(c - c/(a*x))/a - 5*c**3*(c - c/(a*x))**(sympy.S(3)/2)/(3*a) - c**2*(c - c/(a*x))**(sympy.S(5)/2)/a - 5*c*(c - c/(a*x))**(sympy.S(7)/2)/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_515():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(7)/2) + 3*c**(sympy.S(7)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 3*c**3*sqrt(c - c/(a*x))/a - c**2*(c - c/(a*x))**(sympy.S(3)/2)/a - 3*c*(c - c/(a*x))**(sympy.S(5)/2)/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_516():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(5)/2) + c**(sympy.S(5)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - c**2*sqrt(c - c/(a*x))/a - c*(c - c/(a*x))**(sympy.S(3)/2)/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_517():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(3)/2) - c**(sympy.S(3)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a + c*sqrt(c - c/(a*x))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_518():
    f = sqrt(c - c/(a*x))*exp(2*atanh(a*x))
    F = -x*sqrt(c - c/(a*x)) - 3*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_519():
    f = exp(2*atanh(a*x))/sqrt(c - c/(a*x))
    F = -x/sqrt(c - c/(a*x)) + 5/(a*sqrt(c - c/(a*x))) - 5*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_520():
    f = exp(2*atanh(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = -x/(c - c/(a*x))**(sympy.S(3)/2) + 7/(3*a*(c - c/(a*x))**(sympy.S(3)/2)) + 7/(a*c*sqrt(c - c/(a*x))) - 7*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_521():
    f = exp(2*atanh(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = -x/(c - c/(a*x))**(sympy.S(5)/2) + 9/(5*a*(c - c/(a*x))**(sympy.S(5)/2)) + 3/(a*c*(c - c/(a*x))**(sympy.S(3)/2)) + 9/(a*c**2*sqrt(c - c/(a*x))) - 9*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_522():
    f = exp(2*atanh(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = -x/(c - c/(a*x))**(sympy.S(7)/2) + 11/(7*a*(c - c/(a*x))**(sympy.S(7)/2)) + 11/(5*a*c*(c - c/(a*x))**(sympy.S(5)/2)) + 11/(3*a*c**2*(c - c/(a*x))**(sympy.S(3)/2)) + 11/(a*c**3*sqrt(c - c/(a*x))) - 11*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_523():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(3*atanh(a*x))
    F = 3*a**(sympy.S(7)/2)*x**(sympy.S(9)/2)*(c - c/(a*x))**(sympy.S(9)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(9)/2) - 3*a**3*x**4*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(9)/2) + 3*a**2*x**3*(c - c/(a*x))**(sympy.S(9)/2)*(-17*a*x + 6)*(a*x + 1)**(sympy.S(3)/2)/(35*(-a*x + 1)**(sympy.S(9)/2)) + 6*a*x**2*(c - c/(a*x))**(sympy.S(9)/2)*(a*x + 1)**(sympy.S(3)/2)/(35*(-a*x + 1)**(sympy.S(5)/2)) - 2*x*(c - c/(a*x))**(sympy.S(9)/2)*(a*x + 1)**(sympy.S(3)/2)/(7*(-a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_524():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(3*atanh(a*x))
    F = -a**(sympy.S(5)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(7)/2) - a**3*x**4*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(7)/2) + 2*a**2*x**3*(c - c/(a*x))**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(3)/2)/(3*(-a*x + 1)**(sympy.S(7)/2)) + 4*a*x**2*(c - c/(a*x))**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(5)/2)/(3*(-a*x + 1)**(sympy.S(7)/2)) - 2*x*(c - c/(a*x))**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(5)/2)/(5*(-a*x + 1)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_525():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(3*atanh(a*x))
    F = -a**(sympy.S(3)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(5)/2) - a**2*x**3*(c - c/(a*x))**(sympy.S(5)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(5)/2) + 2*a*x**2*(c - c/(a*x))**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(3)/2)/(3*(-a*x + 1)**(sympy.S(5)/2)) - 2*x*(c - c/(a*x))**(sympy.S(5)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(3*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_526():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(3*atanh(a*x))
    F = 3*sqrt(a)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(3)/2) + 3*a*x**2*(c - c/(a*x))**(sympy.S(3)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(3)/2) - 2*x*(c - c/(a*x))**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)/(-a*x + 1)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_527():
    f = sqrt(c - c/(a*x))*exp(3*atanh(a*x))
    F = -x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1) - 5*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(sqrt(a)*sqrt(-a*x + 1)) + 4*sqrt(2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(sqrt(a)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_528():
    f = exp(3*atanh(a*x))/sqrt(c - c/(a*x))
    F = 2*sqrt(-a*x + 1)*sqrt(a*x + 1)/(a*sqrt(c - c/(a*x))) + (a*x + 1)**(sympy.S(3)/2)/(a*sqrt(c - c/(a*x))*sqrt(-a*x + 1)) + 7*sqrt(-a*x + 1)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(3)/2)*sqrt(x)*sqrt(c - c/(a*x))) - 5*sqrt(2)*sqrt(-a*x + 1)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(a**(sympy.S(3)/2)*sqrt(x)*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_529():
    f = exp(3*atanh(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = (a*x + 1)**(sympy.S(3)/2)/(2*a*(c - c/(a*x))**(sympy.S(3)/2)*sqrt(-a*x + 1)) - 21*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(8*a**2*x*(c - c/(a*x))**(sympy.S(3)/2)) - 9*sqrt(-a*x + 1)*(a*x + 1)**(sympy.S(3)/2)/(8*a**2*x*(c - c/(a*x))**(sympy.S(3)/2)) - 9*(-a*x + 1)**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2)) + 51*sqrt(2)*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(8*a**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_530():
    f = exp(3*atanh(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = (a*x + 1)**(sympy.S(3)/2)/(3*a*(c - c/(a*x))**(sympy.S(5)/2)*sqrt(-a*x + 1)) - 13*sqrt(-a*x + 1)*(a*x + 1)**(sympy.S(3)/2)/(24*a**2*x*(c - c/(a*x))**(sympy.S(5)/2)) + 103*(-a*x + 1)**(sympy.S(5)/2)*sqrt(a*x + 1)/(32*a**3*x**2*(c - c/(a*x))**(sympy.S(5)/2)) + 43*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)/(32*a**3*x**2*(c - c/(a*x))**(sympy.S(5)/2)) + 11*(-a*x + 1)**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)) - 249*sqrt(2)*(-a*x + 1)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(32*a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_531():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(-atanh(a*x))
    F = 11*a**(sympy.S(7)/2)*x**(sympy.S(9)/2)*(c - c/(a*x))**(sympy.S(9)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(9)/2) + a**3*x**4*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)*(521*a*x + 2718)/(105*(-a*x + 1)**(sympy.S(9)/2)) - 94*a**2*x**3*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)/(21*(-a*x + 1)**(sympy.S(5)/2)) + 6*a*x**2*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)/(5*(-a*x + 1)**(sympy.S(3)/2)) - 2*x*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(a*x + 1)/(7*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_532():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(-atanh(a*x))
    F = -9*a**(sympy.S(5)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(7)/2) - a**2*x**3*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)*(7*a*x + 66)/(5*(-a*x + 1)**(sympy.S(7)/2)) + 2*a*x**2*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(3)/2) - 2*x*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)/(5*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_533():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(-atanh(a*x))
    F = 7*a**(sympy.S(3)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(5)/2) + a*x**2*(c - c/(a*x))**(sympy.S(5)/2)*(-a*x + 18)*sqrt(a*x + 1)/(3*(-a*x + 1)**(sympy.S(5)/2)) - 2*x*(c - c/(a*x))**(sympy.S(5)/2)*sqrt(a*x + 1)/(3*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_534():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(-atanh(a*x))
    F = -5*sqrt(a)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(3)/2) + a*x**2*(c - c/(a*x))**(sympy.S(3)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(3)/2) - 2*x*(c - c/(a*x))**(sympy.S(3)/2)*sqrt(a*x + 1)/(-a*x + 1)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_535():
    f = sqrt(c - c/(a*x))*exp(-atanh(a*x))
    F = -x*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(-a*x + 1) + 3*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(sqrt(a)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_536():
    f = exp(-atanh(a*x))/sqrt(c - c/(a*x))
    F = sqrt(-a*x + 1)*sqrt(a*x + 1)/(a*sqrt(c - c/(a*x))) - sqrt(-a*x + 1)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(3)/2)*sqrt(x)*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_537():
    f = exp(-atanh(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = -(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(a**2*x*(c - c/(a*x))**(sympy.S(3)/2)) - (-a*x + 1)**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2)) + sqrt(2)*(-a*x + 1)**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(a**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_538():
    f = exp(-atanh(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = (-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(2*a**2*x*(c - c/(a*x))**(sympy.S(5)/2)) + 3*(-a*x + 1)**(sympy.S(5)/2)*sqrt(a*x + 1)/(2*a**3*x**2*(c - c/(a*x))**(sympy.S(5)/2)) + 3*(-a*x + 1)**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)) - 9*sqrt(2)*(-a*x + 1)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(4*a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_539():
    f = exp(-atanh(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = (-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(4*a**2*x*(c - c/(a*x))**(sympy.S(7)/2)) - 15*(-a*x + 1)**(sympy.S(5)/2)*sqrt(a*x + 1)/(16*a**3*x**2*(c - c/(a*x))**(sympy.S(7)/2)) - 35*(-a*x + 1)**(sympy.S(7)/2)*sqrt(a*x + 1)/(16*a**4*x**3*(c - c/(a*x))**(sympy.S(7)/2)) - 5*(-a*x + 1)**(sympy.S(7)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(9)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2)) + 115*sqrt(2)*(-a*x + 1)**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(32*a**(sympy.S(9)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_540():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(-2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(9)/2) + 13*c**(sympy.S(9)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 64*sqrt(2)*c**(sympy.S(9)/2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a + 51*c**4*sqrt(c - c/(a*x))/a + 19*c**3*(c - c/(a*x))**(sympy.S(3)/2)/(3*a) + 3*c**2*(c - c/(a*x))**(sympy.S(5)/2)/(5*a) - 5*c*(c - c/(a*x))**(sympy.S(7)/2)/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_541():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(-2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(7)/2) + 11*c**(sympy.S(7)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 32*sqrt(2)*c**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a + 21*c**3*sqrt(c - c/(a*x))/a + 5*c**2*(c - c/(a*x))**(sympy.S(3)/2)/(3*a) - 3*c*(c - c/(a*x))**(sympy.S(5)/2)/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_542():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(-2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(5)/2) + 9*c**(sympy.S(5)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 16*sqrt(2)*c**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a + 7*c**2*sqrt(c - c/(a*x))/a - c*(c - c/(a*x))**(sympy.S(3)/2)/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_543():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(-2*atanh(a*x))
    F = -x*(c - c/(a*x))**(sympy.S(3)/2) + 7*c**(sympy.S(3)/2)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 8*sqrt(2)*c**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a + c*sqrt(c - c/(a*x))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_544():
    f = sqrt(c - c/(a*x))*exp(-2*atanh(a*x))
    F = -x*sqrt(c - c/(a*x)) + 5*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_545():
    f = exp(-2*atanh(a*x))/sqrt(c - c/(a*x))
    F = -x*sqrt(c - c/(a*x))/c + 3*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*sqrt(c)) - 2*sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_546():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = -x*sqrt(c - c/(a*x))/c**2 + atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(a*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_547():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = -x/(c**2*sqrt(c - c/(a*x))) + 2/(a*c**2*sqrt(c - c/(a*x))) - atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(5)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(2*a*c**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_548():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = -x/(c**2*(c - c/(a*x))**(sympy.S(3)/2)) + 4/(3*a*c**2*(c - c/(a*x))**(sympy.S(3)/2)) + 7/(2*a*c**3*sqrt(c - c/(a*x))) - 3*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(7)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(4*a*c**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_549():
    f = exp(-2*atanh(a*x))/(c - c/(a*x))**(sympy.S(9)/2)
    F = -x/(c**2*(c - c/(a*x))**(sympy.S(5)/2)) + 6/(5*a*c**2*(c - c/(a*x))**(sympy.S(5)/2)) + 11/(6*a*c**3*(c - c/(a*x))**(sympy.S(3)/2)) + 21/(4*a*c**4*sqrt(c - c/(a*x))) - 5*atanh(sqrt(c - c/(a*x))/sqrt(c))/(a*c**(sympy.S(9)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/(8*a*c**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_550():
    f = (c - c/(a*x))**(sympy.S(9)/2)*exp(-3*atanh(a*x))
    F = -15*a**(sympy.S(7)/2)*x**(sympy.S(9)/2)*(c - c/(a*x))**(sympy.S(9)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(9)/2) + 5*a**4*x**5*(c - c/(a*x))**(sympy.S(9)/2)*(-109*a*x + 587)/(7*(-a*x + 1)**(sympy.S(9)/2)*sqrt(a*x + 1)) + 70*a**3*x**4*(c - c/(a*x))**(sympy.S(9)/2)/((-a*x + 1)**(sympy.S(5)/2)*sqrt(a*x + 1)) - 50*a**2*x**3*(c - c/(a*x))**(sympy.S(9)/2)/(7*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)) + 10*a*x**2*(c - c/(a*x))**(sympy.S(9)/2)/(7*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*x*(c - c/(a*x))**(sympy.S(9)/2)*sqrt(-a*x + 1)/(7*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_551():
    f = (c - c/(a*x))**(sympy.S(7)/2)*exp(-3*atanh(a*x))
    F = 13*a**(sympy.S(5)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(7)/2) - a**3*x**4*(c - c/(a*x))**(sympy.S(7)/2)*(-427*a*x + 2525)/(15*(-a*x + 1)**(sympy.S(7)/2)*sqrt(a*x + 1)) - 398*a**2*x**3*(c - c/(a*x))**(sympy.S(7)/2)/(15*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)) + 38*a*x**2*(c - c/(a*x))**(sympy.S(7)/2)/(15*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*x*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(-a*x + 1)/(5*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_552():
    f = (c - c/(a*x))**(sympy.S(5)/2)*exp(-3*atanh(a*x))
    F = -11*a**(sympy.S(3)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(5)/2) + a**2*x**3*(c - c/(a*x))**(sympy.S(5)/2)*(-25*a*x + 191)/(3*(-a*x + 1)**(sympy.S(5)/2)*sqrt(a*x + 1)) + 26*a*x**2*(c - c/(a*x))**(sympy.S(5)/2)/(3*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*x*(c - c/(a*x))**(sympy.S(5)/2)*sqrt(-a*x + 1)/(3*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_553():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(-3*atanh(a*x))
    F = 9*sqrt(a)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(-a*x + 1)**(sympy.S(3)/2) - a*x**2*(c - c/(a*x))**(sympy.S(3)/2)*(-a*x + 23)/((-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)) - 2*x*(c - c/(a*x))**(sympy.S(3)/2)*sqrt(-a*x + 1)/sqrt(a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_554():
    f = sqrt(c - c/(a*x))*exp(-3*atanh(a*x))
    F = x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1) + 8*x*sqrt(c - c/(a*x))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 7*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(sqrt(a)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_555():
    f = exp(-3*atanh(a*x))/sqrt(c - c/(a*x))
    F = -x*(-a*x + 1)/(sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)) - 5*sqrt(-a*x + 1)/(a*sqrt(c - c/(a*x))*sqrt(a*x + 1)) + 5*sqrt(-a*x + 1)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(3)/2)*sqrt(x)*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_556():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = -2*(-a*x + 1)**(sympy.S(3)/2)/(a*(c - c/(a*x))**(sympy.S(3)/2)*sqrt(a*x + 1)) + 3*(-a*x + 1)**(sympy.S(3)/2)*sqrt(a*x + 1)/(a**2*x*(c - c/(a*x))**(sympy.S(3)/2)) - 3*(-a*x + 1)**(sympy.S(3)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(5)/2)*x**(sympy.S(3)/2)*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_557():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))**(sympy.S(5)/2)
    F = (-a*x + 1)**(sympy.S(5)/2)/(a**2*x*(c - c/(a*x))**(sympy.S(5)/2)*sqrt(a*x + 1)) - 2*(-a*x + 1)**(sympy.S(5)/2)*sqrt(a*x + 1)/(a**3*x**2*(c - c/(a*x))**(sympy.S(5)/2)) + (-a*x + 1)**(sympy.S(5)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2)) + sqrt(2)*(-a*x + 1)**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(2*a**(sympy.S(7)/2)*x**(sympy.S(5)/2)*(c - c/(a*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_558():
    f = exp(-3*atanh(a*x))/(c - c/(a*x))**(sympy.S(7)/2)
    F = (-a*x + 1)**(sympy.S(5)/2)/(2*a**2*x*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)) - (-a*x + 1)**(sympy.S(7)/2)/(4*a**3*x**2*(c - c/(a*x))**(sympy.S(7)/2)*sqrt(a*x + 1)) + 7*(-a*x + 1)**(sympy.S(7)/2)*sqrt(a*x + 1)/(4*a**4*x**3*(c - c/(a*x))**(sympy.S(7)/2)) + (-a*x + 1)**(sympy.S(7)/2)*asinh(sqrt(a)*sqrt(x))/(a**(sympy.S(9)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2)) - 11*sqrt(2)*(-a*x + 1)**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(8*a**(sympy.S(9)/2)*x**(sympy.S(7)/2)*(c - c/(a*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_559():
    f = exp(3*atanh(a*x))/(x**3*(c - c/(a*x)))
    F = -8*a**2*(a*x + 1)/(3*c*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 4*a**2*(4*a*x + 3)/(3*c*sqrt(-a**2*x**2 + 1)) + 4*a**2*atanh(sqrt(-a**2*x**2 + 1))/c + a*sqrt(-a**2*x**2 + 1)/(c*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_560():
    f = x**m*sqrt(c - c/(a*x))*exp(atanh(a*x))
    F = 2*x**(m + 1)*sqrt(c - c/(a*x))*hyper((sympy.S(-1)/2, m + sympy.S.Half), (m + sympy.S(3)/2,), -a*x)/((2*m + 1)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_561():
    f = x**2*sqrt(c - c/(a*x))*exp(atanh(a*x))
    F = x**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(3*sqrt(-a*x + 1)) + x**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(12*a*sqrt(-a*x + 1)) - x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(8*a**2*sqrt(-a*x + 1)) + sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(8*a**(sympy.S(5)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_562():
    f = x*sqrt(c - c/(a*x))*exp(atanh(a*x))
    F = x**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(2*sqrt(-a*x + 1)) + x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(4*a*sqrt(-a*x + 1)) - sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(4*a**(sympy.S(3)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_563():
    f = sqrt(c - c/(a*x))*exp(atanh(a*x))
    F = x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1) + sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(sqrt(a)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_564():
    f = sqrt(c - c/(a*x))*exp(atanh(a*x))/x
    F = 2*sqrt(a)*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/sqrt(-a*x + 1) - 2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_565():
    f = sqrt(c - c/(a*x))*exp(atanh(a*x))/x**2
    F = -2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(3*x*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_566():
    f = sqrt(c - c/(a*x))*exp(atanh(a*x))/x**3
    F = 4*a*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(15*x*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(5*x**2*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_567():
    f = sqrt(c - c/(a*x))*exp(atanh(a*x))/x**4
    F = -16*a**2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(105*x*sqrt(-a*x + 1)) + 8*a*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(35*x**2*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(7*x**3*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_568():
    f = sqrt(c - c/(a*x))*exp(atanh(a*x))/x**5
    F = 32*a**3*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(315*x*sqrt(-a*x + 1)) - 16*a**2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(105*x**2*sqrt(-a*x + 1)) + 4*a*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(21*x**3*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(9*x**4*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_569():
    f = x**3*sqrt(c - c/(a*x))*exp(2*atanh(a*x))
    F = -x**4*sqrt(c - c/(a*x))/4 - 5*x**3*sqrt(c - c/(a*x))/(8*a) - 25*x**2*sqrt(c - c/(a*x))/(32*a**2) - 75*x*sqrt(c - c/(a*x))/(64*a**3) - 75*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_570():
    f = x**2*sqrt(c - c/(a*x))*exp(2*atanh(a*x))
    F = -x**3*sqrt(c - c/(a*x))/3 - 11*x**2*sqrt(c - c/(a*x))/(12*a) - 11*x*sqrt(c - c/(a*x))/(8*a**2) - 11*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_571():
    f = x*sqrt(c - c/(a*x))*exp(2*atanh(a*x))
    F = -x**2*sqrt(c - c/(a*x))/2 - 7*x*sqrt(c - c/(a*x))/(4*a) - 7*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_572():
    f = sqrt(c - c/(a*x))*exp(2*atanh(a*x))
    F = -x*sqrt(c - c/(a*x)) - 3*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_573():
    f = sqrt(c - c/(a*x))*exp(2*atanh(a*x))/x
    F = -2*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c)) - 2*sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_574():
    f = sqrt(c - c/(a*x))*exp(2*atanh(a*x))/x**2
    F = -4*a*sqrt(c - c/(a*x)) + 2*a*(c - c/(a*x))**(sympy.S(3)/2)/(3*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_575():
    f = sqrt(c - c/(a*x))*exp(2*atanh(a*x))/x**3
    F = -4*a**2*sqrt(c - c/(a*x)) + 2*a**2*(c - c/(a*x))**(sympy.S(3)/2)/c - 2*a**2*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_576():
    f = sqrt(c - c/(a*x))*exp(2*atanh(a*x))/x**4
    F = -4*a**3*sqrt(c - c/(a*x)) + 10*a**3*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) - 8*a**3*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2) + 2*a**3*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_577():
    f = sqrt(c - c/(a*x))*exp(2*atanh(a*x))/x**5
    F = -4*a**4*sqrt(c - c/(a*x)) + 14*a**4*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) - 18*a**4*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2) + 10*a**4*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3) - 2*a**4*(c - c/(a*x))**(sympy.S(9)/2)/(9*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_578():
    f = x**3*sqrt(c - c/(a*x))*exp(3*atanh(a*x))
    F = -x**3*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(4*a*sqrt(-a*x + 1)) - 11*x**2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(24*a**2*sqrt(-a*x + 1)) - 21*x*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(32*a**3*sqrt(-a*x + 1)) - 107*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(64*a**3*sqrt(-a*x + 1)) - 363*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(64*a**(sympy.S(7)/2)*sqrt(-a*x + 1)) + 4*sqrt(2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(a**(sympy.S(7)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_579():
    f = x**2*sqrt(c - c/(a*x))*exp(3*atanh(a*x))
    F = -x**2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(3*a*sqrt(-a*x + 1)) - 3*x*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(4*a**2*sqrt(-a*x + 1)) - 13*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(8*a**2*sqrt(-a*x + 1)) - 45*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(8*a**(sympy.S(5)/2)*sqrt(-a*x + 1)) + 4*sqrt(2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(a**(sympy.S(5)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_580():
    f = x*sqrt(c - c/(a*x))*exp(3*atanh(a*x))
    F = -x*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(2*a*sqrt(-a*x + 1)) - 7*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(4*a*sqrt(-a*x + 1)) - 23*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(4*a**(sympy.S(3)/2)*sqrt(-a*x + 1)) + 4*sqrt(2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(a**(sympy.S(3)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_581():
    f = sqrt(c - c/(a*x))*exp(3*atanh(a*x))
    F = -x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1) - 5*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(sqrt(a)*sqrt(-a*x + 1)) + 4*sqrt(2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/(sqrt(a)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_582():
    f = sqrt(c - c/(a*x))*exp(3*atanh(a*x))/x
    F = -2*sqrt(a)*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/sqrt(-a*x + 1) + 4*sqrt(2)*sqrt(a)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/sqrt(-a*x + 1) - 2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_583():
    f = sqrt(c - c/(a*x))*exp(3*atanh(a*x))/x**2
    F = 4*sqrt(2)*a**(sympy.S(3)/2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/sqrt(-a*x + 1) - 4*a*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1) - 2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(3*x*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_584():
    f = sqrt(c - c/(a*x))*exp(3*atanh(a*x))/x**3
    F = 4*sqrt(2)*a**(sympy.S(5)/2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/sqrt(-a*x + 1) - 4*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1) - 2*a*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(3)/2)/(3*x*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*(a*x + 1)**(sympy.S(5)/2)/(5*x**2*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_585():
    f = sqrt(c - c/(a*x))*exp(3*atanh(a*x))/x**4
    F = 4*sqrt(2)*a**(sympy.S(7)/2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/sqrt(-a*x + 1) - 104*a**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(21*sqrt(-a*x + 1)) - 32*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(21*x*sqrt(-a*x + 1)) - 6*a*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(7*x**2*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(7*x**3*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_586():
    f = sqrt(c - c/(a*x))*exp(3*atanh(a*x))/x**5
    F = 4*sqrt(2)*a**(sympy.S(9)/2)*sqrt(x)*sqrt(c - c/(a*x))*atanh(sqrt(2)*sqrt(a)*sqrt(x)/sqrt(a*x + 1))/sqrt(-a*x + 1) - 1576*a**4*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(315*sqrt(-a*x + 1)) - 472*a**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(315*x*sqrt(-a*x + 1)) - 92*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(105*x**2*sqrt(-a*x + 1)) - 38*a*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(63*x**3*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(9*x**4*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_587():
    f = x**m*sqrt(c - c/(a*x))*exp(-atanh(a*x))
    F = -x**(m + 1)*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/((m + 1)*(-a*x + 1)) + x**(m + 1)*sqrt(c - c/(a*x))*(4*m + 3)*hyper((sympy.S.Half, m + sympy.S.Half), (m + sympy.S(3)/2,), -a*x)/((m + 1)*(2*m + 1)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_588():
    f = x**2*sqrt(c - c/(a*x))*exp(-atanh(a*x))
    F = -x**3*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(-3*a*x + 3) + 11*x**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(12*a*sqrt(-a*x + 1)) - 11*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(8*a**2*sqrt(-a*x + 1)) + 11*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(8*a**(sympy.S(5)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_589():
    f = x*sqrt(c - c/(a*x))*exp(-atanh(a*x))
    F = -x**2*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(-2*a*x + 2) + 7*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(4*a*sqrt(-a*x + 1)) - 7*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(4*a**(sympy.S(3)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_590():
    f = sqrt(c - c/(a*x))*exp(-atanh(a*x))
    F = -x*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(-a*x + 1) + 3*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(sqrt(a)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_591():
    f = sqrt(c - c/(a*x))*exp(-atanh(a*x))/x
    F = -2*sqrt(a)*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/sqrt(-a*x + 1) - 2*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_592():
    f = sqrt(c - c/(a*x))*exp(-atanh(a*x))/x**2
    F = 10*a*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(3*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(3*x*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_593():
    f = sqrt(c - c/(a*x))*exp(-atanh(a*x))/x**3
    F = -12*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(5*sqrt(-a*x + 1)) + 6*a*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(5*x*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(5*x**2*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_594():
    f = sqrt(c - c/(a*x))*exp(-atanh(a*x))/x**4
    F = 208*a**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(105*sqrt(-a*x + 1)) - 104*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(105*x*sqrt(-a*x + 1)) + 26*a*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(35*x**2*sqrt(-a*x + 1)) - 2*sqrt(c - c/(a*x))*sqrt(-a**2*x**2 + 1)/(7*x**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_595():
    f = x**3*sqrt(c - c/(a*x))*exp(-2*atanh(a*x))
    F = -x**4*sqrt(c - c/(a*x))/4 + 17*x**3*sqrt(c - c/(a*x))/(24*a) - 107*x**2*sqrt(c - c/(a*x))/(96*a**2) + 149*x*sqrt(c - c/(a*x))/(64*a**3) - 363*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(64*a**4) + 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_596():
    f = x**2*sqrt(c - c/(a*x))*exp(-2*atanh(a*x))
    F = -x**3*sqrt(c - c/(a*x))/3 + 13*x**2*sqrt(c - c/(a*x))/(12*a) - 19*x*sqrt(c - c/(a*x))/(8*a**2) + 45*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(8*a**3) - 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_597():
    f = x*sqrt(c - c/(a*x))*exp(-2*atanh(a*x))
    F = -x**2*sqrt(c - c/(a*x))/2 + 9*x*sqrt(c - c/(a*x))/(4*a) - 23*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/(4*a**2) + 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_598():
    f = sqrt(c - c/(a*x))*exp(-2*atanh(a*x))
    F = -x*sqrt(c - c/(a*x)) + 5*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c))/a - 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c)))/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_599():
    f = sqrt(c - c/(a*x))*exp(-2*atanh(a*x))/x
    F = -2*sqrt(c)*atanh(sqrt(c - c/(a*x))/sqrt(c)) + 4*sqrt(2)*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) - 2*sqrt(c - c/(a*x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_600():
    f = sqrt(c - c/(a*x))*exp(-2*atanh(a*x))/x**2
    F = -4*sqrt(2)*a*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) + 4*a*sqrt(c - c/(a*x)) + 2*a*(c - c/(a*x))**(sympy.S(3)/2)/(3*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_601():
    f = sqrt(c - c/(a*x))*exp(-2*atanh(a*x))/x**3
    F = 4*sqrt(2)*a**2*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) - 4*a**2*sqrt(c - c/(a*x)) - 2*a**2*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) - 2*a**2*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_602():
    f = sqrt(c - c/(a*x))*exp(-2*atanh(a*x))/x**4
    F = -4*sqrt(2)*a**3*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) + 4*a**3*sqrt(c - c/(a*x)) + 2*a**3*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) + 2*a**3*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_603():
    f = sqrt(c - c/(a*x))*exp(-2*atanh(a*x))/x**5
    F = 4*sqrt(2)*a**4*sqrt(c)*atanh(sqrt(2)*sqrt(c - c/(a*x))/(2*sqrt(c))) - 4*a**4*sqrt(c - c/(a*x)) - 2*a**4*(c - c/(a*x))**(sympy.S(3)/2)/(3*c) - 2*a**4*(c - c/(a*x))**(sympy.S(5)/2)/(5*c**2) + 2*a**4*(c - c/(a*x))**(sympy.S(7)/2)/(7*c**3) - 2*a**4*(c - c/(a*x))**(sympy.S(9)/2)/(9*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_604():
    f = x**3*sqrt(c - c/(a*x))*exp(-3*atanh(a*x))
    F = x**4*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(4*sqrt(-a*x + 1)) + 8*x**4*sqrt(c - c/(a*x))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 223*x**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(24*a*sqrt(-a*x + 1)) + 1115*x**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(96*a**2*sqrt(-a*x + 1)) - 1115*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(64*a**3*sqrt(-a*x + 1)) + 1115*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(64*a**(sympy.S(7)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_605():
    f = x**2*sqrt(c - c/(a*x))*exp(-3*atanh(a*x))
    F = x**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(3*sqrt(-a*x + 1)) + 8*x**3*sqrt(c - c/(a*x))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 119*x**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(12*a*sqrt(-a*x + 1)) + 119*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(8*a**2*sqrt(-a*x + 1)) - 119*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(8*a**(sympy.S(5)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_606():
    f = x*sqrt(c - c/(a*x))*exp(-3*atanh(a*x))
    F = x**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(2*sqrt(-a*x + 1)) + 8*x**2*sqrt(c - c/(a*x))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 47*x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(4*a*sqrt(-a*x + 1)) + 47*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(4*a**(sympy.S(3)/2)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_607():
    f = sqrt(c - c/(a*x))*exp(-3*atanh(a*x))
    F = x*sqrt(c - c/(a*x))*sqrt(a*x + 1)/sqrt(-a*x + 1) + 8*x*sqrt(c - c/(a*x))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 7*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/(sqrt(a)*sqrt(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_608():
    f = sqrt(c - c/(a*x))*exp(-3*atanh(a*x))/x
    F = 2*sqrt(a)*sqrt(x)*sqrt(c - c/(a*x))*asinh(sqrt(a)*sqrt(x))/sqrt(-a*x + 1) - 10*a*x*sqrt(c - c/(a*x))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*sqrt(c - c/(a*x))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_609():
    f = sqrt(c - c/(a*x))*exp(-3*atanh(a*x))/x**2
    F = 46*a**2*x*sqrt(c - c/(a*x))/(3*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 20*a*sqrt(c - c/(a*x))/(3*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*sqrt(c - c/(a*x))/(3*x*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_610():
    f = sqrt(c - c/(a*x))*exp(-3*atanh(a*x))/x**3
    F = -316*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(15*sqrt(-a*x + 1)) + 158*a**2*sqrt(c - c/(a*x))/(15*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 32*a*sqrt(c - c/(a*x))/(15*x*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*sqrt(c - c/(a*x))/(5*x**2*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_611():
    f = sqrt(c - c/(a*x))*exp(-3*atanh(a*x))/x**4
    F = 2672*a**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(105*sqrt(-a*x + 1)) - 1336*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(105*x*sqrt(-a*x + 1)) + 334*a**2*sqrt(c - c/(a*x))/(35*x*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 44*a*sqrt(c - c/(a*x))/(35*x**2*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*sqrt(c - c/(a*x))/(7*x**3*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_612():
    f = sqrt(c - c/(a*x))*exp(-3*atanh(a*x))/x**5
    F = -1312*a**4*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(45*sqrt(-a*x + 1)) + 656*a**3*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(45*x*sqrt(-a*x + 1)) - 164*a**2*sqrt(c - c/(a*x))*sqrt(a*x + 1)/(15*x**2*sqrt(-a*x + 1)) + 82*a**2*sqrt(c - c/(a*x))/(9*x**2*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 8*a*sqrt(c - c/(a*x))/(9*x**3*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*sqrt(c - c/(a*x))/(9*x**4*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_613():
    f = (c - c/(a*x))**p*exp(n*atanh(a*x))
    F = x*(c - c/(a*x))**p*appellf1(1 - p, -n/2, n/2 - p, 2 - p, -a*x, a*x)/((1 - p)*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_614():
    f = (c - c/(a*x))**p*exp(-2*p*atanh(a*x))
    F = x*(c - c/(a*x))**p*appellf1(1 - p, -2*p, p, 2 - p, a*x, -a*x)/((1 - p)*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_615():
    f = (c - c/(a*x))**p*exp(2*p*atanh(a*x))
    F = x*(c - c/(a*x))**p*hyper((-p, 1 - p), (2 - p,), -a*x)/((1 - p)*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_616():
    f = (c - c/(a*x))**2*exp(n*atanh(a*x))
    F = 2**(n/2)*c**2*(-a*x + 1)**(2 - n/2)*hyper((1 - n/2, 2 - n/2), (3 - n/2,), -a*x/2 + sympy.S.Half)/(a*(4 - n)) + 4*c**2*(a*x + 1)**(n/2)*hyper((2, n/2), (n/2 + 1,), (a*x + 1)/(-a*x + 1))/(a*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_617():
    f = (c - c/(a*x))*exp(n*atanh(a*x))
    F = 2**(n/2)*c*(1 - n)*(-a*x + 1)**(2 - n/2)*hyper((1 - n/2, 2 - n/2), (3 - n/2,), -a*x/2 + sympy.S.Half)/(a*(2 - n)*(4 - n)) - 2*c*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)*hyper((1, n/2 - 1), (n/2,), (a*x + 1)/(-a*x + 1))/(a*(2 - n)) + c*(-a*x + 1)**(2 - n/2)*(a*x + 1)**(n/2 - 1)/(a*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_618():
    f = exp(n*atanh(a*x))/(c - c/(a*x))
    F = -2**(n/2 + 1)*(n + 1)*(-a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*n*(2 - n)) - (a*x + 1)**(n/2 + 1)/(a*c*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_619():
    f = exp(n*atanh(a*x))/(c - c/(a*x))**2
    F = -2**(n/2 + 1)*(n + 2)*hyper((-n/2, -n/2), (1 - n/2,), -a*x/2 + sympy.S.Half)/(a*c**2*n*(-a*x + 1)**(n/2)) - x*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 + 1)/c**2 + (n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 + 1)/(a*c**2*(n + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_620():
    f = (c - c/(a*x))**(sympy.S(3)/2)*exp(n*atanh(a*x))
    F = -2*x*(c - c/(a*x))**(sympy.S(3)/2)*appellf1(sympy.S(-1)/2, -n/2, n/2 + sympy.S(-3)/2, sympy.S.Half, -a*x, a*x)/(-a*x + 1)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_621():
    f = sqrt(c - c/(a*x))*exp(n*atanh(a*x))
    F = 2*x*sqrt(c - c/(a*x))*appellf1(sympy.S.Half, -n/2, n/2 + sympy.S(-1)/2, sympy.S(3)/2, -a*x, a*x)/sqrt(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_622():
    f = exp(n*atanh(a*x))/sqrt(c - c/(a*x))
    F = 2*x*sqrt(-a*x + 1)*appellf1(sympy.S(3)/2, -n/2, n/2 + sympy.S.Half, sympy.S(5)/2, -a*x, a*x)/(3*sqrt(c - c/(a*x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_623():
    f = exp(n*atanh(a*x))/(c - c/(a*x))**(sympy.S(3)/2)
    F = 2*x*(-a*x + 1)**(sympy.S(3)/2)*appellf1(sympy.S(5)/2, -n/2, n/2 + sympy.S(3)/2, sympy.S(7)/2, -a*x, a*x)/(5*(c - c/(a*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_624():
    f = (c - c/(a**2*x**2))**4*exp(atanh(a*x))
    F = c**4*asin(a*x)/a + 35*c**4*atanh(sqrt(-a**2*x**2 + 1))/(16*a) + c**4*(-35*a*x + 16)*sqrt(-a**2*x**2 + 1)/(16*a**2*x) - c**4*(35*a*x + 16)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(48*a**4*x**3) + c**4*(35*a*x + 24)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(120*a**6*x**5) - c**4*(7*a*x + 6)*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(42*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_625():
    f = (c - c/(a**2*x**2))**3*exp(atanh(a*x))
    F = c**3*asin(a*x)/a + 15*c**3*atanh(sqrt(-a**2*x**2 + 1))/(8*a) + c**3*(-15*a*x + 8)*sqrt(-a**2*x**2 + 1)/(8*a**2*x) - c**3*(15*a*x + 8)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(24*a**4*x**3) + c**3*(5*a*x + 4)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(20*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_626():
    f = (c - c/(a**2*x**2))**2*exp(atanh(a*x))
    F = c**2*asin(a*x)/a + 3*c**2*atanh(sqrt(-a**2*x**2 + 1))/(2*a) + c**2*(-3*a*x + 2)*sqrt(-a**2*x**2 + 1)/(2*a**2*x) - c**2*(3*a*x + 2)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_627():
    f = (c - c/(a**2*x**2))*exp(atanh(a*x))
    F = c*asin(a*x)/a + c*atanh(sqrt(-a**2*x**2 + 1))/a + c*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_628():
    f = exp(atanh(a*x))/(c - c/(a**2*x**2))
    F = -(a*x + 1)/(a*c*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*x**2 + 1)/(a*c) + asin(a*x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_629():
    f = exp(atanh(a*x))/(c - c/(a**2*x**2))**2
    F = a**2*x**3*(a*x + 1)/(3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(4*a*x + 3)/(3*c**2*sqrt(-a**2*x**2 + 1)) - 8*sqrt(-a**2*x**2 + 1)/(3*a*c**2) + asin(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_630():
    f = exp(atanh(a*x))/(c - c/(a**2*x**2))**3
    F = -a**4*x**5*(a*x + 1)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + a**2*x**3*(6*a*x + 5)/(15*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(8*a*x + 5)/(5*c**3*sqrt(-a**2*x**2 + 1)) - 16*sqrt(-a**2*x**2 + 1)/(5*a*c**3) + asin(a*x)/(a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_631():
    f = exp(atanh(a*x))/(c - c/(a**2*x**2))**4
    F = a**6*x**7*(a*x + 1)/(7*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - a**4*x**5*(8*a*x + 7)/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + a**2*x**3*(48*a*x + 35)/(105*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(64*a*x + 35)/(35*c**4*sqrt(-a**2*x**2 + 1)) - 128*sqrt(-a**2*x**2 + 1)/(35*a*c**4) + asin(a*x)/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_632():
    f = (c - c/(a**2*x**2))**5*exp(2*atanh(a*x))
    F = -c**5*x - 2*c**5*log(x)/a - 3*c**5/(a**2*x) - 4*c**5/(a**3*x**2) + 2*c**5/(3*a**4*x**3) + 3*c**5/(a**5*x**4) + 2*c**5/(5*a**6*x**5) - 4*c**5/(3*a**7*x**6) - 3*c**5/(7*a**8*x**7) + c**5/(4*a**9*x**8) + c**5/(9*a**10*x**9)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_633():
    f = (c - c/(a**2*x**2))**4*exp(2*atanh(a*x))
    F = -c**4*x - 2*c**4*log(x)/a - 2*c**4/(a**2*x) - 3*c**4/(a**3*x**2) + 3*c**4/(2*a**5*x**4) + 2*c**4/(5*a**6*x**5) - c**4/(3*a**7*x**6) - c**4/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_634():
    f = (c - c/(a**2*x**2))**3*exp(2*atanh(a*x))
    F = -c**3*x - 2*c**3*log(x)/a - c**3/(a**2*x) - 2*c**3/(a**3*x**2) - c**3/(3*a**4*x**3) + c**3/(2*a**5*x**4) + c**3/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_635():
    f = (c - c/(a**2*x**2))**2*exp(2*atanh(a*x))
    F = -c**2*x - 2*c**2*log(x)/a - c**2/(a**3*x**2) - c**2/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_636():
    f = (c - c/(a**2*x**2))*exp(2*atanh(a*x))
    F = -c*x - 2*c*log(x)/a + c/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_637():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))
    F = -x/c - 2*log(-a*x + 1)/(a*c) - 1/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_638():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))**2
    F = -x/c**2 - 17*log(-a*x + 1)/(8*a*c**2) + log(a*x + 1)/(8*a*c**2) - 7/(4*a*c**2*(-a*x + 1)) + 1/(4*a*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_639():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))**3
    F = -x/c**3 - 9*log(-a*x + 1)/(4*a*c**3) + log(a*x + 1)/(4*a*c**3) + 1/(16*a*c**3*(a*x + 1)) - 39/(16*a*c**3*(-a*x + 1)) + 5/(8*a*c**3*(-a*x + 1)**2) - 1/(12*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_640():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))**4
    F = -x/c**4 - 303*log(-a*x + 1)/(128*a*c**4) + 47*log(a*x + 1)/(128*a*c**4) + 11/(64*a*c**4*(a*x + 1)) - 1/(64*a*c**4*(a*x + 1)**2) - 99/(32*a*c**4*(-a*x + 1)) + 35/(32*a*c**4*(-a*x + 1)**2) - 13/(48*a*c**4*(-a*x + 1)**3) + 1/(32*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_641():
    f = (c - c/(a**2*x**2))**4*exp(3*atanh(a*x))
    F = -3*c**4*asin(a*x)/a - 15*c**4*atanh(sqrt(-a**2*x**2 + 1))/(16*a) - 3*c**4*(-5*a*x + 16)*sqrt(-a**2*x**2 + 1)/(16*a**2*x) + c**4*(5*a*x + 16)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(16*a**4*x**3) - c**4*(5*a*x + 24)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(40*a**6*x**5) - c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(2*a**7*x**6) - c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_642():
    f = (c - c/(a**2*x**2))**3*exp(3*atanh(a*x))
    F = -3*c**3*asin(a*x)/a - 3*c**3*atanh(sqrt(-a**2*x**2 + 1))/(8*a) - 3*c**3*(-a*x + 8)*sqrt(-a**2*x**2 + 1)/(8*a**2*x) + c**3*(a*x + 8)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(8*a**4*x**3) + 3*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(4*a**5*x**4) + c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_643():
    f = (c - c/(a**2*x**2))**2*exp(3*atanh(a*x))
    F = -3*c**2*asin(a*x)/a + c**2*atanh(sqrt(-a**2*x**2 + 1))/(2*a) - c**2*(a*x + 6)*sqrt(-a**2*x**2 + 1)/(2*a**2*x) - 3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**3*x**2) - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_644():
    f = (c - c/(a**2*x**2))*exp(3*atanh(a*x))
    F = c*sqrt(-a**2*x**2 + 1)/a - 3*c*asin(a*x)/a + 3*c*atanh(sqrt(-a**2*x**2 + 1))/a + c*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_645():
    f = exp(3*atanh(a*x))/(c - c/(a**2*x**2))
    F = -(a*x + 1)**3/(3*a*c*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 2*(a*x + 1)**2/(a*c*sqrt(-a**2*x**2 + 1)) + 3*sqrt(-a**2*x**2 + 1)/(a*c) - 3*asin(a*x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_646():
    f = exp(3*atanh(a*x))/(c - c/(a**2*x**2))**2
    F = (a*x + 1)**3/(5*a*c**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 6*(a*x + 1)**2/(5*a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (24*a*x + 24)/(5*a*c**2*sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(a*c**2) - 3*asin(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_647():
    f = exp(3*atanh(a*x))/(c - c/(a**2*x**2))**3
    F = -(a*x + 1)**3/(7*a*c**3*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + 38*(a*x + 1)**2/(35*a*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - (137*a*x + 137)/(35*a*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (181*a*x + 245)/(35*a*c**3*sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(a*c**3) - 3*asin(a*x)/(a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_648():
    f = exp(3*atanh(a*x))/(c - c/(a**2*x**2))**4
    F = (a*x + 1)**3/(9*a*c**4*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - 22*(a*x + 1)**2/(21*a*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + (478*a*x + 478)/(105*a*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - (1658*a*x + 2310)/(315*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (1724*a*x + 2520)/(315*a*c**4*sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(a*c**4) - 3*asin(a*x)/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_649():
    f = (c - c/(a**2*x**2))**5*exp(4*atanh(a*x))
    F = c**5*x + 4*c**5*log(x)/a - 3*c**5/(a**2*x) + 4*c**5/(a**3*x**2) + 14*c**5/(3*a**4*x**3) - 14*c**5/(5*a**6*x**5) - 4*c**5/(3*a**7*x**6) + 3*c**5/(7*a**8*x**7) + c**5/(2*a**9*x**8) + c**5/(9*a**10*x**9)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_650():
    f = (c - c/(a**2*x**2))**4*exp(4*atanh(a*x))
    F = c**4*x + 4*c**4*log(x)/a - 4*c**4/(a**2*x) + 2*c**4/(a**3*x**2) + 10*c**4/(3*a**4*x**3) + c**4/(a**5*x**4) - 4*c**4/(5*a**6*x**5) - 2*c**4/(3*a**7*x**6) - c**4/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_651():
    f = (c - c/(a**2*x**2))**3*exp(4*atanh(a*x))
    F = c**3*x + 4*c**3*log(x)/a - 5*c**3/(a**2*x) + 5*c**3/(3*a**4*x**3) + c**3/(a**5*x**4) + c**3/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_652():
    f = (c - c/(a**2*x**2))**2*exp(4*atanh(a*x))
    F = c**2*x + 4*c**2*log(x)/a - 6*c**2/(a**2*x) - 2*c**2/(a**3*x**2) - c**2/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_653():
    f = (c - c/(a**2*x**2))*exp(4*atanh(a*x))
    F = c*x - 4*c*log(x)/a + 8*c*log(-a*x + 1)/a + c/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_654():
    f = exp(4*atanh(a*x))/(c - c/(a**2*x**2))
    F = x/c + 4*log(-a*x + 1)/(a*c) + 5/(a*c*(-a*x + 1)) - 1/(a*c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_655():
    f = exp(4*atanh(a*x))/(c - c/(a**2*x**2))**2
    F = x/c**2 + 4*log(-a*x + 1)/(a*c**2) + 6/(a*c**2*(-a*x + 1)) - 2/(a*c**2*(-a*x + 1)**2) + 1/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_656():
    f = exp(4*atanh(a*x))/(c - c/(a**2*x**2))**3
    F = x/c**3 + 129*log(-a*x + 1)/(32*a*c**3) - log(a*x + 1)/(32*a*c**3) + 111/(16*a*c**3*(-a*x + 1)) - 49/(16*a*c**3*(-a*x + 1)**2) + 11/(12*a*c**3*(-a*x + 1)**3) - 1/(8*a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_657():
    f = exp(4*atanh(a*x))/(c - c/(a**2*x**2))**4
    F = x/c**4 + 261*log(-a*x + 1)/(64*a*c**4) - 5*log(a*x + 1)/(64*a*c**4) - 1/(64*a*c**4*(a*x + 1)) + 501/(64*a*c**4*(-a*x + 1)) - 67/(16*a*c**4*(-a*x + 1)**2) + 83/(48*a*c**4*(-a*x + 1)**3) - 7/(16*a*c**4*(-a*x + 1)**4) + 1/(20*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_658():
    f = (c - c/(a**2*x**2))**4*exp(-atanh(a*x))
    F = c**4*asin(a*x)/a - 35*c**4*atanh(sqrt(-a**2*x**2 + 1))/(16*a) + c**4*(35*a*x + 16)*sqrt(-a**2*x**2 + 1)/(16*a**2*x) - c**4*(-35*a*x + 16)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(48*a**4*x**3) + c**4*(-35*a*x + 24)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(120*a**6*x**5) - c**4*(-7*a*x + 6)*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(42*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_659():
    f = (c - c/(a**2*x**2))**3*exp(-atanh(a*x))
    F = c**3*asin(a*x)/a - 15*c**3*atanh(sqrt(-a**2*x**2 + 1))/(8*a) + c**3*(15*a*x + 8)*sqrt(-a**2*x**2 + 1)/(8*a**2*x) - c**3*(-15*a*x + 8)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(24*a**4*x**3) + c**3*(-5*a*x + 4)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(20*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_660():
    f = (c - c/(a**2*x**2))**2*exp(-atanh(a*x))
    F = c**2*asin(a*x)/a - 3*c**2*atanh(sqrt(-a**2*x**2 + 1))/(2*a) + c**2*(3*a*x + 2)*sqrt(-a**2*x**2 + 1)/(2*a**2*x) - c**2*(-3*a*x + 2)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(6*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_661():
    f = (c - c/(a**2*x**2))*exp(-atanh(a*x))
    F = c*asin(a*x)/a - c*atanh(sqrt(-a**2*x**2 + 1))/a + c*(a*x + 1)*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_662():
    f = exp(-atanh(a*x))/(c - c/(a**2*x**2))
    F = (-a*x + 1)/(a*c*sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(a*c) + asin(a*x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_663():
    f = exp(-atanh(a*x))/(c - c/(a**2*x**2))**2
    F = a**2*x**3*(-a*x + 1)/(3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(-4*a*x + 3)/(3*c**2*sqrt(-a**2*x**2 + 1)) + 8*sqrt(-a**2*x**2 + 1)/(3*a*c**2) + asin(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_664():
    f = exp(-atanh(a*x))/(c - c/(a**2*x**2))**3
    F = -a**4*x**5*(-a*x + 1)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + a**2*x**3*(-6*a*x + 5)/(15*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(-8*a*x + 5)/(5*c**3*sqrt(-a**2*x**2 + 1)) + 16*sqrt(-a**2*x**2 + 1)/(5*a*c**3) + asin(a*x)/(a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_665():
    f = exp(-atanh(a*x))/(c - c/(a**2*x**2))**4
    F = a**6*x**7*(-a*x + 1)/(7*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - a**4*x**5*(-8*a*x + 7)/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + a**2*x**3*(-48*a*x + 35)/(105*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(-64*a*x + 35)/(35*c**4*sqrt(-a**2*x**2 + 1)) + 128*sqrt(-a**2*x**2 + 1)/(35*a*c**4) + asin(a*x)/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_666():
    f = (c - c/(a**2*x**2))**4*exp(-2*atanh(a*x))
    F = -c**4*x + 2*c**4*log(x)/a - 2*c**4/(a**2*x) + 3*c**4/(a**3*x**2) - 3*c**4/(2*a**5*x**4) + 2*c**4/(5*a**6*x**5) + c**4/(3*a**7*x**6) - c**4/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_667():
    f = (c - c/(a**2*x**2))**3*exp(-2*atanh(a*x))
    F = -c**3*x + 2*c**3*log(x)/a - c**3/(a**2*x) + 2*c**3/(a**3*x**2) - c**3/(3*a**4*x**3) - c**3/(2*a**5*x**4) + c**3/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_668():
    f = (c - c/(a**2*x**2))**2*exp(-2*atanh(a*x))
    F = -c**2*x + 2*c**2*log(x)/a + c**2/(a**3*x**2) - c**2/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_669():
    f = (c - c/(a**2*x**2))*exp(-2*atanh(a*x))
    F = -c*x + 2*c*log(x)/a + c/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_670():
    f = exp(-2*atanh(a*x))/(c - c/(a**2*x**2))
    F = -x/c + 2*log(a*x + 1)/(a*c) + 1/(a*c*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_671():
    f = exp(-2*atanh(a*x))/(c - c/(a**2*x**2))**2
    F = -x/c**2 - log(-a*x + 1)/(8*a*c**2) + 17*log(a*x + 1)/(8*a*c**2) + 7/(4*a*c**2*(a*x + 1)) - 1/(4*a*c**2*(a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_672():
    f = exp(-2*atanh(a*x))/(c - c/(a**2*x**2))**3
    F = -x/c**3 - log(-a*x + 1)/(4*a*c**3) + 9*log(a*x + 1)/(4*a*c**3) + 39/(16*a*c**3*(a*x + 1)) - 5/(8*a*c**3*(a*x + 1)**2) + 1/(12*a*c**3*(a*x + 1)**3) - 1/(16*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_673():
    f = exp(-2*atanh(a*x))/(c - c/(a**2*x**2))**4
    F = -x/c**4 - 47*log(-a*x + 1)/(128*a*c**4) + 303*log(a*x + 1)/(128*a*c**4) + 99/(32*a*c**4*(a*x + 1)) - 35/(32*a*c**4*(a*x + 1)**2) + 13/(48*a*c**4*(a*x + 1)**3) - 1/(32*a*c**4*(a*x + 1)**4) - 11/(64*a*c**4*(-a*x + 1)) + 1/(64*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_674():
    f = (c - c/(a**2*x**2))**4*exp(-3*atanh(a*x))
    F = -3*c**4*asin(a*x)/a + 15*c**4*atanh(sqrt(-a**2*x**2 + 1))/(16*a) - 3*c**4*(5*a*x + 16)*sqrt(-a**2*x**2 + 1)/(16*a**2*x) + c**4*(-5*a*x + 16)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(16*a**4*x**3) - c**4*(-5*a*x + 24)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(40*a**6*x**5) + c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(2*a**7*x**6) - c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(7*a**8*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_675():
    f = (c - c/(a**2*x**2))**3*exp(-3*atanh(a*x))
    F = -3*c**3*asin(a*x)/a + 3*c**3*atanh(sqrt(-a**2*x**2 + 1))/(8*a) - 3*c**3*(a*x + 8)*sqrt(-a**2*x**2 + 1)/(8*a**2*x) + c**3*(-a*x + 8)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(8*a**4*x**3) - 3*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(4*a**5*x**4) + c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a**6*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_676():
    f = (c - c/(a**2*x**2))**2*exp(-3*atanh(a*x))
    F = -3*c**2*asin(a*x)/a - c**2*atanh(sqrt(-a**2*x**2 + 1))/(2*a) - c**2*(-a*x + 6)*sqrt(-a**2*x**2 + 1)/(2*a**2*x) + 3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**3*x**2) - c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_677():
    f = (c - c/(a**2*x**2))*exp(-3*atanh(a*x))
    F = -c*sqrt(-a**2*x**2 + 1)/a - 3*c*asin(a*x)/a - 3*c*atanh(sqrt(-a**2*x**2 + 1))/a + c*sqrt(-a**2*x**2 + 1)/(a**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_678():
    f = exp(-3*atanh(a*x))/(c - c/(a**2*x**2))
    F = (-a*x + 1)**3/(3*a*c*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 2*(-a*x + 1)**2/(a*c*sqrt(-a**2*x**2 + 1)) - 3*sqrt(-a**2*x**2 + 1)/(a*c) - 3*asin(a*x)/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_679():
    f = exp(-3*atanh(a*x))/(c - c/(a**2*x**2))**2
    F = -(-24*a*x + 24)/(5*a*c**2*sqrt(-a**2*x**2 + 1)) - (-a*x + 1)**3/(5*a*c**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 6*(-a*x + 1)**2/(5*a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - sqrt(-a**2*x**2 + 1)/(a*c**2) - 3*asin(a*x)/(a*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_680():
    f = exp(-3*atanh(a*x))/(c - c/(a**2*x**2))**3
    F = -(-181*a*x + 245)/(35*a*c**3*sqrt(-a**2*x**2 + 1)) + (-137*a*x + 137)/(35*a*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (-a*x + 1)**3/(7*a*c**3*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - 38*(-a*x + 1)**2/(35*a*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - sqrt(-a**2*x**2 + 1)/(a*c**3) - 3*asin(a*x)/(a*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_681():
    f = exp(-3*atanh(a*x))/(c - c/(a**2*x**2))**4
    F = -(-1724*a*x + 2520)/(315*a*c**4*sqrt(-a**2*x**2 + 1)) + (-1658*a*x + 2310)/(315*a*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - (-478*a*x + 478)/(105*a*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - (-a*x + 1)**3/(9*a*c**4*(-a**2*x**2 + 1)**(sympy.S(9)/2)) + 22*(-a*x + 1)**2/(21*a*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - sqrt(-a**2*x**2 + 1)/(a*c**4) - 3*asin(a*x)/(a*c**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_682():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(atanh(a*x))
    F = a**9*x**10*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + 4*a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + 2*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 2*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(9)/2)) + 4*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(9)/2)) + 2*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(3*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(7*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(8*(-a**2*x**2 + 1)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_683():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(atanh(a*x))
    F = -a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) - a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(7)/2) - 3*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) - 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + 3*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(6*(-a**2*x**2 + 1)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_684():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(atanh(a*x))
    F = a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_685():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(atanh(a*x))
    F = -a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - x*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_686():
    f = sqrt(c - c/(a**2*x**2))*exp(atanh(a*x))
    F = a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_687():
    f = exp(atanh(a*x))/sqrt(c - c/(a**2*x**2))
    F = -sqrt(-a**2*x**2 + 1)/(a*sqrt(c - c/(a**2*x**2))) - sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_688():
    f = exp(atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + 5*(-a**2*x**2 + 1)**(sympy.S(3)/2)*log(-a*x + 1)/(4*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - (-a**2*x**2 + 1)**(sympy.S(3)/2)*log(a*x + 1)/(4*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_689():
    f = exp(atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = -(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 23*(-a**2*x**2 + 1)**(sympy.S(5)/2)*log(-a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 7*(-a**2*x**2 + 1)**(sympy.S(5)/2)*log(a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)) - (-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_690():
    f = exp(atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = (-a**2*x**2 + 1)**(sympy.S(7)/2)/(a**7*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 51*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(-a*x + 1)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 19*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(a*x + 1)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 5*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(16*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)**2) + 3*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(2*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)) - 11*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)**2) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(24*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_691():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(2*atanh(a*x))
    F = -501*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(128*(-a*x + 1)**4*(a*x + 1)**4) + 2*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(9)/2)*(a*x + 1)**(sympy.S(9)/2)) + 245*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(128*(-a*x + 1)**(sympy.S(9)/2)*(a*x + 1)**(sympy.S(9)/2)) + 373*a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(192*(-a*x + 1)**4*(a*x + 1)**3) + 501*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(640*(-a*x + 1)**4*(a*x + 1)**2) + 661*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(1680*(-a*x + 1)**4*(a*x + 1)) + 295*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(1344*(-a*x + 1)**4) - 127*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(9)/2)*(a*x + 1)/(420*(-a*x + 1)**4) + 71*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(9)/2)*(a*x + 1)/(336*(-a*x + 1)**3) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(9)/2)*(a*x + 1)/(28*(-a*x + 1)**2) - x*(c - c/(a**2*x**2))**(sympy.S(9)/2)*(a*x + 1)/(-8*a*x + 8)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_692():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(2*atanh(a*x))
    F = 57*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(16*(-a*x + 1)**3*(a*x + 1)**3) - 2*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) - 25*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(16*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) - 41*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(24*(-a*x + 1)**3*(a*x + 1)**2) - 57*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(80*(-a*x + 1)**3*(a*x + 1)) - 11*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(30*(-a*x + 1)**3) + 13*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)/(40*(-a*x + 1)**3) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)/(15*(-a*x + 1)**2) - x*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)/(-6*a*x + 6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_693():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = -25*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(a*x + 1)**2) + 2*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) + 9*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) + 17*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(12*(-a*x + 1)**2*(a*x + 1)) + 5*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)/(6*(-a*x + 1)**2) - x*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)/(-4*a*x + 4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_694():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = -2*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) + 5*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(a*x + 1)) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a*x + 1) - x*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(a*x + 1)/(-2*a*x + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_695():
    f = sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))
    F = -x*sqrt(c - c/(a**2*x**2)) + 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_696():
    f = exp(2*atanh(a*x))/sqrt(c - c/(a**2*x**2))
    F = (-2*a*x + 2)*(a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))) - 2*sqrt(-a*x + 1)*sqrt(a*x + 1)*asin(a*x)/(a**2*x*sqrt(c - c/(a**2*x**2))) + (a*x + 1)**2/(a**2*x*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_697():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = (a*x + 1)**2/(3*a**2*x*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - (-4*a*x + 10)*(-a*x + 1)*(a*x + 1)**2/(3*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + 2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)*asin(a*x)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_698():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = (a*x + 1)**2/(5*a**2*x*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - (-2*a*x + 2)*(a*x + 1)**2/(3*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 58*(-a*x + 1)**2*(a*x + 1)**2/(15*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 2*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)*asin(a*x)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 2*(-a*x + 1)**3*(a*x + 1)**2*(43*a*x + 28)/(15*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_699():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = (a*x + 1)**2/(7*a**2*x*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - (-2*a*x + 2)*(a*x + 1)**2/(5*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 124*(-a*x + 1)**2*(a*x + 1)**2/(105*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 782*(-a*x + 1)**3*(a*x + 1)**2/(105*a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 142*(-a*x + 1)**4*(a*x + 1)**2/(35*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 2*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)*asin(a*x)/(a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 2*(-a*x + 1)**4*(a*x + 1)**3*(107*a*x + 72)/(35*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_700():
    f = exp(2*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(9)/2)
    F = (a*x + 1)**2/(9*a**2*x*(c - c/(a**2*x**2))**(sympy.S(9)/2)) - (-2*a*x + 2)*(a*x + 1)**2/(7*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(9)/2)) + 214*(-a*x + 1)**2*(a*x + 1)**2/(315*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(9)/2)) - 646*(-a*x + 1)**3*(a*x + 1)**2/(315*a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(9)/2)) + 302*(-a*x + 1)**4*(a*x + 1)**2/(21*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(9)/2)) + 2458*(-a*x + 1)**5*(a*x + 1)**2/(315*a**7*x**6*(c - c/(a**2*x**2))**(sympy.S(9)/2)) + 1334*(-a*x + 1)**5*(a*x + 1)**3/(315*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(9)/2)) - 2*(-a*x + 1)**(sympy.S(9)/2)*(a*x + 1)**(sympy.S(9)/2)*asin(a*x)/(a**10*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)) + 2*(-a*x + 1)**5*(a*x + 1)**4*(1019*a*x + 704)/(315*a**10*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_701():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(3*atanh(a*x))
    F = -a**9*x**10*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 3*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 4*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 2*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(9)/2)) + 8*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - 3*a*x**2*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(7*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(8*(-a**2*x**2 + 1)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_702():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(3*atanh(a*x))
    F = a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + 3*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(7)/2) - a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + 5*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + 5*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(3*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - 3*a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(6*(-a**2*x**2 + 1)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_703():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(3*atanh(a*x))
    F = -a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - x*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_704():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(3*atanh(a*x))
    F = a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) + 3*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - 3*a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - x*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_705():
    f = sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))
    F = -a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_706():
    f = exp(3*atanh(a*x))/sqrt(c - c/(a**2*x**2))
    F = sqrt(-a**2*x**2 + 1)/(a*sqrt(c - c/(a**2*x**2))) + 3*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))) + 2*sqrt(-a**2*x**2 + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_707():
    f = exp(3*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = -(-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)*log(-a*x + 1)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_708():
    f = exp(3*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = (-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 49*(-a**2*x**2 + 1)**(sympy.S(5)/2)*log(-a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - (-a**2*x**2 + 1)**(sympy.S(5)/2)*log(a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 31*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1)) - 9*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1)**2) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(6*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_709():
    f = exp(3*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = -(-a**2*x**2 + 1)**(sympy.S(7)/2)/(a**7*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 201*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(-a*x + 1)/(64*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 9*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(a*x + 1)/(64*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)) - 75*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(16*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)) + 59*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)**2) - (-a**2*x**2 + 1)**(sympy.S(7)/2)/(2*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(16*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_710():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(-atanh(a*x))
    F = -a**9*x**10*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 4*a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + 2*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + 2*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - 4*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(9)/2)) + 2*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(3*(-a**2*x**2 + 1)**(sympy.S(9)/2)) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(7*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(8*(-a**2*x**2 + 1)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_711():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(-atanh(a*x))
    F = a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) - a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + 3*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) - 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + 3*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(6*(-a**2*x**2 + 1)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_712():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(-atanh(a*x))
    F = -a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_713():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(-atanh(a*x))
    F = a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(3)/2) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - x*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_714():
    f = sqrt(c - c/(a**2*x**2))*exp(-atanh(a*x))
    F = -a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_715():
    f = exp(-atanh(a*x))/sqrt(c - c/(a**2*x**2))
    F = sqrt(-a**2*x**2 + 1)/(a*sqrt(c - c/(a**2*x**2))) - sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_716():
    f = exp(-atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = -(-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - (-a**2*x**2 + 1)**(sympy.S(3)/2)*log(-a*x + 1)/(4*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + 5*(-a**2*x**2 + 1)**(sympy.S(3)/2)*log(a*x + 1)/(4*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_717():
    f = exp(-atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = (-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 7*(-a**2*x**2 + 1)**(sympy.S(5)/2)*log(-a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 23*(-a**2*x**2 + 1)**(sympy.S(5)/2)*log(a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - (-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)**2) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_718():
    f = exp(-atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = -(-a**2*x**2 + 1)**(sympy.S(7)/2)/(a**7*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 19*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(-a*x + 1)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 51*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(a*x + 1)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 3*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(2*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)) - 11*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)**2) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(24*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)**3) - 5*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(16*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_719():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(-2*atanh(a*x))
    F = 11*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(128*(-a*x + 1)**4*(a*x + 1)**4) - 2*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(9)/2)*(a*x + 1)**(sympy.S(9)/2)) + 245*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(128*(-a*x + 1)**(sympy.S(9)/2)*(a*x + 1)**(sympy.S(9)/2)) + 39*a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(64*(-a*x + 1)**4*(a*x + 1)**3) - 11*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(640*(-a*x + 1)**4*(a*x + 1)**2) - 103*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(160*(-a*x + 1)**4*(a*x + 1)) + 629*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(960*(-a*x + 1)**3*(a*x + 1)) - 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(5*(-a*x + 1)**2*(a*x + 1)) + 47*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(9)/2)/((-336*a*x + 336)*(a*x + 1)) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(28*a*x + 28) - x*(c - c/(a**2*x**2))**(sympy.S(9)/2)*(-a*x + 1)/(8*a*x + 8)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_720():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(-2*atanh(a*x))
    F = -7*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(16*(-a*x + 1)**3*(a*x + 1)**3) + 2*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) - 25*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(16*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)) - 3*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(8*(-a*x + 1)**3*(a*x + 1)**2) + 19*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(16*(-a*x + 1)**3*(a*x + 1)) - 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(3*(-a*x + 1)**2*(a*x + 1)) + 23*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)/((-120*a*x + 120)*(a*x + 1)) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(15*a*x + 15) - x*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1)/(6*a*x + 6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_721():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(-2*atanh(a*x))
    F = 7*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(8*(-a*x + 1)**2*(a*x + 1)**2) - 2*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) + 9*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)) - 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/((-a*x + 1)**2*(a*x + 1)) + 7*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/((-24*a*x + 24)*(a*x + 1)) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(6*a*x + 6) - x*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(-a*x + 1)/(4*a*x + 4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_722():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(-2*atanh(a*x))
    F = 2*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*asin(a*x)/((-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)) + 5*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)/((-2*a*x + 2)*(a*x + 1)) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(a*x + 1) - x*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)/(2*a*x + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_723():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))
    F = -x*sqrt(c - c/(a**2*x**2)) - 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_724():
    f = exp(-2*atanh(a*x))/sqrt(c - c/(a**2*x**2))
    F = (-2*a*x + 2)*(a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))) + 2*sqrt(-a*x + 1)*sqrt(a*x + 1)*asin(a*x)/(a**2*x*sqrt(c - c/(a**2*x**2))) + (-a*x + 1)**2/(a**2*x*sqrt(c - c/(a**2*x**2)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_725():
    f = exp(-2*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = (-a*x + 1)**2/(3*a**2*x*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - 2*(-a*x + 1)**(sympy.S(3)/2)*(a*x + 1)**(sympy.S(3)/2)*asin(a*x)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - 2*(-a*x + 1)**2*(a*x + 1)*(2*a*x + 5)/(3*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_726():
    f = exp(-2*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = (-a*x + 1)**2/(a**2*x*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 2*(-a*x + 1)**3/(5*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - 2*(-a*x + 1)**3*(a*x + 1)/(15*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 2*(-a*x + 1)**(sympy.S(5)/2)*(a*x + 1)**(sympy.S(5)/2)*asin(a*x)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 2*(-a*x + 1)**3*(a*x + 1)**2*(13*a*x + 28)/(15*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_727():
    f = exp(-2*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = (-a*x + 1)**2/(3*a**2*x*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 10*(-a*x + 1)**3/(3*a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 12*(-a*x + 1)**4/(7*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 82*(-a*x + 1)**4*(a*x + 1)/(105*a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 2*(-a*x + 1)**4*(a*x + 1)**2/(35*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 2*(-a*x + 1)**(sympy.S(7)/2)*(a*x + 1)**(sympy.S(7)/2)*asin(a*x)/(a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 2*(-a*x + 1)**4*(a*x + 1)**3*(37*a*x + 72)/(35*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_728():
    f = (c - c/(a**2*x**2))**(sympy.S(9)/2)*exp(-3*atanh(a*x))
    F = a**9*x**10*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 3*a**8*x**9*(c - c/(a**2*x**2))**(sympy.S(9)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(9)/2) - 4*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + 2*a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2) + 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - 8*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(9)/2)) + 3*a*x**2*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(7*(-a**2*x**2 + 1)**(sympy.S(9)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(9)/2)/(8*(-a**2*x**2 + 1)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_729():
    f = (c - c/(a**2*x**2))**(sympy.S(7)/2)*exp(-3*atanh(a*x))
    F = -a**7*x**8*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + 3*a**6*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2) + 5*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - 5*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(3*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + 3*a*x**2*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(5*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(7)/2)/(6*(-a**2*x**2 + 1)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_730():
    f = (c - c/(a**2*x**2))**(sympy.S(5)/2)*exp(-3*atanh(a*x))
    F = a**5*x**6*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - 3*a**4*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - 2*a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) + a*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2) - x*(c - c/(a**2*x**2))**(sympy.S(5)/2)/(4*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_731():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(-3*atanh(a*x))
    F = -a**3*x**4*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) + 3*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*log(x)/(-a**2*x**2 + 1)**(sympy.S(3)/2) + 3*a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2) - x*(c - c/(a**2*x**2))**(sympy.S(3)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_732():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))
    F = a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_733():
    f = exp(-3*atanh(a*x))/sqrt(c - c/(a**2*x**2))
    F = -sqrt(-a**2*x**2 + 1)/(a*sqrt(c - c/(a**2*x**2))) + 3*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))) + 2*sqrt(-a**2*x**2 + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_734():
    f = exp(-3*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = (-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)*log(a*x + 1)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(a*x + 1)) + (-a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**4*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_735():
    f = exp(-3*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = -(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - (-a**2*x**2 + 1)**(sympy.S(5)/2)*log(-a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 49*(-a**2*x**2 + 1)**(sympy.S(5)/2)*log(a*x + 1)/(16*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 31*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)) - 9*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(8*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)**2) + (-a**2*x**2 + 1)**(sympy.S(5)/2)/(6*a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_736():
    f = exp(-3*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(7)/2)
    F = (-a**2*x**2 + 1)**(sympy.S(7)/2)/(a**7*x**6*(c - c/(a**2*x**2))**(sympy.S(7)/2)) + 9*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(-a*x + 1)/(64*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 201*(-a**2*x**2 + 1)**(sympy.S(7)/2)*log(a*x + 1)/(64*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)) - 75*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(16*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)) + 59*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)**2) - (-a**2*x**2 + 1)**(sympy.S(7)/2)/(2*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)**3) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(16*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(a*x + 1)**4) + (-a**2*x**2 + 1)**(sympy.S(7)/2)/(32*a**8*x**7*(c - c/(a**2*x**2))**(sympy.S(7)/2)*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_737():
    f = x**m*sqrt(c - c/(a**2*x**2))*exp(atanh(a*x))
    F = a*x**(m + 2)*sqrt(c - c/(a**2*x**2))/((m + 1)*sqrt(-a**2*x**2 + 1)) + x**(m + 1)*sqrt(c - c/(a**2*x**2))/(m*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_738():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(atanh(a*x))
    F = a*x**4*sqrt(c - c/(a**2*x**2))/(3*sqrt(-a**2*x**2 + 1)) + x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_739():
    f = x*sqrt(c - c/(a**2*x**2))*exp(atanh(a*x))
    F = a*x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1)) + x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_740():
    f = sqrt(c - c/(a**2*x**2))*exp(atanh(a*x))
    F = a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_741():
    f = sqrt(c - c/(a**2*x**2))*exp(atanh(a*x))/x
    F = a*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_742():
    f = sqrt(c - c/(a**2*x**2))*exp(atanh(a*x))/x**2
    F = -sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(2*x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_743():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))
    F = -x**2*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(4*a**2) - x*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(6*a**3) - 7*x*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(24*a**3) - 7*x*sqrt(c - c/(a**2*x**2))/(8*a**3) + 7*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(8*a**3*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_744():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))
    F = -x*sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(3*a**2) - x*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(3*a**2) - x*sqrt(c - c/(a**2*x**2))/a**2 + x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(a**2*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_745():
    f = x*sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))
    F = -x*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(2*a) - 3*x*sqrt(c - c/(a**2*x**2))/(2*a) + 3*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(2*a*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_746():
    f = sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))
    F = -x*sqrt(c - c/(a**2*x**2)) + 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_747():
    f = sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))/x
    F = a*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - 2*a*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - sqrt(c - c/(a**2*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_748():
    f = sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))/x**2
    F = -3*a**2*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 3*a*sqrt(c - c/(a**2*x**2))/2 - sqrt(c - c/(a**2*x**2))*(a*x + 1)/(2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_749():
    f = sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))/x**3
    F = -a**3*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - a**2*sqrt(c - c/(a**2*x**2)) - a*sqrt(c - c/(a**2*x**2))*(a*x + 1)/(3*x) - sqrt(c - c/(a**2*x**2))*(a*x + 1)**2/(3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_750():
    f = sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))/x**4
    F = -7*a**4*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 4*a**3*sqrt(c - c/(a**2*x**2))/3 - 7*a**2*sqrt(c - c/(a**2*x**2))/(8*x) - 2*a*sqrt(c - c/(a**2*x**2))/(3*x**2) - sqrt(c - c/(a**2*x**2))/(4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_751():
    f = sqrt(c - c/(a**2*x**2))*exp(2*atanh(a*x))/x**5
    F = -3*a**5*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(4*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 6*a**4*sqrt(c - c/(a**2*x**2))/5 - 3*a**3*sqrt(c - c/(a**2*x**2))/(4*x) - 3*a**2*sqrt(c - c/(a**2*x**2))/(5*x**2) - a*sqrt(c - c/(a**2*x**2))/(2*x**3) - sqrt(c - c/(a**2*x**2))/(5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_752():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))
    F = -a*x**5*sqrt(c - c/(a**2*x**2))/(4*sqrt(-a**2*x**2 + 1)) - x**4*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - 2*x**3*sqrt(c - c/(a**2*x**2))/(a*sqrt(-a**2*x**2 + 1)) - 4*x**2*sqrt(c - c/(a**2*x**2))/(a**2*sqrt(-a**2*x**2 + 1)) - 4*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_753():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))
    F = -a*x**4*sqrt(c - c/(a**2*x**2))/(3*sqrt(-a**2*x**2 + 1)) - 3*x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1)) - 4*x**2*sqrt(c - c/(a**2*x**2))/(a*sqrt(-a**2*x**2 + 1)) - 4*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_754():
    f = x*sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))
    F = -a*x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1)) - 3*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - 4*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/(a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_755():
    f = sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))
    F = -a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_756():
    f = sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))/x
    F = 3*a*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*a*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_757():
    f = sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))/x**2
    F = 4*a**2*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**2*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - 3*a*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - sqrt(c - c/(a**2*x**2))/(2*x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_758():
    f = sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))/x**3
    F = 4*a**3*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**3*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - 3*a*sqrt(c - c/(a**2*x**2))/(2*x*sqrt(-a**2*x**2 + 1)) - sqrt(c - c/(a**2*x**2))/(3*x**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_759():
    f = sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))/x**4
    F = 4*a**4*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**4*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**3*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - 2*a**2*sqrt(c - c/(a**2*x**2))/(x*sqrt(-a**2*x**2 + 1)) - a*sqrt(c - c/(a**2*x**2))/(x**2*sqrt(-a**2*x**2 + 1)) - sqrt(c - c/(a**2*x**2))/(4*x**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_760():
    f = sqrt(c - c/(a**2*x**2))*exp(3*atanh(a*x))/x**5
    F = 4*a**5*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**5*x*sqrt(c - c/(a**2*x**2))*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**4*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - 2*a**3*sqrt(c - c/(a**2*x**2))/(x*sqrt(-a**2*x**2 + 1)) - 4*a**2*sqrt(c - c/(a**2*x**2))/(3*x**2*sqrt(-a**2*x**2 + 1)) - 3*a*sqrt(c - c/(a**2*x**2))/(4*x**3*sqrt(-a**2*x**2 + 1)) - sqrt(c - c/(a**2*x**2))/(5*x**4*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_761():
    f = x**m*sqrt(c - c/(a**2*x**2))*exp(-atanh(a*x))
    F = -a*x**(m + 2)*sqrt(c - c/(a**2*x**2))/((m + 1)*sqrt(-a**2*x**2 + 1)) + x**(m + 1)*sqrt(c - c/(a**2*x**2))/(m*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_762():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(-atanh(a*x))
    F = -a*x**4*sqrt(c - c/(a**2*x**2))/(3*sqrt(-a**2*x**2 + 1)) + x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_763():
    f = x*sqrt(c - c/(a**2*x**2))*exp(-atanh(a*x))
    F = -a*x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1)) + x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_764():
    f = sqrt(c - c/(a**2*x**2))*exp(-atanh(a*x))
    F = -a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_765():
    f = sqrt(c - c/(a**2*x**2))*exp(-atanh(a*x))/x
    F = -a*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_766():
    f = sqrt(c - c/(a**2*x**2))*exp(-atanh(a*x))/x**2
    F = -sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(2*x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_767():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))
    F = -x**2*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(4*a**2) + x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(6*a**3) + 7*x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(24*a**3) + 7*x*sqrt(c - c/(a**2*x**2))/(8*a**3) + 7*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(8*a**3*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_768():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))
    F = -x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(3*a**2) - x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(3*a**2) - x*sqrt(c - c/(a**2*x**2))/a**2 - x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(a**2*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_769():
    f = x*sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))
    F = x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(2*a) + 3*x*sqrt(c - c/(a**2*x**2))/(2*a) + 3*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(2*a*sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_770():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))
    F = -x*sqrt(c - c/(a**2*x**2)) - 2*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_771():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))/x
    F = a*x*sqrt(c - c/(a**2*x**2))*asin(a*x)/(sqrt(-a*x + 1)*sqrt(a*x + 1)) + 2*a*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - sqrt(c - c/(a**2*x**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_772():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))/x**2
    F = -3*a**2*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(2*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 3*a*sqrt(c - c/(a**2*x**2))/2 - sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_773():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))/x**3
    F = a**3*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(sqrt(-a*x + 1)*sqrt(a*x + 1)) - a**2*sqrt(c - c/(a**2*x**2)) + a*sqrt(c - c/(a**2*x**2))*(-a*x + 1)/(3*x) - sqrt(c - c/(a**2*x**2))*(-a*x + 1)**2/(3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_774():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))/x**4
    F = -7*a**4*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(8*sqrt(-a*x + 1)*sqrt(a*x + 1)) + 4*a**3*sqrt(c - c/(a**2*x**2))/3 - 7*a**2*sqrt(c - c/(a**2*x**2))/(8*x) + 2*a*sqrt(c - c/(a**2*x**2))/(3*x**2) - sqrt(c - c/(a**2*x**2))/(4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_775():
    f = sqrt(c - c/(a**2*x**2))*exp(-2*atanh(a*x))/x**5
    F = 3*a**5*x*sqrt(c - c/(a**2*x**2))*atanh(sqrt(-a*x + 1)*sqrt(a*x + 1))/(4*sqrt(-a*x + 1)*sqrt(a*x + 1)) - 6*a**4*sqrt(c - c/(a**2*x**2))/5 + 3*a**3*sqrt(c - c/(a**2*x**2))/(4*x) - 3*a**2*sqrt(c - c/(a**2*x**2))/(5*x**2) + a*sqrt(c - c/(a**2*x**2))/(2*x**3) - sqrt(c - c/(a**2*x**2))/(5*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_776():
    f = x**3*sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))
    F = a*x**5*sqrt(c - c/(a**2*x**2))/(4*sqrt(-a**2*x**2 + 1)) - x**4*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + 2*x**3*sqrt(c - c/(a**2*x**2))/(a*sqrt(-a**2*x**2 + 1)) - 4*x**2*sqrt(c - c/(a**2*x**2))/(a**2*sqrt(-a**2*x**2 + 1)) + 4*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_777():
    f = x**2*sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))
    F = a*x**4*sqrt(c - c/(a**2*x**2))/(3*sqrt(-a**2*x**2 + 1)) - 3*x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1)) + 4*x**2*sqrt(c - c/(a**2*x**2))/(a*sqrt(-a**2*x**2 + 1)) - 4*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_778():
    f = x*sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))
    F = a*x**3*sqrt(c - c/(a**2*x**2))/(2*sqrt(-a**2*x**2 + 1)) - 3*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + 4*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/(a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_779():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))
    F = a*x**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_780():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))/x
    F = -3*a*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) + 4*a*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*x**2 + 1) - sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_781():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))/x**2
    F = 4*a**2*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**2*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*x**2 + 1) + 3*a*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - sqrt(c - c/(a**2*x**2))/(2*x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_782():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))/x**3
    F = -4*a**3*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) + 4*a**3*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**2*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + 3*a*sqrt(c - c/(a**2*x**2))/(2*x*sqrt(-a**2*x**2 + 1)) - sqrt(c - c/(a**2*x**2))/(3*x**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_783():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))/x**4
    F = 4*a**4*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**4*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*x**2 + 1) + 4*a**3*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) - 2*a**2*sqrt(c - c/(a**2*x**2))/(x*sqrt(-a**2*x**2 + 1)) + a*sqrt(c - c/(a**2*x**2))/(x**2*sqrt(-a**2*x**2 + 1)) - sqrt(c - c/(a**2*x**2))/(4*x**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_784():
    f = sqrt(c - c/(a**2*x**2))*exp(-3*atanh(a*x))/x**5
    F = -4*a**5*x*sqrt(c - c/(a**2*x**2))*log(x)/sqrt(-a**2*x**2 + 1) + 4*a**5*x*sqrt(c - c/(a**2*x**2))*log(a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**4*sqrt(c - c/(a**2*x**2))/sqrt(-a**2*x**2 + 1) + 2*a**3*sqrt(c - c/(a**2*x**2))/(x*sqrt(-a**2*x**2 + 1)) - 4*a**2*sqrt(c - c/(a**2*x**2))/(3*x**2*sqrt(-a**2*x**2 + 1)) + 3*a*sqrt(c - c/(a**2*x**2))/(4*x**3*sqrt(-a**2*x**2 + 1)) - sqrt(c - c/(a**2*x**2))/(5*x**4*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_785():
    f = (c - c/(a**2*x**2))**p*exp(-2*p*atanh(a*x))
    F = x*(c - c/(a**2*x**2))**p*hyper((-2*p, 1 - 2*p), (2 - 2*p,), a*x)/((1 - 2*p)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_786():
    f = (c - c/(a**2*x**2))**p*exp(2*p*atanh(a*x))
    F = x*(c - c/(a**2*x**2))**p*hyper((-2*p, 1 - 2*p), (2 - 2*p,), -a*x)/((1 - 2*p)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_787():
    f = (c - c/(a**2*x**2))**2*exp(n*atanh(a*x))
    F = 2**(n/2 - 1)*c**2*n*(-a*x + 1)**(3 - n/2)*hyper((2 - n/2, 3 - n/2), (4 - n/2,), -a*x/2 + sympy.S.Half)/(a*(n**2 - 10*n + 24)) - c**2*n*(10 - n**2)*(-a*x + 1)**(2 - n/2)*(a*x + 1)**(n/2 - 2)*hyper((1, n/2 - 2), (n/2 - 1,), (a*x + 1)/(-a*x + 1))/(3*a*(4 - n)) - 4*c**2*(-a*x + 1)**(3 - n/2)*(a*x + 1)**(n/2 - 2)/(a*(4 - n)) - c**2*(-a*x + 1)**(3 - n/2)*(a*x + 1)**(n/2 - 2)*(n**2 + 5*n + 14)/(6*a**2*x) - c**2*(n + 10)*(-a*x + 1)**(3 - n/2)*(a*x + 1)**(n/2 - 2)/(6*a**3*x**2) - c**2*(-a*x + 1)**(3 - n/2)*(a*x + 1)**(n/2 - 2)/(3*a**4*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_788():
    f = (c - c/(a**2*x**2))*exp(n*atanh(a*x))
    F = -2**(n/2 + 1)*c*(-a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(a*(2 - n)) + 4*c*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (-a*x + 1)/(a*x + 1))/(a*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_789():
    f = exp(n*atanh(a*x))/(c - c/(a**2*x**2))
    F = -2**(n/2 + 1)*hyper((-n/2, -n/2), (1 - n/2,), -a*x/2 + sympy.S.Half)/(a*c*(-a*x + 1)**(n/2)) + x*(a*x + 1)**(n/2)/(c*(-a*x + 1)**(n/2)) - (1 - n)*(a*x + 1)**(n/2)/(a*c*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_790():
    f = exp(n*atanh(a*x))/(c - c/(a**2*x**2))**2
    F = -2**(n/2)*n*(-a*x + 1)**(1 - n/2)*hyper((1 - n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(a*c**2*(2 - n)) - a**2*x**3*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/c**2 + x*(n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/c**2 + (1 - n)*(n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(a*c**2*(2 - n)) - (2 - n**2)*(n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2)/(a*c**2*(4 - n**2)) - (a*x + 1)**(n/2 - 1)/(a*c**2*(-a*x + 1)**(n/2)) + (-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)/(a*c**2*(2 - n)) - (2 - n**2)*(n + 3)*(a*x + 1)**(n/2)/(a*c**2*n*(4 - n**2)*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_791():
    f = (c - c/(a**2*x**2))**(sympy.S(3)/2)*exp(n*atanh(a*x))
    F = 2**(n/2 + sympy.S(-1)/2)*a**2*n*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)**(sympy.S(5)/2 - n/2)*hyper((sympy.S(3)/2 - n/2, sympy.S(5)/2 - n/2), (sympy.S(7)/2 - n/2,), -a*x/2 + sympy.S.Half)/((3 - n)*(5 - n)*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - a**2*x**3*(3 - n**2)*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*hyper((1, n/2 + sympy.S(-3)/2), (n/2 + sympy.S(-1)/2,), (a*x + 1)/(-a*x + 1))/((3 - n)*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 3*a**2*x**3*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)**(sympy.S(5)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)/((3 - n)*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - a*x**2*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(n + 4)*(-a*x + 1)**(sympy.S(5)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(c - c/(a**2*x**2))**(sympy.S(3)/2)*(-a*x + 1)**(sympy.S(5)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)/(2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_792():
    f = sqrt(c - c/(a**2*x**2))*exp(n*atanh(a*x))
    F = 2**(n/2 + sympy.S.Half)*n*x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**(sympy.S(3)/2 - n/2)*hyper((sympy.S.Half - n/2, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(sqrt(-a**2*x**2 + 1)*(n**2 - 4*n + 3)) + 2*x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/((1 - n)*sqrt(-a**2*x**2 + 1)) - x*sqrt(c - c/(a**2*x**2))*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)/((1 - n)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_793():
    f = exp(n*atanh(a*x))/sqrt(c - c/(a**2*x**2))
    F = -2**(n/2 + sympy.S(3)/2)*n*(-a*x + 1)**(sympy.S.Half - n/2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(3)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a**2*x*(1 - n**2)*sqrt(c - c/(a**2*x**2))) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(a**2*x*sqrt(c - c/(a**2*x**2))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_794():
    f = exp(n*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(3)/2)
    F = -2**(n/2 + sympy.S(-1)/2)*n*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(-a**2*x**2 + 1)**(sympy.S(3)/2)*hyper((sympy.S(3)/2 - n/2, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a**4*x**3*(3 - n)*(c - c/(a**2*x**2))**(sympy.S(3)/2)) - (-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(a**2*x*(c - c/(a**2*x**2))**(sympy.S(3)/2)) + (-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*(-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a*n*x*(2*n + 3) + n**2 + 2*n + 2)/(a**4*x**3*(1 - n**2)*(c - c/(a**2*x**2))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_795():
    f = exp(n*atanh(a*x))/(c - c/(a**2*x**2))**(sympy.S(5)/2)
    F = -2**(n/2 + sympy.S(3)/2)*n*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)*hyper((-n/2 + sympy.S(-1)/2, -n/2 + sympy.S(-1)/2), (sympy.S.Half - n/2,), -a*x/2 + sympy.S.Half)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 1)) - (-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**2*x*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + (n + 4)*(-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**3*x**2*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 3)) - (3*n + 12)*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**5*x**4*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 3)) + 2*n*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 1)*(n**2 - 4*n + 3)) + n*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 1)) - 3*n*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 1)) + 3*n*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S.Half)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 1)) - 2*n*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(1 - n**2)*(c - c/(a**2*x**2))**(sympy.S(5)/2)) + 3*n*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(1 - n**2)*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - (6 - 3*n)*(n + 4)*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(a**6*x**5*(9 - n**2)*(c - c/(a**2*x**2))**(sympy.S(5)/2)) - (3*n + 12)*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)*(-n**2 + 2*n + 1)/(a**6*x**5*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n**4 - 10*n**2 + 9)) + (3*n + 12)*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*(-a**2*x**2 + 1)**(sympy.S(5)/2)*(-n**2 + 2*n + 1)/(a**6*x**5*(3 - n)*(c - c/(a**2*x**2))**(sympy.S(5)/2)*(n + 1)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_796():
    f = (c - c/(a**2*x**2))**p*exp(n*atanh(a*x))
    F = x*(c - c/(a**2*x**2))**p*appellf1(1 - 2*p, -n/2 - p, n/2 - p, 2 - 2*p, -a*x, a*x)/((1 - 2*p)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_797():
    f = (c - c/(a**2*x**2))**p*exp(4*atanh(a*x))
    F = a**4*x**5*(c - c/(a**2*x**2))**p*hyper((2 - p, sympy.S(5)/2 - p), (sympy.S(7)/2 - p,), a**2*x**2)/((5 - 2*p)*(-a*x + 1)**p*(a*x + 1)**p) + 2*a**3*x**4*(c - c/(a**2*x**2))**p*hyper((2 - p, 2 - p), (3 - p,), a**2*x**2)/((2 - p)*(-a*x + 1)**p*(a*x + 1)**p) + 6*a**2*x**3*(c - c/(a**2*x**2))**p*hyper((sympy.S(3)/2 - p, 2 - p), (sympy.S(5)/2 - p,), a**2*x**2)/((3 - 2*p)*(-a*x + 1)**p*(a*x + 1)**p) + 2*a*x**2*(c - c/(a**2*x**2))**p/((1 - p)*(-a*x + 1)*(a*x + 1)) + x*(c - c/(a**2*x**2))**p*hyper((sympy.S.Half - p, 2 - p), (sympy.S(3)/2 - p,), a**2*x**2)/((1 - 2*p)*(-a*x + 1)**p*(a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_798():
    f = (c - c/(a**2*x**2))**p*exp(3*atanh(a*x))
    F = 3*a**2*x**3*(c - c/(a**2*x**2))**p*hyper((sympy.S(3)/2 - p, sympy.S(3)/2 - p), (sympy.S(5)/2 - p,), a**2*x**2)/((3 - 2*p)*(-a**2*x**2 + 1)**p) - a*x**2*(c - c/(a**2*x**2))**p/sqrt(-a**2*x**2 + 1) + a*x**2*(5 - 2*p)*(c - c/(a**2*x**2))**p*hyper((1 - p, sympy.S(3)/2 - p), (2 - p,), a**2*x**2)/((2 - 2*p)*(-a**2*x**2 + 1)**p) + x*(c - c/(a**2*x**2))**p/((1 - 2*p)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_799():
    f = (c - c/(a**2*x**2))**p*exp(2*atanh(a*x))
    F = a**2*x**3*(c - c/(a**2*x**2))**p*hyper((1 - p, sympy.S(3)/2 - p), (sympy.S(5)/2 - p,), a**2*x**2)/((3 - 2*p)*(-a*x + 1)**p*(a*x + 1)**p) + a*x**2*(c - c/(a**2*x**2))**p*hyper((1 - p, 1 - p), (2 - p,), a**2*x**2)/((1 - p)*(-a*x + 1)**p*(a*x + 1)**p) + x*(c - c/(a**2*x**2))**p*hyper((sympy.S.Half - p, 1 - p), (sympy.S(3)/2 - p,), a**2*x**2)/((1 - 2*p)*(-a*x + 1)**p*(a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_800():
    f = (c - c/(a**2*x**2))**p*exp(atanh(a*x))
    F = a*x**2*(c - c/(a**2*x**2))**p*hyper((sympy.S.Half - p, 1 - p), (2 - p,), a**2*x**2)/((2 - 2*p)*(-a**2*x**2 + 1)**p) + x*(c - c/(a**2*x**2))**p*hyper((sympy.S.Half - p, sympy.S.Half - p), (sympy.S(3)/2 - p,), a**2*x**2)/((1 - 2*p)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_801():
    f = (c - c/(a**2*x**2))**p*exp(-atanh(a*x))
    F = -a*x**2*(c - c/(a**2*x**2))**p*hyper((sympy.S.Half - p, 1 - p), (2 - p,), a**2*x**2)/((2 - 2*p)*(-a**2*x**2 + 1)**p) + x*(c - c/(a**2*x**2))**p*hyper((sympy.S.Half - p, sympy.S.Half - p), (sympy.S(3)/2 - p,), a**2*x**2)/((1 - 2*p)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_802():
    f = (c - c/(a**2*x**2))**p*exp(-2*atanh(a*x))
    F = a**2*x**3*(c - c/(a**2*x**2))**p*hyper((1 - p, sympy.S(3)/2 - p), (sympy.S(5)/2 - p,), a**2*x**2)/((3 - 2*p)*(-a*x + 1)**p*(a*x + 1)**p) - a*x**2*(c - c/(a**2*x**2))**p*hyper((1 - p, 1 - p), (2 - p,), a**2*x**2)/((1 - p)*(-a*x + 1)**p*(a*x + 1)**p) + x*(c - c/(a**2*x**2))**p*hyper((sympy.S.Half - p, 1 - p), (sympy.S(3)/2 - p,), a**2*x**2)/((1 - 2*p)*(-a*x + 1)**p*(a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_803():
    f = (c - c/(a**2*x**2))**p*exp(-3*atanh(a*x))
    F = 3*a**2*x**3*(c - c/(a**2*x**2))**p*hyper((sympy.S(3)/2 - p, sympy.S(3)/2 - p), (sympy.S(5)/2 - p,), a**2*x**2)/((3 - 2*p)*(-a**2*x**2 + 1)**p) + a*x**2*(c - c/(a**2*x**2))**p/sqrt(-a**2*x**2 + 1) - a*x**2*(5 - 2*p)*(c - c/(a**2*x**2))**p*hyper((1 - p, sympy.S(3)/2 - p), (2 - p,), a**2*x**2)/((2 - 2*p)*(-a**2*x**2 + 1)**p) + x*(c - c/(a**2*x**2))**p/((1 - 2*p)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_804():
    f = x*sqrt(x + 1)*exp(atanh(x))*sin(x)
    F = (Integer(3) * sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.cos(x)) + (Integer(-1) * (((Integer(1) + (Integer(-1) * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos(x))) + (Integer(-1) * (Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1))) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + (Integer(-1) * (Integer(3) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.sin(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_805():
    f = sqrt(x + 1)*exp(atanh(x))*sin(x)
    F = (sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.cos(x)) + (Integer(-1) * (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x)))))) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + (Integer(-1) * (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_806():
    f = x*sqrt(1 - x)*exp(atanh(x))*sin(x)
    F = (sympy.sqrt((Integer(1) + x)) * sympy.cos(x)) + (Integer(-1) * (((Integer(1) + x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos(x))) + (Integer(-1) * (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1))) + (Integer(-1) * (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1)))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((Integer(1) + x)) * sympy.sin(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_807():
    f = sqrt(1 - x)*exp(atanh(x))*sin(x)
    F = ((Integer(-1) * sympy.sqrt((Integer(1) + x))) * sympy.cos(x)) + (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x))))) + (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_808():
    f = x*(x + 1)**(sympy.S(3)/2)*exp(atanh(x))*sin(x)
    F = ((Integer(17) * (Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.cos(x)) + (Integer(-1) * (Integer(5) * ((Integer(1) + (Integer(-1) * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos(x))) + (((Integer(1) + (Integer(-1) * x)))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.cos(x)) + ((Integer(15) * (Integer(4))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x)))))) + (Integer(-1) * (Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (Integer(-1) * ((Integer(15) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x)))))) + ((Integer(15) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1))) + (Integer(-1) * (Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + ((Integer(15) * (Integer(4))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1))) + (Integer(-1) * (Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + (Integer(-1) * ((Integer(15) * (Integer(2))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.sin(x))) + ((Integer(5) * (Integer(2))**(Integer(-1))) * ((Integer(1) + (Integer(-1) * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin(x))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_809():
    f = (x + 1)**(sympy.S(3)/2)*exp(atanh(x))*sin(x)
    F = (Integer(4) * sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.cos(x)) + (Integer(-1) * (((Integer(1) + (Integer(-1) * x)))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos(x))) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x)))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1))) + (Integer(-1) * (Integer(4) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + (Integer(-1) * (Integer(2) * sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.sin(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_810():
    f = x*(1 - x)**(sympy.S(3)/2)*exp(atanh(x))*sin(x)
    F = ((Integer(-1) * (Integer(7) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)) * sympy.cos(x)) + (Integer(-1) * (Integer(3) * ((Integer(1) + x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos(x))) + (((Integer(1) + x))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.cos(x)) + ((Integer(7) * (Integer(4))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x))))) + (Integer(-1) * ((Integer(9) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))))) + ((Integer(9) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1))) + ((Integer(7) * (Integer(4))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1))) + ((Integer(9) * (Integer(2))**(Integer(-1))) * sympy.sqrt((Integer(1) + x)) * sympy.sin(x)) + (Integer(-1) * ((Integer(5) * (Integer(2))**(Integer(-1))) * ((Integer(1) + x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_811():
    f = (1 - x)**(sympy.S(3)/2)*exp(atanh(x))*sin(x)
    F = (Integer(-2) * sympy.sqrt((Integer(1) + x)) * sympy.cos(x)) + (((Integer(1) + x))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.cos(x)) + (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x))))) + ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x))))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1)))) + (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + x)))) * sympy.sin(Integer(1))) + (Integer(-1) * ((Integer(3) * (Integer(2))**(Integer(-1))) * sympy.sqrt((Integer(1) + x)) * sympy.sin(x)))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_812():
    f = x*exp(atanh(x))*sin(x)/sqrt(x + 1)
    F = (sympy.sqrt((Integer(1) + (Integer(-1) * x))) * sympy.cos(x)) + (Integer(-1) * (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.cos(Integer(1)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))))) + (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x)))))) + (Integer(-1) * (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1)))) + (Integer(-1) * (sympy.sqrt((sympy.pi * (Integer(2))**(Integer(-1)))) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_813():
    f = exp(atanh(x))*sin(x)/sqrt(x + 1)
    F = (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.cos(Integer(1)) * sympy.Function('FresnelS')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x)))))) + (Integer(-1) * (sympy.sqrt((Integer(2) * sympy.pi)) * sympy.Function('FresnelC')((sympy.sqrt((Integer(2) * (sympy.pi)**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * x))))) * sympy.sin(Integer(1))))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_814():
    f = x**3*exp(atanh(a + b*x))
    F = -x**2*sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)/(4*b**2) - sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)*(18*a**2 - 10*a + b*x*(2 - 12*a) + 7)/(24*b**4) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(-8*a**3 + 12*a**2 - 12*a + 3)/(8*b**4) + (-8*a**3 + 12*a**2 - 12*a + 3)*asin(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_815():
    f = x**2*exp(atanh(a + b*x))
    F = -x*sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)/(3*b**2) - (1 - 4*a)*sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)/(6*b**3) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(2*a**2 - 2*a + 1)/(2*b**3) + (2*a**2 - 2*a + 1)*asin(a + b*x)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_816():
    f = x*exp(atanh(a + b*x))
    F = -(1 - 2*a)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(2*b**2) + (1 - 2*a)*asin(a + b*x)/(2*b**2) - sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_817():
    f = exp(atanh(a + b*x))
    F = -sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/b + asin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_818():
    f = exp(atanh(a + b*x))/x
    F = asin(a + b*x) - (2*a + 2)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/sqrt(1 - a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_819():
    f = exp(atanh(a + b*x))/x**2
    F = -2*b*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)*sqrt(1 - a**2)) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(x*(1 - a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_820():
    f = exp(atanh(a + b*x))/x**3
    F = -b**2*(2*a + 1)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)**2*sqrt(1 - a**2)*(a + 1)) - b*(2*a + 1)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(2*x*(1 - a)**2*(a + 1)) - sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)/(x**2*(2 - 2*a**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_821():
    f = exp(atanh(a + b*x))/x**4
    F = -b**3*(2*a**2 + 2*a + 1)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)*(1 - a**2)**(sympy.S(5)/2)) - b**2*(a + 4)*(2*a + 1)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(6*x*(1 - a)**3*(a + 1)**2) - b*(2*a + 3)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(6*x**2*(1 - a)**2*(a + 1)) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(x**3*(3 - 3*a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_822():
    f = x**4*exp(2*atanh(a + b*x))
    F = -x**5/5 - x**4/(2*b) - x**3*(2 - 2*a)/(3*b**2) - x**2*(1 - a)**2/b**3 - 2*x*(1 - a)**3/b**4 - 2*(1 - a)**4*log(-a - b*x + 1)/b**5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_823():
    f = x**3*exp(2*atanh(a + b*x))
    F = -x**4/4 - 2*x**3/(3*b) - x**2*(1 - a)/b**2 - 2*x*(1 - a)**2/b**3 - 2*(1 - a)**3*log(-a - b*x + 1)/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_824():
    f = x**2*exp(2*atanh(a + b*x))
    F = -x**3/3 - x**2/b - x*(2 - 2*a)/b**2 - 2*(1 - a)**2*log(-a - b*x + 1)/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_825():
    f = x*exp(2*atanh(a + b*x))
    F = -x**2/2 - 2*x/b - (2 - 2*a)*log(-a - b*x + 1)/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_826():
    f = exp(2*atanh(a + b*x))
    F = -x - 2*log(-a - b*x + 1)/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_827():
    f = exp(2*atanh(a + b*x))/x
    F = (a + 1)*log(x)/(1 - a) - 2*log(-a - b*x + 1)/(1 - a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_828():
    f = exp(2*atanh(a + b*x))/x**2
    F = 2*b*log(x)/(1 - a)**2 - 2*b*log(-a - b*x + 1)/(1 - a)**2 - (a + 1)/(x*(1 - a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_829():
    f = exp(2*atanh(a + b*x))/x**3
    F = 2*b**2*log(x)/(1 - a)**3 - 2*b**2*log(-a - b*x + 1)/(1 - a)**3 - 2*b/(x*(1 - a)**2) - (a + 1)/(x**2*(2 - 2*a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_830():
    f = exp(2*atanh(a + b*x))/x**4
    F = 2*b**3*log(x)/(1 - a)**4 - 2*b**3*log(-a - b*x + 1)/(1 - a)**4 - 2*b**2/(x*(1 - a)**3) - b/(x**2*(1 - a)**2) - (a + 1)/(x**3*(3 - 3*a))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_831():
    f = x**3*exp(3*atanh(a + b*x))
    F = 2*x**3*(a + b*x + 1)**(sympy.S(3)/2)/(b*sqrt(-a - b*x + 1)) + 9*x**2*sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)/(4*b**2) + sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)*(22*a**2 - 54*a + b*x*(22 - 20*a) + 29)/(8*b**4) + sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(-24*a**3 + 108*a**2 - 132*a + 51)/(8*b**4) - (-24*a**3 + 108*a**2 - 132*a + 51)*asin(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_832():
    f = x**2*exp(3*atanh(a + b*x))
    F = (1 - a)**2*(a + b*x + 1)**(sympy.S(5)/2)/(b**3*sqrt(-a - b*x + 1)) + sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(5)/2)/(3*b**3) + sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)*(6*a**2 - 18*a + 11)/(6*b**3) + sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(6*a**2 - 18*a + 11)/(2*b**3) - (6*a**2 - 18*a + 11)*asin(a + b*x)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_833():
    f = x*exp(3*atanh(a + b*x))
    F = (1 - a)*(a + b*x + 1)**(sympy.S(5)/2)/(b**2*sqrt(-a - b*x + 1)) + (3 - 2*a)*sqrt(-a - b*x + 1)*(a + b*x + 1)**(sympy.S(3)/2)/(2*b**2) + (9 - 6*a)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(2*b**2) - (9 - 6*a)*asin(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_834():
    f = exp(3*atanh(a + b*x))
    F = 3*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/b - 3*asin(a + b*x)/b + 2*(a + b*x + 1)**(sympy.S(3)/2)/(b*sqrt(-a - b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_835():
    f = exp(3*atanh(a + b*x))/x
    F = -asin(a + b*x) + 4*sqrt(a + b*x + 1)/((1 - a)*sqrt(-a - b*x + 1)) - 2*(a + 1)**2*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)*sqrt(1 - a**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_836():
    f = exp(3*atanh(a + b*x))/x**2
    F = 6*b*sqrt(a + b*x + 1)/((1 - a)**2*sqrt(-a - b*x + 1)) - b*(6*a + 6)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)**2*sqrt(1 - a**2)) - (a + b*x + 1)**(sympy.S(3)/2)/(x*(1 - a)*sqrt(-a - b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_837():
    f = exp(3*atanh(a + b*x))/x**3
    F = b**2*(6*a + 9)*sqrt(a + b*x + 1)/((1 - a)**3*(a + 1)*sqrt(-a - b*x + 1)) - b**2*(6*a + 9)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)**3*sqrt(1 - a**2)) - b*(2*a + 3)*(a + b*x + 1)**(sympy.S(3)/2)/(2*x*(1 - a)**2*(a + 1)*sqrt(-a - b*x + 1)) - (a + b*x + 1)**(sympy.S(5)/2)/(x**2*(2 - 2*a**2)*sqrt(-a - b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_838():
    f = exp(3*atanh(a + b*x))/x**4
    F = b**3*sqrt(a + b*x + 1)*(2*a**2 + 51*a + 52)/(6*(1 - a)**4*(a + 1)*sqrt(-a - b*x + 1)) - b**3*(6*a**2 + 18*a + 11)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)**4*sqrt(1 - a**2)*(a + 1)) - b**2*(16*a + 19)*sqrt(a + b*x + 1)/(6*x*(1 - a)**3*(a + 1)*sqrt(-a - b*x + 1)) - 7*b*sqrt(a + b*x + 1)/(6*x**2*(1 - a)**2*sqrt(-a - b*x + 1)) - (a + 1)*sqrt(a + b*x + 1)/(x**3*(3 - 3*a)*sqrt(-a - b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_839():
    f = x**3*exp(-atanh(a + b*x))
    F = -x**2*(-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)/(4*b**2) - (-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)*(18*a**2 + 10*a - b*x*(12*a + 2) + 7)/(24*b**4) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(8*a**3 + 12*a**2 + 12*a + 3)/(8*b**4) - (8*a**3 + 12*a**2 + 12*a + 3)*asin(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_840():
    f = x**2*exp(-atanh(a + b*x))
    F = -x*(-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)/(3*b**2) + (4*a + 1)*(-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)/(6*b**3) + sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(2*a**2 + 2*a + 1)/(2*b**3) + (2*a**2 + 2*a + 1)*asin(a + b*x)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_841():
    f = x*exp(-atanh(a + b*x))
    F = -(2*a + 1)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(2*b**2) - (2*a + 1)*asin(a + b*x)/(2*b**2) - (-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_842():
    f = exp(-atanh(a + b*x))
    F = sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/b + asin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_843():
    f = exp(-atanh(a + b*x))/x
    F = -asin(a + b*x) - (2 - 2*a)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/sqrt(1 - a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_844():
    f = exp(-atanh(a + b*x))/x**2
    F = 2*b*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/(sqrt(1 - a**2)*(a + 1)) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(x*(a + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_845():
    f = exp(-atanh(a + b*x))/x**3
    F = -b**2*(1 - 2*a)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)*sqrt(1 - a**2)*(a + 1)**2) + b*(1 - 2*a)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(x*(2 - 2*a)*(a + 1)**2) - (-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)/(x**2*(2 - 2*a**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_846():
    f = exp(-atanh(a + b*x))/x**4
    F = b**3*(2*a**2 - 2*a + 1)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a**2)**(sympy.S(5)/2)*(a + 1)) - b**2*(1 - 2*a)*(4 - a)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(6*x*(1 - a)**2*(a + 1)**3) + b*(3 - 2*a)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(x**2*(6 - 6*a)*(a + 1)**2) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(x**3*(3*a + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_847():
    f = x**4*exp(-2*atanh(a + b*x))
    F = -x**5/5 + x**4/(2*b) - x**3*(2*a + 2)/(3*b**2) + x**2*(a + 1)**2/b**3 - 2*x*(a + 1)**3/b**4 + 2*(a + 1)**4*log(a + b*x + 1)/b**5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_848():
    f = x**3*exp(-2*atanh(a + b*x))
    F = -x**4/4 + 2*x**3/(3*b) - x**2*(a + 1)/b**2 + 2*x*(a + 1)**2/b**3 - 2*(a + 1)**3*log(a + b*x + 1)/b**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_849():
    f = x**2*exp(-2*atanh(a + b*x))
    F = -x**3/3 + x**2/b - x*(2*a + 2)/b**2 + 2*(a + 1)**2*log(a + b*x + 1)/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_850():
    f = x*exp(-2*atanh(a + b*x))
    F = -x**2/2 + 2*x/b - (2*a + 2)*log(a + b*x + 1)/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_851():
    f = exp(-2*atanh(a + b*x))
    F = -x + 2*log(a + b*x + 1)/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_852():
    f = exp(-2*atanh(a + b*x))/x
    F = (1 - a)*log(x)/(a + 1) - 2*log(a + b*x + 1)/(a + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_853():
    f = exp(-2*atanh(a + b*x))/x**2
    F = -2*b*log(x)/(a + 1)**2 + 2*b*log(a + b*x + 1)/(a + 1)**2 - (1 - a)/(x*(a + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_854():
    f = exp(-2*atanh(a + b*x))/x**3
    F = 2*b**2*log(x)/(a + 1)**3 - 2*b**2*log(a + b*x + 1)/(a + 1)**3 + 2*b/(x*(a + 1)**2) - (1 - a)/(x**2*(2*a + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_855():
    f = exp(-2*atanh(a + b*x))/x**4
    F = -2*b**3*log(x)/(a + 1)**4 + 2*b**3*log(a + b*x + 1)/(a + 1)**4 - 2*b**2/(x*(a + 1)**3) + b/(x**2*(a + 1)**2) - (1 - a)/(x**3*(3*a + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_856():
    f = x**3*exp(-3*atanh(a + b*x))
    F = -2*x**3*(-a - b*x + 1)**(sympy.S(3)/2)/(b*sqrt(a + b*x + 1)) + 9*x**2*(-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)/(4*b**2) + (-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)*(22*a**2 + 54*a - b*x*(20*a + 22) + 29)/(8*b**4) + sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(24*a**3 + 108*a**2 + 132*a + 51)/(8*b**4) + (24*a**3 + 108*a**2 + 132*a + 51)*asin(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_857():
    f = x**2*exp(-3*atanh(a + b*x))
    F = -(a + 1)**2*(-a - b*x + 1)**(sympy.S(5)/2)/(b**3*sqrt(a + b*x + 1)) - (-a - b*x + 1)**(sympy.S(5)/2)*sqrt(a + b*x + 1)/(3*b**3) - (-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)*(6*a**2 + 18*a + 11)/(6*b**3) - sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)*(6*a**2 + 18*a + 11)/(2*b**3) - (6*a**2 + 18*a + 11)*asin(a + b*x)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_858():
    f = x*exp(-3*atanh(a + b*x))
    F = (a + 1)*(-a - b*x + 1)**(sympy.S(5)/2)/(b**2*sqrt(a + b*x + 1)) + (2*a + 3)*(-a - b*x + 1)**(sympy.S(3)/2)*sqrt(a + b*x + 1)/(2*b**2) + (6*a + 9)*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(2*b**2) + (6*a + 9)*asin(a + b*x)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_859():
    f = exp(-3*atanh(a + b*x))
    F = -2*(-a - b*x + 1)**(sympy.S(3)/2)/(b*sqrt(a + b*x + 1)) - 3*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/b - 3*asin(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_860():
    f = exp(-3*atanh(a + b*x))/x
    F = -2*(1 - a)**2*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/(sqrt(1 - a**2)*(a + 1)) + asin(a + b*x) + 4*sqrt(-a - b*x + 1)/((a + 1)*sqrt(a + b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_861():
    f = exp(-3*atanh(a + b*x))/x**2
    F = -6*b*sqrt(-a - b*x + 1)/((a + 1)**2*sqrt(a + b*x + 1)) + b*(6 - 6*a)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/(sqrt(1 - a**2)*(a + 1)**2) - (-a - b*x + 1)**(sympy.S(3)/2)/(x*(a + 1)*sqrt(a + b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_862():
    f = exp(-3*atanh(a + b*x))/x**3
    F = -b**2*(9 - 6*a)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/(sqrt(1 - a**2)*(a + 1)**3) + b**2*(9 - 6*a)*sqrt(-a - b*x + 1)/((1 - a)*(a + 1)**3*sqrt(a + b*x + 1)) + b*(3 - 2*a)*(-a - b*x + 1)**(sympy.S(3)/2)/(x*(2 - 2*a)*(a + 1)**2*sqrt(a + b*x + 1)) - (-a - b*x + 1)**(sympy.S(5)/2)/(x**2*(2 - 2*a**2)*sqrt(a + b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_863():
    f = exp(-3*atanh(a + b*x))/x**4
    F = -b**3*sqrt(-a - b*x + 1)*(2*a**2 - 51*a + 52)/((6 - 6*a)*(a + 1)**4*sqrt(a + b*x + 1)) + b**3*(6*a**2 - 18*a + 11)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)*sqrt(1 - a**2)*(a + 1)**4) - b**2*(19 - 16*a)*sqrt(-a - b*x + 1)/(x*(6 - 6*a)*(a + 1)**3*sqrt(a + b*x + 1)) + 7*b*sqrt(-a - b*x + 1)/(6*x**2*(a + 1)**2*sqrt(a + b*x + 1)) - (1 - a)*sqrt(-a - b*x + 1)/(x**3*(3*a + 3)*sqrt(a + b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_864():
    f = exp(atanh(b*x + 1))/(b*x + 2)
    F = asin(b*x + 1)/b
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_865():
    f = x**3*exp(atanh(a + b*x))/(-a**2 - 2*a*b*x - b**2*x**2 + 1)
    F = x**2*(1 - a)*sqrt(a + b*x + 1)/(b**2*sqrt(-a - b*x + 1)) + (b*x*(3 - 2*a) + (1 - 2*a)*(4 - a))*sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/(2*b**4) - (6*a**2 - 6*a + 3)*asin(a + b*x)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_866():
    f = x**2*exp(atanh(a + b*x))/(-a**2 - 2*a*b*x - b**2*x**2 + 1)
    F = -(1 - 2*a)*asin(a + b*x)/b**3 + (1 - a)**2*sqrt(a + b*x + 1)/(b**3*sqrt(-a - b*x + 1)) + sqrt(-a - b*x + 1)*sqrt(a + b*x + 1)/b**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_867():
    f = x*exp(atanh(a + b*x))/(-a**2 - 2*a*b*x - b**2*x**2 + 1)
    F = (1 - a)*sqrt(a + b*x + 1)/(b**2*sqrt(-a - b*x + 1)) - asin(a + b*x)/b**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_868():
    f = exp(atanh(a + b*x))/(-a**2 - 2*a*b*x - b**2*x**2 + 1)
    F = sqrt(a + b*x + 1)/(b*sqrt(-a - b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_869():
    f = exp(atanh(a + b*x))/(x*(-a**2 - 2*a*b*x - b**2*x**2 + 1))
    F = sqrt(a + b*x + 1)/((1 - a)*sqrt(-a - b*x + 1)) - 2*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)*sqrt(1 - a**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_870():
    f = exp(atanh(a + b*x))/(x**2*(-a**2 - 2*a*b*x - b**2*x**2 + 1))
    F = b*(a + 2)*sqrt(a + b*x + 1)/((1 - a)**2*(a + 1)*sqrt(-a - b*x + 1)) - b*(4*a + 2)*atanh(sqrt(1 - a)*sqrt(a + b*x + 1)/(sqrt(a + 1)*sqrt(-a - b*x + 1)))/((1 - a)**2*sqrt(1 - a**2)*(a + 1)) - sqrt(a + b*x + 1)/(x*(1 - a**2)*sqrt(-a - b*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_871():
    f = x**m*exp(n*atanh(a + b*x))
    F = x**(m + 1)*(-b*x/(1 - a) + 1)**(n/2)*(a + b*x + 1)**(n/2)*appellf1(m + 1, -n/2, n/2, m + 2, -b*x/(a + 1), b*x/(1 - a))/((m + 1)*(b*x/(a + 1) + 1)**(n/2)*(-a - b*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_872():
    f = x**3*exp(n*atanh(a + b*x))
    F = 2**(n/2 - 2)*(-a - b*x + 1)**(1 - n/2)*(24*a**3 - 36*a**2*n + 12*a*(n**2 + 2) - n*(n**2 + 8))*hyper((-n/2, 1 - n/2), (2 - n/2,), -a/2 - b*x/2 + sympy.S.Half)/(3*b**4*(2 - n)) - x**2*(-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 + 1)/(4*b**2) - (-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 + 1)*(18*a**2 - 10*a*n - 2*b*x*(6*a - n) + n**2 + 6)/(24*b**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_873():
    f = x**2*exp(n*atanh(a + b*x))
    F = -2**(n/2)*(-a - b*x + 1)**(1 - n/2)*(6*a**2 - 6*a*n + n**2 + 2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a/2 - b*x/2 + sympy.S.Half)/(3*b**3*(2 - n)) - x*(-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 + 1)/(3*b**2) + (4*a - n)*(-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 + 1)/(6*b**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_874():
    f = x*exp(n*atanh(a + b*x))
    F = 2**(n/2)*(2*a - n)*(-a - b*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a/2 - b*x/2 + sympy.S.Half)/(b**2*(2 - n)) - (-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 + 1)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_875():
    f = exp(n*atanh(a + b*x))
    F = -2**(n/2 + 1)*(-a - b*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -a/2 - b*x/2 + sympy.S.Half)/(b*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_876():
    f = exp(n*atanh(a + b*x))/x
    F = -2**(n/2 + 1)*hyper((-n/2, -n/2), (1 - n/2,), -a/2 - b*x/2 + sympy.S.Half)/(n*(-a - b*x + 1)**(n/2)) + 2*(a + b*x + 1)**(n/2)*hyper((1, -n/2), (1 - n/2,), (a + 1)*(-a - b*x + 1)/((1 - a)*(a + b*x + 1)))/(n*(-a - b*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_877():
    f = exp(n*atanh(a + b*x))/x**2
    F = -4*b*(-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (a + 1)*(-a - b*x + 1)/((1 - a)*(a + b*x + 1)))/((1 - a)**2*(2 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_878():
    f = exp(n*atanh(a + b*x))/x**3
    F = -2*b**2*(2*a + n)*(-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (a + 1)*(-a - b*x + 1)/((1 - a)*(a + b*x + 1)))/((1 - a)**3*(2 - n)*(a + 1)) - (-a - b*x + 1)**(1 - n/2)*(a + b*x + 1)**(n/2 + 1)/(x**2*(2 - 2*a**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_879():
    f = (-a**2*c*x**2 + c)**4*exp(atanh(a*x))
    F = c**4*x*(-a**2*x**2 + 1)**(sympy.S(7)/2)/8 + 7*c**4*x*(-a**2*x**2 + 1)**(sympy.S(5)/2)/48 + 35*c**4*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/192 + 35*c**4*x*sqrt(-a**2*x**2 + 1)/128 - c**4*(-a**2*x**2 + 1)**(sympy.S(9)/2)/(9*a) + 35*c**4*asin(a*x)/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_880():
    f = (-a**2*c*x**2 + c)**3*exp(atanh(a*x))
    F = c**3*x*(-a**2*x**2 + 1)**(sympy.S(5)/2)/6 + 5*c**3*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/24 + 5*c**3*x*sqrt(-a**2*x**2 + 1)/16 - c**3*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(7*a) + 5*c**3*asin(a*x)/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_881():
    f = (-a**2*c*x**2 + c)**2*exp(atanh(a*x))
    F = c**2*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/4 + 3*c**2*x*sqrt(-a**2*x**2 + 1)/8 - c**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a) + 3*c**2*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_882():
    f = (-a**2*c*x**2 + c)*exp(atanh(a*x))
    F = c*x*sqrt(-a**2*x**2 + 1)/2 - c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a) + c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_883():
    f = x**4*exp(atanh(a*x))/(-a**2*c*x**2 + c)
    F = x**3*(a*x + 1)/(a**2*c*sqrt(-a**2*x**2 + 1)) + 4*x**2*sqrt(-a**2*x**2 + 1)/(3*a**3*c) + (9*a*x + 16)*sqrt(-a**2*x**2 + 1)/(6*a**5*c) - 3*asin(a*x)/(2*a**5*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_884():
    f = x**3*exp(atanh(a*x))/(-a**2*c*x**2 + c)
    F = x**2*(a*x + 1)/(a**2*c*sqrt(-a**2*x**2 + 1)) + (3*a*x + 4)*sqrt(-a**2*x**2 + 1)/(2*a**4*c) - 3*asin(a*x)/(2*a**4*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_885():
    f = x**2*exp(atanh(a*x))/(-a**2*c*x**2 + c)
    F = (a*x + 1)/(a**3*c*sqrt(-a**2*x**2 + 1)) + sqrt(-a**2*x**2 + 1)/(a**3*c) - asin(a*x)/(a**3*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_886():
    f = x*exp(atanh(a*x))/(-a**2*c*x**2 + c)
    F = (a*x + 1)/(a**2*c*sqrt(-a**2*x**2 + 1)) - asin(a*x)/(a**2*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_887():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)
    F = exp(atanh(a*x))/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_888():
    f = exp(atanh(a*x))/(x*(-a**2*c*x**2 + c))
    F = (a*x + 1)/(c*sqrt(-a**2*x**2 + 1)) - atanh(sqrt(-a**2*x**2 + 1))/c
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_889():
    f = exp(atanh(a*x))/(x**2*(-a**2*c*x**2 + c))
    F = -a*atanh(sqrt(-a**2*x**2 + 1))/c + (a*x + 1)/(c*x*sqrt(-a**2*x**2 + 1)) - 2*sqrt(-a**2*x**2 + 1)/(c*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_890():
    f = exp(atanh(a*x))/(x**3*(-a**2*c*x**2 + c))
    F = -3*a**2*atanh(sqrt(-a**2*x**2 + 1))/(2*c) - 2*a*sqrt(-a**2*x**2 + 1)/(c*x) + (a*x + 1)/(c*x**2*sqrt(-a**2*x**2 + 1)) - 3*sqrt(-a**2*x**2 + 1)/(2*c*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_891():
    f = exp(atanh(a*x))/(x**4*(-a**2*c*x**2 + c))
    F = -3*a**3*atanh(sqrt(-a**2*x**2 + 1))/(2*c) - 8*a**2*sqrt(-a**2*x**2 + 1)/(3*c*x) - 3*a*sqrt(-a**2*x**2 + 1)/(2*c*x**2) + (a*x + 1)/(c*x**3*sqrt(-a**2*x**2 + 1)) - 4*sqrt(-a**2*x**2 + 1)/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_892():
    f = x**6*exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = x**5*(a*x + 1)/(3*a**2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x**3*(6*a*x + 5)/(3*a**4*c**2*sqrt(-a**2*x**2 + 1)) - 8*x**2*sqrt(-a**2*x**2 + 1)/(3*a**5*c**2) - (15*a*x + 32)*sqrt(-a**2*x**2 + 1)/(6*a**7*c**2) + 5*asin(a*x)/(2*a**7*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_893():
    f = x**5*exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = x**4*(a*x + 1)/(3*a**2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x**2*(5*a*x + 4)/(3*a**4*c**2*sqrt(-a**2*x**2 + 1)) - (15*a*x + 16)*sqrt(-a**2*x**2 + 1)/(6*a**6*c**2) + 5*asin(a*x)/(2*a**6*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_894():
    f = x**4*exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = x**3*(a*x + 1)/(3*a**2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - x*(4*a*x + 3)/(3*a**4*c**2*sqrt(-a**2*x**2 + 1)) - 8*sqrt(-a**2*x**2 + 1)/(3*a**5*c**2) + asin(a*x)/(a**5*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_895():
    f = x**3*exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = x**2*(a*x + 1)/(3*a**2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - (3*a*x + 2)/(3*a**4*c**2*sqrt(-a**2*x**2 + 1)) + asin(a*x)/(a**4*c**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_896():
    f = x**2*exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = x**2*(a*x + 1)/(3*a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 2/(3*a**3*c**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_897():
    f = x*exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -x/(3*a*c**2*sqrt(-a**2*x**2 + 1)) + (a*x + 1)/(3*a**2*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_898():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = 2*x/(3*c**2*sqrt(-a**2*x**2 + 1)) + (a*x + 1)/(3*a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_899():
    f = exp(atanh(a*x))/(x*(-a**2*c*x**2 + c)**2)
    F = (a*x + 1)/(3*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (2*a*x + 3)/(3*c**2*sqrt(-a**2*x**2 + 1)) - atanh(sqrt(-a**2*x**2 + 1))/c**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_900():
    f = exp(atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**2)
    F = -a*atanh(sqrt(-a**2*x**2 + 1))/c**2 + (a*x + 1)/(3*c**2*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (3*a*x + 4)/(3*c**2*x*sqrt(-a**2*x**2 + 1)) - 8*sqrt(-a**2*x**2 + 1)/(3*c**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_901():
    f = exp(atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**2)
    F = -5*a**2*atanh(sqrt(-a**2*x**2 + 1))/(2*c**2) - 8*a*sqrt(-a**2*x**2 + 1)/(3*c**2*x) + (a*x + 1)/(3*c**2*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (4*a*x + 5)/(3*c**2*x**2*sqrt(-a**2*x**2 + 1)) - 5*sqrt(-a**2*x**2 + 1)/(2*c**2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_902():
    f = exp(atanh(a*x))/(x**4*(-a**2*c*x**2 + c)**2)
    F = -5*a**3*atanh(sqrt(-a**2*x**2 + 1))/(2*c**2) - 16*a**2*sqrt(-a**2*x**2 + 1)/(3*c**2*x) - 5*a*sqrt(-a**2*x**2 + 1)/(2*c**2*x**2) + (a*x + 1)/(3*c**2*x**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (5*a*x + 6)/(3*c**2*x**3*sqrt(-a**2*x**2 + 1)) - 8*sqrt(-a**2*x**2 + 1)/(3*c**2*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_903():
    f = x**7*exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = x**6*(a*x + 1)/(5*a**2*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - x**4*(7*a*x + 6)/(15*a**4*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + x**2*(35*a*x + 24)/(15*a**6*c**3*sqrt(-a**2*x**2 + 1)) + (35*a*x + 32)*sqrt(-a**2*x**2 + 1)/(10*a**8*c**3) - 7*asin(a*x)/(2*a**8*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_904():
    f = x**6*exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = x**5*(a*x + 1)/(5*a**2*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - x**3*(6*a*x + 5)/(15*a**4*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + x*(8*a*x + 5)/(5*a**6*c**3*sqrt(-a**2*x**2 + 1)) + 16*sqrt(-a**2*x**2 + 1)/(5*a**7*c**3) - asin(a*x)/(a**7*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_905():
    f = x**5*exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = x**4*(a*x + 1)/(5*a**2*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - x**2*(5*a*x + 4)/(15*a**4*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (15*a*x + 8)/(15*a**6*c**3*sqrt(-a**2*x**2 + 1)) - asin(a*x)/(a**6*c**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_906():
    f = x**4*exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = x**4*(a*x + 1)/(5*a*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 4/(5*a**5*c**3*sqrt(-a**2*x**2 + 1)) - 4/(15*a**5*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_907():
    f = x**2*exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = x**2*(a*x + 1)/(5*a*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 2*x/(15*a**2*c**3*sqrt(-a**2*x**2 + 1)) - (-2*a*x + 2)/(15*a**3*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_908():
    f = x*exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = -2*x/(15*a*c**3*sqrt(-a**2*x**2 + 1)) - x/(15*a*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (a*x + 1)/(5*a**2*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_909():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = 8*x/(15*c**3*sqrt(-a**2*x**2 + 1)) + 4*x/(15*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (a*x + 1)/(5*a*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_910():
    f = exp(atanh(a*x))/(x*(-a**2*c*x**2 + c)**3)
    F = (a*x + 1)/(5*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + (4*a*x + 5)/(15*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (8*a*x + 15)/(15*c**3*sqrt(-a**2*x**2 + 1)) - atanh(sqrt(-a**2*x**2 + 1))/c**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_911():
    f = exp(atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**3)
    F = -a*atanh(sqrt(-a**2*x**2 + 1))/c**3 + (a*x + 1)/(5*c**3*x*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + (5*a*x + 6)/(15*c**3*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (5*a*x + 8)/(5*c**3*x*sqrt(-a**2*x**2 + 1)) - 16*sqrt(-a**2*x**2 + 1)/(5*c**3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_912():
    f = exp(atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**3)
    F = -7*a**2*atanh(sqrt(-a**2*x**2 + 1))/(2*c**3) - 16*a*sqrt(-a**2*x**2 + 1)/(5*c**3*x) + (a*x + 1)/(5*c**3*x**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + (6*a*x + 7)/(15*c**3*x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + (24*a*x + 35)/(15*c**3*x**2*sqrt(-a**2*x**2 + 1)) - 7*sqrt(-a**2*x**2 + 1)/(2*c**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_913():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = 16*x/(35*c**4*sqrt(-a**2*x**2 + 1)) + 8*x/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 6*x/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + (a*x + 1)/(7*a*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_914():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)**5
    F = 128*x/(315*c**5*sqrt(-a**2*x**2 + 1)) + 64*x/(315*c**5*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 16*x/(105*c**5*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 8*x/(63*c**5*(-a**2*x**2 + 1)**(sympy.S(7)/2)) + (a*x + 1)/(9*a*c**5*(-a**2*x**2 + 1)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_915():
    f = x**4*exp(atanh(a*x))/sqrt(-a**2*x**2 + 1)
    F = -x**4/(4*a) - x**3/(3*a**2) - x**2/(2*a**3) - x/a**4 - log(-a*x + 1)/a**5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_916():
    f = x**3*exp(atanh(a*x))/sqrt(-a**2*x**2 + 1)
    F = -x**3/(3*a) - x**2/(2*a**2) - x/a**3 - log(-a*x + 1)/a**4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_917():
    f = x**2*exp(atanh(a*x))/sqrt(-a**2*x**2 + 1)
    F = -x**2/(2*a) - x/a**2 - log(-a*x + 1)/a**3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_918():
    f = x*exp(atanh(a*x))/sqrt(-a**2*x**2 + 1)
    F = -x/a - log(-a*x + 1)/a**2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_919():
    f = exp(atanh(a*x))/sqrt(-a**2*x**2 + 1)
    F = -log(-a*x + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_920():
    f = exp(atanh(a*x))/(x*sqrt(-a**2*x**2 + 1))
    F = log(x) - log(-a*x + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_921():
    f = exp(atanh(a*x))/(x**2*sqrt(-a**2*x**2 + 1))
    F = a*log(x) - a*log(-a*x + 1) - 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_922():
    f = exp(atanh(a*x))/(x**3*sqrt(-a**2*x**2 + 1))
    F = a**2*log(x) - a**2*log(-a*x + 1) - a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_923():
    f = exp(atanh(a*x))/(x**4*sqrt(-a**2*x**2 + 1))
    F = a**3*log(x) - a**3*log(-a*x + 1) - a**2/x - a/(2*x**2) - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_924():
    f = x**4*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = x**2/(2*a**3) + x/a**4 + 7*log(-a*x + 1)/(4*a**5) + log(a*x + 1)/(4*a**5) + 1/(2*a**5*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_925():
    f = x**3*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = x/a**3 + 5*log(-a*x + 1)/(4*a**4) - log(a*x + 1)/(4*a**4) + 1/(2*a**4*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_926():
    f = x**2*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = 3*log(-a*x + 1)/(4*a**3) + log(a*x + 1)/(4*a**3) + 1/(2*a**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_927():
    f = x*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = -atanh(a*x)/(2*a**2) + 1/(2*a**2*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_928():
    f = exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = atanh(a*x)/(2*a) + 1/(2*a*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_929():
    f = exp(atanh(a*x))/(x*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    F = log(x) - 3*log(-a*x + 1)/4 - log(a*x + 1)/4 + 1/(-2*a*x + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_930():
    f = exp(atanh(a*x))/(x**2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    F = a*log(x) - 5*a*log(-a*x + 1)/4 + a*log(a*x + 1)/4 + a/(-2*a*x + 2) - 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_931():
    f = exp(atanh(a*x))/(x**3*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    F = 2*a**2*log(x) - 7*a**2*log(-a*x + 1)/4 - a**2*log(a*x + 1)/4 + a**2/(-2*a*x + 2) - a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_932():
    f = exp(atanh(a*x))/(x**4*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    F = 2*a**3*log(x) - 9*a**3*log(-a*x + 1)/4 + a**3*log(a*x + 1)/4 + a**3/(-2*a*x + 2) - 2*a**2/x - a/(2*x**2) - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_933():
    f = x**6*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = -x**2/(2*a**5) - x/a**6 - 39*log(-a*x + 1)/(16*a**7) - 9*log(a*x + 1)/(16*a**7) - 1/(8*a**7*(a*x + 1)) - 5/(4*a**7*(-a*x + 1)) + 1/(8*a**7*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_934():
    f = x**5*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = -x/a**5 - 23*log(-a*x + 1)/(16*a**6) + 7*log(a*x + 1)/(16*a**6) + 1/(8*a**6*(a*x + 1)) - 1/(a**6*(-a*x + 1)) + 1/(8*a**6*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_935():
    f = x**4*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = -11*log(-a*x + 1)/(16*a**5) - 5*log(a*x + 1)/(16*a**5) - 1/(8*a**5*(a*x + 1)) - 3/(4*a**5*(-a*x + 1)) + 1/(8*a**5*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_936():
    f = x**3*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = 3*atanh(a*x)/(8*a**4) + 1/(8*a**4*(a*x + 1)) - 1/(2*a**4*(-a*x + 1)) + 1/(8*a**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_937():
    f = x**2*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = -atanh(a*x)/(8*a**3) - 1/(8*a**3*(a*x + 1)) - 1/(4*a**3*(-a*x + 1)) + 1/(8*a**3*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_938():
    f = x*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = -atanh(a*x)/(8*a**2) + 1/(8*a**2*(a*x + 1)) + 1/(8*a**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_939():
    f = exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = 3*atanh(a*x)/(8*a) - 1/(8*a*(a*x + 1)) + 1/(4*a*(-a*x + 1)) + 1/(8*a*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_940():
    f = exp(atanh(a*x))/(x*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    F = log(x) - 11*log(-a*x + 1)/16 - 5*log(a*x + 1)/16 + 1/(8*a*x + 8) + 1/(8*(-a*x + 1)**2) + 1/(-2*a*x + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_941():
    f = exp(atanh(a*x))/(x**2*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    F = a*log(x) - 23*a*log(-a*x + 1)/16 + 7*a*log(a*x + 1)/16 - a/(8*a*x + 8) + a/(8*(-a*x + 1)**2) + 3*a/(-4*a*x + 4) - 1/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_942():
    f = exp(atanh(a*x))/(x**3*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    F = 3*a**2*log(x) - 39*a**2*log(-a*x + 1)/16 - 9*a**2*log(a*x + 1)/16 + a**2/(8*a*x + 8) + a**2/(-a*x + 1) + a**2/(8*(-a*x + 1)**2) - a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_943():
    f = exp(atanh(a*x))/(x**4*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    F = 3*a**3*log(x) - 59*a**3*log(-a*x + 1)/16 + 11*a**3*log(a*x + 1)/16 - a**3/(8*a*x + 8) + a**3/(8*(-a*x + 1)**2) + 5*a**3/(-4*a*x + 4) - 3*a**2/x - a/(2*x**2) - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_944():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(atanh(a*x))
    F = a*x**4*sqrt(-a**2*c*x**2 + c)/(4*sqrt(-a**2*x**2 + 1)) + x**3*sqrt(-a**2*c*x**2 + c)/(3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_945():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(atanh(a*x))
    F = a*x**3*sqrt(-a**2*c*x**2 + c)/(3*sqrt(-a**2*x**2 + 1)) + x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_946():
    f = sqrt(-a**2*c*x**2 + c)*exp(atanh(a*x))
    F = a*x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1)) + x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_947():
    f = sqrt(-a**2*c*x**2 + c)*exp(atanh(a*x))/x
    F = a*x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) + sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_948():
    f = sqrt(-a**2*c*x**2 + c)*exp(atanh(a*x))/x**2
    F = a*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_949():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(atanh(a*x))
    F = -c*(a*x + 1)**4*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(-a**2*x**2 + 1)) + 2*c*(a*x + 1)**3*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_950():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(atanh(a*x))
    F = c**2*(a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(6*a*sqrt(-a**2*x**2 + 1)) - 4*c**2*(a*x + 1)**5*sqrt(-a**2*c*x**2 + c)/(5*a*sqrt(-a**2*x**2 + 1)) + c**2*(a*x + 1)**4*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_951():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(atanh(a*x))
    F = -c**3*(a*x + 1)**8*sqrt(-a**2*c*x**2 + c)/(8*a*sqrt(-a**2*x**2 + 1)) + 6*c**3*(a*x + 1)**7*sqrt(-a**2*c*x**2 + c)/(7*a*sqrt(-a**2*x**2 + 1)) - 2*c**3*(a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1)) + 8*c**3*(a*x + 1)**5*sqrt(-a**2*c*x**2 + c)/(5*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_952():
    f = x**4*exp(atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -x**4*sqrt(-a**2*x**2 + 1)/(4*a*sqrt(-a**2*c*x**2 + c)) - x**3*sqrt(-a**2*x**2 + 1)/(3*a**2*sqrt(-a**2*c*x**2 + c)) - x**2*sqrt(-a**2*x**2 + 1)/(2*a**3*sqrt(-a**2*c*x**2 + c)) - x*sqrt(-a**2*x**2 + 1)/(a**4*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a**5*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_953():
    f = x**3*exp(atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -x**3*sqrt(-a**2*x**2 + 1)/(3*a*sqrt(-a**2*c*x**2 + c)) - x**2*sqrt(-a**2*x**2 + 1)/(2*a**2*sqrt(-a**2*c*x**2 + c)) - x*sqrt(-a**2*x**2 + 1)/(a**3*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a**4*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_954():
    f = x**2*exp(atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -x**2*sqrt(-a**2*x**2 + 1)/(2*a*sqrt(-a**2*c*x**2 + c)) - x*sqrt(-a**2*x**2 + 1)/(a**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_955():
    f = x*exp(atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -x*sqrt(-a**2*x**2 + 1)/(a*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_956():
    f = exp(atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_957():
    f = exp(atanh(a*x))/(x*sqrt(-a**2*c*x**2 + c))
    F = sqrt(-a**2*x**2 + 1)*log(x)/sqrt(-a**2*c*x**2 + c) - sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/sqrt(-a**2*c*x**2 + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_958():
    f = exp(atanh(a*x))/(x**2*sqrt(-a**2*c*x**2 + c))
    F = a*sqrt(-a**2*x**2 + 1)*log(x)/sqrt(-a**2*c*x**2 + c) - a*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/sqrt(-a**2*c*x**2 + c) - sqrt(-a**2*x**2 + 1)/(x*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_959():
    f = exp(atanh(a*x))/(x**3*sqrt(-a**2*c*x**2 + c))
    F = a**2*sqrt(-a**2*x**2 + 1)*log(x)/sqrt(-a**2*c*x**2 + c) - a**2*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/sqrt(-a**2*c*x**2 + c) - a*sqrt(-a**2*x**2 + 1)/(x*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(2*x**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_960():
    f = exp(atanh(a*x))/(x**4*sqrt(-a**2*c*x**2 + c))
    F = a**3*sqrt(-a**2*x**2 + 1)*log(x)/sqrt(-a**2*c*x**2 + c) - a**3*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/sqrt(-a**2*c*x**2 + c) - a**2*sqrt(-a**2*x**2 + 1)/(x*sqrt(-a**2*c*x**2 + c)) - a*sqrt(-a**2*x**2 + 1)/(2*x**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(3*x**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_961():
    f = x**5*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x**3*sqrt(-a**2*x**2 + 1)/(3*a**3*c*sqrt(-a**2*c*x**2 + c)) + x**2*sqrt(-a**2*x**2 + 1)/(2*a**4*c*sqrt(-a**2*c*x**2 + c)) + 2*x*sqrt(-a**2*x**2 + 1)/(a**5*c*sqrt(-a**2*c*x**2 + c)) + 9*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*a**6*c*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*a**6*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*a**6*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_962():
    f = x**4*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x**2*sqrt(-a**2*x**2 + 1)/(2*a**3*c*sqrt(-a**2*c*x**2 + c)) + x*sqrt(-a**2*x**2 + 1)/(a**4*c*sqrt(-a**2*c*x**2 + c)) + 7*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*a**5*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*a**5*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*a**5*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_963():
    f = x**3*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x*sqrt(-a**2*x**2 + 1)/(a**3*c*sqrt(-a**2*c*x**2 + c)) + 5*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*a**4*c*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*a**4*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*a**4*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_964():
    f = x**2*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = 3*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*a**3*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*a**3*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*a**3*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_965():
    f = x*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -sqrt(-a**2*x**2 + 1)*atanh(a*x)/(2*a**2*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*a**2*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_966():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = sqrt(-a**2*x**2 + 1)*atanh(a*x)/(2*a*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*a*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_967():
    f = exp(atanh(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = sqrt(-a**2*x**2 + 1)*log(x)/(c*sqrt(-a**2*c*x**2 + c)) - 3*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_968():
    f = exp(atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = a*sqrt(-a**2*x**2 + 1)*log(x)/(c*sqrt(-a**2*c*x**2 + c)) - 5*a*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) + a*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) + a*sqrt(-a**2*x**2 + 1)/(2*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(c*x*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_969():
    f = exp(atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = 2*a**2*sqrt(-a**2*x**2 + 1)*log(x)/(c*sqrt(-a**2*c*x**2 + c)) - 7*a**2*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) - a**2*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) + a**2*sqrt(-a**2*x**2 + 1)/(2*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) - a*sqrt(-a**2*x**2 + 1)/(c*x*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(2*c*x**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_970():
    f = exp(atanh(a*x))/(x**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = 2*a**3*sqrt(-a**2*x**2 + 1)*log(x)/(c*sqrt(-a**2*c*x**2 + c)) - 9*a**3*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) + a**3*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*c*sqrt(-a**2*c*x**2 + c)) + a**3*sqrt(-a**2*x**2 + 1)/(2*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) - 2*a**2*sqrt(-a**2*x**2 + 1)/(c*x*sqrt(-a**2*c*x**2 + c)) - a*sqrt(-a**2*x**2 + 1)/(2*c*x**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(3*c*x**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_971():
    f = x**6*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -x**2*sqrt(-a**2*x**2 + 1)/(2*a**5*c**2*sqrt(-a**2*c*x**2 + c)) - x*sqrt(-a**2*x**2 + 1)/(a**6*c**2*sqrt(-a**2*c*x**2 + c)) - 39*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(16*a**7*c**2*sqrt(-a**2*c*x**2 + c)) - 9*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(16*a**7*c**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a**7*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - 5*sqrt(-a**2*x**2 + 1)/(4*a**7*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**7*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_972():
    f = x**5*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -x*sqrt(-a**2*x**2 + 1)/(a**5*c**2*sqrt(-a**2*c*x**2 + c)) - 23*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(16*a**6*c**2*sqrt(-a**2*c*x**2 + c)) + 7*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(16*a**6*c**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**6*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(a**6*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**6*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_973():
    f = x**4*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -11*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(16*a**5*c**2*sqrt(-a**2*c*x**2 + c)) - 5*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(16*a**5*c**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a**5*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - 3*sqrt(-a**2*x**2 + 1)/(4*a**5*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**5*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_974():
    f = x**3*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = 3*sqrt(-a**2*x**2 + 1)*atanh(a*x)/(8*a**4*c**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**4*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(2*a**4*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**4*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_975():
    f = x**2*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -sqrt(-a**2*x**2 + 1)*atanh(a*x)/(8*a**3*c**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a**3*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(4*a**3*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**3*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_976():
    f = x*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -sqrt(-a**2*x**2 + 1)*atanh(a*x)/(8*a**2*c**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**2*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a**2*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_977():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = 3*sqrt(-a**2*x**2 + 1)*atanh(a*x)/(8*a*c**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(4*a*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_978():
    f = exp(atanh(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = sqrt(-a**2*x**2 + 1)*log(x)/(c**2*sqrt(-a**2*c*x**2 + c)) - 11*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(16*c**2*sqrt(-a**2*c*x**2 + c)) - 5*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(16*c**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_979():
    f = exp(atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = a*sqrt(-a**2*x**2 + 1)*log(x)/(c**2*sqrt(-a**2*c*x**2 + c)) - 23*a*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(16*c**2*sqrt(-a**2*c*x**2 + c)) + 7*a*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(16*c**2*sqrt(-a**2*c*x**2 + c)) - a*sqrt(-a**2*x**2 + 1)/(8*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) + 3*a*sqrt(-a**2*x**2 + 1)/(4*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + a*sqrt(-a**2*x**2 + 1)/(8*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(c**2*x*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_980():
    f = exp(atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = 3*a**2*sqrt(-a**2*x**2 + 1)*log(x)/(c**2*sqrt(-a**2*c*x**2 + c)) - 39*a**2*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(16*c**2*sqrt(-a**2*c*x**2 + c)) - 9*a**2*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(16*c**2*sqrt(-a**2*c*x**2 + c)) + a**2*sqrt(-a**2*x**2 + 1)/(8*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) + a**2*sqrt(-a**2*x**2 + 1)/(c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + a**2*sqrt(-a**2*x**2 + 1)/(8*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) - a*sqrt(-a**2*x**2 + 1)/(c**2*x*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(2*c**2*x**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_981():
    f = exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = 5*sqrt(-a**2*x**2 + 1)*atanh(a*x)/(16*a*c**3*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a*c**3*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(32*a*c**3*(a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) + 3*sqrt(-a**2*x**2 + 1)/(16*a*c**3*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + 3*sqrt(-a**2*x**2 + 1)/(32*a*c**3*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(24*a*c**3*(-a*x + 1)**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_982():
    f = x**m*(-a**2*c*x**2 + c)**2*exp(atanh(a*x))
    F = a*c**2*x**(m + 2)*hyper((sympy.S(-3)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) + c**2*x**(m + 1)*hyper((sympy.S(-3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_983():
    f = x**m*(-a**2*c*x**2 + c)*exp(atanh(a*x))
    F = a*c*x**(m + 2)*hyper((sympy.S(-1)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) + c*x**(m + 1)*hyper((sympy.S(-1)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_984():
    f = x**m*exp(atanh(a*x))/(-a**2*c*x**2 + c)
    F = a*x**(m + 2)*hyper((sympy.S(3)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(c*(m + 2)) + x**(m + 1)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_985():
    f = x**m*exp(atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = a*x**(m + 2)*hyper((sympy.S(5)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(c**2*(m + 2)) + x**(m + 1)*hyper((sympy.S(5)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(c**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_986():
    f = x**m*exp(atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = a*x**(m + 2)*hyper((sympy.S(7)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(c**3*(m + 2)) + x**(m + 1)*hyper((sympy.S(7)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(c**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_987():
    f = x**m*(-a**2*x**2 + 1)**(sympy.S(5)/2)*exp(atanh(a*x))
    F = a**5*x**(m + 6)/(m + 6) + a**4*x**(m + 5)/(m + 5) - 2*a**3*x**(m + 4)/(m + 4) - 2*a**2*x**(m + 3)/(m + 3) + a*x**(m + 2)/(m + 2) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_988():
    f = x**m*(-a**2*x**2 + 1)**(sympy.S(3)/2)*exp(atanh(a*x))
    F = -a**3*x**(m + 4)/(m + 4) - a**2*x**(m + 3)/(m + 3) + a*x**(m + 2)/(m + 2) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_989():
    f = x**m*sqrt(-a**2*x**2 + 1)*exp(atanh(a*x))
    F = a*x**(m + 2)/(m + 2) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_990():
    f = x**m*exp(atanh(a*x))/sqrt(-a**2*x**2 + 1)
    F = x**(m + 1)*hyper((1, m + 1), (m + 2,), a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_991():
    f = x**m*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = a*x**(m + 2)*hyper((2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) + x**(m + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_992():
    f = x**m*exp(atanh(a*x))/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = a*x**(m + 2)*hyper((3, m/2 + 1), (m/2 + 2,), a**2*x**2)/(m + 2) + x**(m + 1)*hyper((3, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_993():
    f = x**m*(-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(atanh(a*x))
    F = a**5*c**2*x**(m + 6)*sqrt(-a**2*c*x**2 + c)/((m + 6)*sqrt(-a**2*x**2 + 1)) + a**4*c**2*x**(m + 5)*sqrt(-a**2*c*x**2 + c)/((m + 5)*sqrt(-a**2*x**2 + 1)) - 2*a**3*c**2*x**(m + 4)*sqrt(-a**2*c*x**2 + c)/((m + 4)*sqrt(-a**2*x**2 + 1)) - 2*a**2*c**2*x**(m + 3)*sqrt(-a**2*c*x**2 + c)/((m + 3)*sqrt(-a**2*x**2 + 1)) + a*c**2*x**(m + 2)*sqrt(-a**2*c*x**2 + c)/((m + 2)*sqrt(-a**2*x**2 + 1)) + c**2*x**(m + 1)*sqrt(-a**2*c*x**2 + c)/((m + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_994():
    f = x**m*(-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(atanh(a*x))
    F = -a**3*c*x**(m + 4)*sqrt(-a**2*c*x**2 + c)/((m + 4)*sqrt(-a**2*x**2 + 1)) - a**2*c*x**(m + 3)*sqrt(-a**2*c*x**2 + c)/((m + 3)*sqrt(-a**2*x**2 + 1)) + a*c*x**(m + 2)*sqrt(-a**2*c*x**2 + c)/((m + 2)*sqrt(-a**2*x**2 + 1)) + c*x**(m + 1)*sqrt(-a**2*c*x**2 + c)/((m + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_995():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(atanh(a*x))
    F = a*x**(m + 2)*sqrt(-a**2*c*x**2 + c)/((m + 2)*sqrt(-a**2*x**2 + 1)) + x**(m + 1)*sqrt(-a**2*c*x**2 + c)/((m + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_996():
    f = x**m*exp(atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = x**(m + 1)*sqrt(-a**2*x**2 + 1)*hyper((1, m + 1), (m + 2,), a*x)/((m + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_997():
    f = x**m*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = a*x**(m + 2)*sqrt(-a**2*x**2 + 1)*hyper((2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(c*(m + 2)*sqrt(-a**2*c*x**2 + c)) + x**(m + 1)*sqrt(-a**2*x**2 + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(c*(m + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_998():
    f = x**m*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = a*x**(m + 2)*sqrt(-a**2*x**2 + 1)*hyper((3, m/2 + 1), (m/2 + 2,), a**2*x**2)/(c**2*(m + 2)*sqrt(-a**2*c*x**2 + c)) + x**(m + 1)*sqrt(-a**2*x**2 + 1)*hyper((3, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(c**2*(m + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_999():
    f = x**m*(-a**2*c*x**2 + c)**p*exp(atanh(a*x))
    F = a*x**(m + 2)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half - p, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*(-a**2*x**2 + 1)**p) + x**(m + 1)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half - p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1000():
    f = x**3*(-a**2*x**2 + 1)**p*exp(atanh(a*x))
    F = a*x**5*hyper((sympy.S(5)/2, sympy.S.Half - p), (sympy.S(7)/2,), a**2*x**2)/5 + (-a**2*x**2 + 1)**(p + sympy.S(3)/2)/(a**4*(2*p + 3)) - (-a**2*x**2 + 1)**(p + sympy.S.Half)/(a**4*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1001():
    f = x**2*(-a**2*x**2 + 1)**p*exp(atanh(a*x))
    F = x**3*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/3 + (-a**2*x**2 + 1)**(p + sympy.S(3)/2)/(a**3*(2*p + 3)) - (-a**2*x**2 + 1)**(p + sympy.S.Half)/(a**3*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1002():
    f = x*(-a**2*x**2 + 1)**p*exp(atanh(a*x))
    F = a*x**3*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/3 - (-a**2*x**2 + 1)**(p + sympy.S.Half)/(a**2*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1003():
    f = (-a**2*x**2 + 1)**p*exp(atanh(a*x))
    F = -2**(p + sympy.S(3)/2)*(-a*x + 1)**(p + sympy.S.Half)*hyper((p + sympy.S.Half, -p + sympy.S(-1)/2), (p + sympy.S(3)/2,), -a*x/2 + sympy.S.Half)/(a*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1004():
    f = (-a**2*x**2 + 1)**p*exp(atanh(a*x))/x
    F = a*x*hyper((sympy.S.Half, sympy.S.Half - p), (sympy.S(3)/2,), a**2*x**2) - (-a**2*x**2 + 1)**(p + sympy.S.Half)*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1005():
    f = (-a**2*x**2 + 1)**p*exp(atanh(a*x))/x**2
    F = -a*(-a**2*x**2 + 1)**(p + sympy.S.Half)*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) - hyper((sympy.S(-1)/2, sympy.S.Half - p), (sympy.S.Half,), a**2*x**2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1006():
    f = (-a**2*x**2 + 1)**p*exp(atanh(a*x))/x**3
    F = -a**2*(-a**2*x**2 + 1)**(p + sympy.S.Half)*hyper((2, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) - a*hyper((sympy.S(-1)/2, sympy.S.Half - p), (sympy.S.Half,), a**2*x**2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1007():
    f = x**3*(-a**2*c*x**2 + c)**p*exp(atanh(a*x))
    F = a*x**5*(-a**2*c*x**2 + c)**p*hyper((sympy.S(5)/2, sympy.S.Half - p), (sympy.S(7)/2,), a**2*x**2)/(5*(-a**2*x**2 + 1)**p) + (-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a**2*c*x**2 + c)**p/(a**4*(2*p + 3)) - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**4*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1008():
    f = x**2*(-a**2*c*x**2 + c)**p*exp(atanh(a*x))
    F = x**3*(-a**2*c*x**2 + c)**p*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/(3*(-a**2*x**2 + 1)**p) + (-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a**2*c*x**2 + c)**p/(a**3*(2*p + 3)) - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**3*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1009():
    f = x*(-a**2*c*x**2 + c)**p*exp(atanh(a*x))
    F = a*x**3*(-a**2*c*x**2 + c)**p*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/(3*(-a**2*x**2 + 1)**p) - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**2*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1010():
    f = (-a**2*c*x**2 + c)**p*exp(atanh(a*x))
    F = -2**(p + sympy.S(3)/2)*(-a*x + 1)**(p + sympy.S.Half)*(-a**2*c*x**2 + c)**p*hyper((p + sympy.S.Half, -p + sympy.S(-1)/2), (p + sympy.S(3)/2,), -a*x/2 + sympy.S.Half)/(a*(2*p + 1)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1011():
    f = (-a**2*c*x**2 + c)**p*exp(atanh(a*x))/x
    F = a*x*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half, sympy.S.Half - p), (sympy.S(3)/2,), a**2*x**2)/(-a**2*x**2 + 1)**p - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1012():
    f = (-a**2*c*x**2 + c)**p*exp(atanh(a*x))/x**2
    F = -a*sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) - (-a**2*c*x**2 + c)**p*hyper((sympy.S(-1)/2, sympy.S.Half - p), (sympy.S.Half,), a**2*x**2)/(x*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1013():
    f = (-a**2*c*x**2 + c)**p*exp(atanh(a*x))/x**3
    F = -a**2*sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p*hyper((2, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) - a*(-a**2*c*x**2 + c)**p*hyper((sympy.S(-1)/2, sympy.S.Half - p), (sympy.S.Half,), a**2*x**2)/(x*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1014():
    f = x**4*(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = a**2*c*x**7/7 + a*c*x**6/3 + c*x**5/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1015():
    f = x**3*(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = a**2*c*x**6/6 + 2*a*c*x**5/5 + c*x**4/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1016():
    f = x**2*(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = a**2*c*x**5/5 + a*c*x**4/2 + c*x**3/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1017():
    f = x*(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = a**2*c*x**4/4 + 2*a*c*x**3/3 + c*x**2/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1018():
    f = (-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = c*(a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1019():
    f = (-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x
    F = a**2*c*x**2/2 + 2*a*c*x + c*log(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1020():
    f = (-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x**2
    F = a**2*c*x + 2*a*c*log(x) - c/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1021():
    f = (-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x**3
    F = a**2*c*log(x) - 2*a*c/x - c/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1022():
    f = (-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x**4
    F = -c*(a*x + 1)**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1023():
    f = x**4*(-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))
    F = -a**4*c**2*x**9/9 - a**3*c**2*x**8/4 + a*c**2*x**6/3 + c**2*x**5/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1024():
    f = x**3*(-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))
    F = -a**4*c**2*x**8/8 - 2*a**3*c**2*x**7/7 + 2*a*c**2*x**5/5 + c**2*x**4/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1025():
    f = x**2*(-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))
    F = -a**4*c**2*x**7/7 - a**3*c**2*x**6/3 + a*c**2*x**4/2 + c**2*x**3/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1026():
    f = x*(-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))
    F = -a**4*c**2*x**6/6 - 2*a**3*c**2*x**5/5 + 2*a*c**2*x**3/3 + c**2*x**2/2
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1027():
    f = (-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))
    F = -c**2*(a*x + 1)**5/(5*a) + c**2*(a*x + 1)**4/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1028():
    f = (-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))/x
    F = -a**4*c**2*x**4/4 - 2*a**3*c**2*x**3/3 + 2*a*c**2*x + c**2*log(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1029():
    f = (-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))/x**2
    F = -a**4*c**2*x**3/3 - a**3*c**2*x**2 + 2*a*c**2*log(x) - c**2/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1030():
    f = (-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))/x**3
    F = -c**2*(a*x + 1)**4/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1031():
    f = (-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))/x**4
    F = -a**4*c**2*x - 2*a**3*c**2*log(x) - a*c**2/x**2 - c**2/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1032():
    f = (-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))/x**5
    F = -a**4*c**2*log(x) + 2*a**3*c**2/x - 2*a*c**2/(3*x**3) - c**2/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1033():
    f = (-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))/x**6
    F = a**4*c**2/x + a**3*c**2/x**2 - a*c**2/(2*x**4) - c**2/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1034():
    f = x**4*(-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))
    F = a**6*c**3*x**11/11 + a**5*c**3*x**10/5 - a**4*c**3*x**9/9 - a**3*c**3*x**8/2 - a**2*c**3*x**7/7 + a*c**3*x**6/3 + c**3*x**5/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1035():
    f = x**3*(-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))
    F = a**6*c**3*x**10/10 + 2*a**5*c**3*x**9/9 - a**4*c**3*x**8/8 - 4*a**3*c**3*x**7/7 - a**2*c**3*x**6/6 + 2*a*c**3*x**5/5 + c**3*x**4/4
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1036():
    f = x**2*(-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))
    F = c**3*(a*x + 1)**9/(9*a**3) - 3*c**3*(a*x + 1)**8/(4*a**3) + 13*c**3*(a*x + 1)**7/(7*a**3) - 2*c**3*(a*x + 1)**6/a**3 + 4*c**3*(a*x + 1)**5/(5*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1037():
    f = x*(-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))
    F = c**3*(a*x + 1)**8/(8*a**2) - 5*c**3*(a*x + 1)**7/(7*a**2) + 4*c**3*(a*x + 1)**6/(3*a**2) - 4*c**3*(a*x + 1)**5/(5*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1038():
    f = (-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))
    F = c**3*(a*x + 1)**7/(7*a) - 2*c**3*(a*x + 1)**6/(3*a) + 4*c**3*(a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1039():
    f = (-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))/x
    F = a**6*c**3*x**6/6 + 2*a**5*c**3*x**5/5 - a**4*c**3*x**4/4 - 4*a**3*c**3*x**3/3 - a**2*c**3*x**2/2 + 2*a*c**3*x + c**3*log(x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1040():
    f = (-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))/x**2
    F = a**6*c**3*x**5/5 + a**5*c**3*x**4/2 - a**4*c**3*x**3/3 - 2*a**3*c**3*x**2 - a**2*c**3*x + 2*a*c**3*log(x) - c**3/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1041():
    f = (-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))/x**3
    F = a**6*c**3*x**4/4 + 2*a**5*c**3*x**3/3 - a**4*c**3*x**2/2 - 4*a**3*c**3*x - a**2*c**3*log(x) - 2*a*c**3/x - c**3/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1042():
    f = (-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))/x**4
    F = a**6*c**3*x**3/3 + a**5*c**3*x**2 - a**4*c**3*x - 4*a**3*c**3*log(x) + a**2*c**3/x - a*c**3/x**2 - c**3/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1043():
    f = (-a**2*c*x**2 + c)**4*exp(2*atanh(a*x))
    F = -c**4*(a*x + 1)**9/(9*a) + 3*c**4*(a*x + 1)**8/(4*a) - 12*c**4*(a*x + 1)**7/(7*a) + 4*c**4*(a*x + 1)**6/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1044():
    f = x**4*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)
    F = x**3/(3*a**2*c) + x**2/(a**3*c) + 3*x/(a**4*c) + 4*log(-a*x + 1)/(a**5*c) + 1/(a**5*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1045():
    f = x**3*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)
    F = x**2/(2*a**2*c) + 2*x/(a**3*c) + 3*log(-a*x + 1)/(a**4*c) + 1/(a**4*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1046():
    f = x**2*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)
    F = x/(a**2*c) + 2*log(-a*x + 1)/(a**3*c) + 1/(a**3*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1047():
    f = x*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)
    F = log(-a*x + 1)/(a**2*c) + 1/(a**2*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1048():
    f = exp(2*atanh(a*x))/(-a**2*c*x**2 + c)
    F = 1/(a*c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1049():
    f = exp(2*atanh(a*x))/(x*(-a**2*c*x**2 + c))
    F = log(x)/c - log(-a*x + 1)/c + 1/(c*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1050():
    f = exp(2*atanh(a*x))/(x**2*(-a**2*c*x**2 + c))
    F = 2*a*log(x)/c - 2*a*log(-a*x + 1)/c + a/(c*(-a*x + 1)) - 1/(c*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1051():
    f = exp(2*atanh(a*x))/(x**3*(-a**2*c*x**2 + c))
    F = 3*a**2*log(x)/c - 3*a**2*log(-a*x + 1)/c + a**2/(c*(-a*x + 1)) - 2*a/(c*x) - 1/(2*c*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1052():
    f = exp(2*atanh(a*x))/(x**4*(-a**2*c*x**2 + c))
    F = 4*a**3*log(x)/c - 4*a**3*log(-a*x + 1)/c + a**3/(c*(-a*x + 1)) - 3*a**2/(c*x) - a/(c*x**2) - 1/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1053():
    f = x**4*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -x/(a**4*c**2) - 17*log(-a*x + 1)/(8*a**5*c**2) + log(a*x + 1)/(8*a**5*c**2) - 7/(4*a**5*c**2*(-a*x + 1)) + 1/(4*a**5*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1054():
    f = x**3*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -7*log(-a*x + 1)/(8*a**4*c**2) - log(a*x + 1)/(8*a**4*c**2) - 5/(4*a**4*c**2*(-a*x + 1)) + 1/(4*a**4*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1055():
    f = x**2*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = atanh(a*x)/(4*a**3*c**2) - 3/(4*a**3*c**2*(-a*x + 1)) + 1/(4*a**3*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1056():
    f = x*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -atanh(a*x)/(4*a**2*c**2) - 1/(4*a**2*c**2*(-a*x + 1)) + 1/(4*a**2*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1057():
    f = exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = atanh(a*x)/(4*a*c**2) + 1/(4*a*c**2*(-a*x + 1)) + 1/(4*a*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1058():
    f = exp(2*atanh(a*x))/(x*(-a**2*c*x**2 + c)**2)
    F = log(x)/c**2 - 7*log(-a*x + 1)/(8*c**2) - log(a*x + 1)/(8*c**2) + 3/(4*c**2*(-a*x + 1)) + 1/(4*c**2*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1059():
    f = exp(2*atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**2)
    F = 2*a*log(x)/c**2 - 17*a*log(-a*x + 1)/(8*c**2) + a*log(a*x + 1)/(8*c**2) + 5*a/(4*c**2*(-a*x + 1)) + a/(4*c**2*(-a*x + 1)**2) - 1/(c**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1060():
    f = exp(2*atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**2)
    F = 4*a**2*log(x)/c**2 - 31*a**2*log(-a*x + 1)/(8*c**2) - a**2*log(a*x + 1)/(8*c**2) + 7*a**2/(4*c**2*(-a*x + 1)) + a**2/(4*c**2*(-a*x + 1)**2) - 2*a/(c**2*x) - 1/(2*c**2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1061():
    f = exp(2*atanh(a*x))/(x**4*(-a**2*c*x**2 + c)**2)
    F = 6*a**3*log(x)/c**2 - 49*a**3*log(-a*x + 1)/(8*c**2) + a**3*log(a*x + 1)/(8*c**2) + 9*a**3/(4*c**2*(-a*x + 1)) + a**3/(4*c**2*(-a*x + 1)**2) - 4*a**2/(c**2*x) - a/(c**2*x**2) - 1/(3*c**2*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1062():
    f = x**5*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = 13*log(-a*x + 1)/(16*a**6*c**3) + 3*log(a*x + 1)/(16*a**6*c**3) + 1/(16*a**6*c**3*(a*x + 1)) + 23/(16*a**6*c**3*(-a*x + 1)) - 1/(2*a**6*c**3*(-a*x + 1)**2) + 1/(12*a**6*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1063():
    f = x**4*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = -atanh(a*x)/(4*a**5*c**3) - 1/(16*a**5*c**3*(a*x + 1)) + 11/(16*a**5*c**3*(-a*x + 1)) - 3/(8*a**5*c**3*(-a*x + 1)**2) + 1/(12*a**5*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1064():
    f = x**3*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = atanh(a*x)/(8*a**4*c**3) + 1/(16*a**4*c**3*(a*x + 1)) + 3/(16*a**4*c**3*(-a*x + 1)) - 1/(4*a**4*c**3*(-a*x + 1)**2) + 1/(12*a**4*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1065():
    f = x**2*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = -(-2*a*x + 1)/(6*a**3*c**3*(-a*x + 1)**3*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1066():
    f = x*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = -atanh(a*x)/(8*a**2*c**3) + 1/(16*a**2*c**3*(a*x + 1)) - 1/(16*a**2*c**3*(-a*x + 1)) + 1/(12*a**2*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1067():
    f = exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = atanh(a*x)/(4*a*c**3) - 1/(16*a*c**3*(a*x + 1)) + 3/(16*a*c**3*(-a*x + 1)) + 1/(8*a*c**3*(-a*x + 1)**2) + 1/(12*a*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1068():
    f = exp(2*atanh(a*x))/(x*(-a**2*c*x**2 + c)**3)
    F = log(x)/c**3 - 13*log(-a*x + 1)/(16*c**3) - 3*log(a*x + 1)/(16*c**3) + 1/(16*c**3*(a*x + 1)) + 11/(16*c**3*(-a*x + 1)) + 1/(4*c**3*(-a*x + 1)**2) + 1/(12*c**3*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1069():
    f = exp(2*atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**3)
    F = 2*a*log(x)/c**3 - 9*a*log(-a*x + 1)/(4*c**3) + a*log(a*x + 1)/(4*c**3) - a/(16*c**3*(a*x + 1)) + 23*a/(16*c**3*(-a*x + 1)) + 3*a/(8*c**3*(-a*x + 1)**2) + a/(12*c**3*(-a*x + 1)**3) - 1/(c**3*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1070():
    f = exp(2*atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**3)
    F = 5*a**2*log(x)/c**3 - 75*a**2*log(-a*x + 1)/(16*c**3) - 5*a**2*log(a*x + 1)/(16*c**3) + a**2/(16*c**3*(a*x + 1)) + 39*a**2/(16*c**3*(-a*x + 1)) + a**2/(2*c**3*(-a*x + 1)**2) + a**2/(12*c**3*(-a*x + 1)**3) - 2*a/(c**3*x) - 1/(2*c**3*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1071():
    f = exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = 15*atanh(a*x)/(64*a*c**4) - 5/(64*a*c**4*(a*x + 1)) - 1/(64*a*c**4*(a*x + 1)**2) + 5/(32*a*c**4*(-a*x + 1)) + 3/(32*a*c**4*(-a*x + 1)**2) + 1/(16*a*c**4*(-a*x + 1)**3) + 1/(32*a*c**4*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1072():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = -x**4*sqrt(-a**2*c*x**2 + c)/5 - x**3*sqrt(-a**2*c*x**2 + c)/(2*a) - 3*x**2*sqrt(-a**2*c*x**2 + c)/(5*a**2) + 3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(4*a**4) - (15*a*x + 24)*sqrt(-a**2*c*x**2 + c)/(20*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1073():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = -x**3*sqrt(-a**2*c*x**2 + c)/4 - 2*x**2*sqrt(-a**2*c*x**2 + c)/(3*a) + 7*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a**3) - (21*a*x + 32)*sqrt(-a**2*c*x**2 + c)/(24*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1074():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = -x**2*sqrt(-a**2*c*x**2 + c)/3 + sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/a**2 - (3*a*x + 5)*sqrt(-a**2*c*x**2 + c)/(3*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1075():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = 3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(2*a) - (a*x + 1)*sqrt(-a**2*c*x**2 + c)/(2*a) - 3*sqrt(-a**2*c*x**2 + c)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1076():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x
    F = 2*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) - sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - sqrt(-a**2*c*x**2 + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1077():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x**2
    F = a*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) - 2*a*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - sqrt(-a**2*c*x**2 + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1078():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x**3
    F = -3*a**2*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/2 - 2*a*sqrt(-a**2*c*x**2 + c)/x - sqrt(-a**2*c*x**2 + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1079():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x**4
    F = -a**3*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - 5*a**2*sqrt(-a**2*c*x**2 + c)/(3*x) - a*sqrt(-a**2*c*x**2 + c)/x**2 - sqrt(-a**2*c*x**2 + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1080():
    f = sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))/x**5
    F = -7*a**4*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/8 - 4*a**3*sqrt(-a**2*c*x**2 + c)/(3*x) - 7*a**2*sqrt(-a**2*c*x**2 + c)/(8*x**2) - 2*a*sqrt(-a**2*c*x**2 + c)/(3*x**3) - sqrt(-a**2*c*x**2 + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1081():
    f = x**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = -x**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/7 - x**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(3*a) - 11*x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(35*a**2) + c*x*sqrt(-a**2*c*x**2 + c)/(8*a**3) + c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a**4) - (105*a*x + 88)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(420*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1082():
    f = x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = -x**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/6 - 2*x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(5*a) + 3*c*x*sqrt(-a**2*c*x**2 + c)/(16*a**2) + 3*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(16*a**3) - (45*a*x + 32)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(120*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1083():
    f = x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = -x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/5 + c*x*sqrt(-a**2*c*x**2 + c)/(4*a) + c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(4*a**2) - (15*a*x + 14)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(30*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1084():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = 5*c*x*sqrt(-a**2*c*x**2 + c)/8 + 5*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a) - (a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a) - 5*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(12*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1085():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x
    F = c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) - c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + c*(a*x + 1)*sqrt(-a**2*c*x**2 + c) - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/3
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1086():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x**2
    F = -a*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/2 - 2*a*c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + a*c*(-a*x + 4)*sqrt(-a**2*c*x**2 + c)/2 - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1087():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x**3
    F = -2*a**2*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) - a**2*c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/2 - a*c*(-a*x + 4)*sqrt(-a**2*c*x**2 + c)/(2*x) - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1088():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x**4
    F = -a**3*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) + a**3*c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - a*c*(a*x + 1)*sqrt(-a**2*c*x**2 + c)/x**2 - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1089():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x**5
    F = 5*a**4*c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/8 - 5*a**2*c*sqrt(-a**2*c*x**2 + c)/(8*x**2) - 2*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(3*x**3) - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1090():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x**6
    F = a**5*c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/4 - a**3*c*sqrt(-a**2*c*x**2 + c)/(4*x**2) - 7*a**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(15*x**3) - a*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(2*x**4) - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1091():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x**7
    F = 3*a**6*c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/16 - 3*a**4*c*sqrt(-a**2*c*x**2 + c)/(16*x**2) - 4*a**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(15*x**3) - 3*a**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(8*x**4) - 2*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(5*x**5) - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/(6*x**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1092():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))/x**8
    F = a**7*c**(sympy.S(3)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/8 - a**5*c*sqrt(-a**2*c*x**2 + c)/(8*x**2) - 22*a**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(105*x**3) - a**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*x**4) - 11*a**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(35*x**5) - a*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(3*x**6) - (-a**2*c*x**2 + c)**(sympy.S(3)/2)/(7*x**7)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1093():
    f = x**3*(-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = -x**4*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/9 - x**3*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(4*a) - 13*x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(63*a**2) + 3*c**2*x*sqrt(-a**2*c*x**2 + c)/(64*a**3) + c*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(32*a**3) + 3*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(64*a**4) - (315*a*x + 208)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(2520*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1094():
    f = x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = -x**3*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/8 - 2*x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(7*a) + 11*c**2*x*sqrt(-a**2*c*x**2 + c)/(128*a**2) + 11*c*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(192*a**2) + 11*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(128*a**3) - (385*a*x + 192)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(1680*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1095():
    f = x*(-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = -x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/7 + c**2*x*sqrt(-a**2*c*x**2 + c)/(8*a) + c*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(12*a) + c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a**2) - (35*a*x + 27)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(105*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1096():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = 7*c**2*x*sqrt(-a**2*c*x**2 + c)/16 + 7*c*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/24 + 7*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(16*a) - (a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a) - 7*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(30*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1097():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))/x
    F = 3*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/4 - c**(sympy.S(5)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + c**2*(3*a*x + 4)*sqrt(-a**2*c*x**2 + c)/4 + c*(3*a*x + 2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/6 - (-a**2*c*x**2 + c)**(sympy.S(5)/2)/5
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1098():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))/x**2
    F = -9*a*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/8 - 2*a*c**(sympy.S(5)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) + a*c**2*(-9*a*x + 16)*sqrt(-a**2*c*x**2 + c)/8 + a*c*(-9*a*x + 8)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/12 - (-a**2*c*x**2 + c)**(sympy.S(5)/2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1099():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))/x**3
    F = -3*a**2*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) + a**2*c**(sympy.S(5)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/2 - a**2*c**2*(6*a*x + 1)*sqrt(-a**2*c*x**2 + c)/2 - a*c*(a*x + 12)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(6*x) - (-a**2*c*x**2 + c)**(sympy.S(5)/2)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1100():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))/x**4
    F = -a**3*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/2 + 3*a**3*c**(sympy.S(5)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - a**2*c**2*(6*a*x + 1)*sqrt(-a**2*c*x**2 + c)/(2*x) - a*c*(-a*x + 6)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(6*x**2) - (-a**2*c*x**2 + c)**(sympy.S(5)/2)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1101():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))/x**5
    F = 2*a**4*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) + 9*a**4*c**(sympy.S(5)/2)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/8 + a**3*c**2*(-9*a*x + 16)*sqrt(-a**2*c*x**2 + c)/(8*x) - a*c*(9*a*x + 16)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(24*x**3) - (-a**2*c*x**2 + c)**(sympy.S(5)/2)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1102():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(2*atanh(a*x))
    F = 45*c**3*x*sqrt(-a**2*c*x**2 + c)/128 + 15*c**2*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/64 + 3*c*x*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/16 + 45*c**(sympy.S(7)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(128*a) - (a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(8*a) - 9*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(56*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1103():
    f = x**3*exp(2*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = x**2*sqrt(-a**2*c*x**2 + c)/(3*a**2*c) + x*sqrt(-a**2*c*x**2 + c)/(a**3*c) + (a*x + 1)**2/(a**4*sqrt(-a**2*c*x**2 + c)) + 11*sqrt(-a**2*c*x**2 + c)/(3*a**4*c) - 3*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a**4*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1104():
    f = x**2*exp(2*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = (a*x + 1)**2/(a**3*sqrt(-a**2*c*x**2 + c)) + (a*x + 6)*sqrt(-a**2*c*x**2 + c)/(2*a**3*c) - 5*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(2*a**3*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1105():
    f = x*exp(2*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = (a*x + 1)**2/(a**2*sqrt(-a**2*c*x**2 + c)) + 2*sqrt(-a**2*c*x**2 + c)/(a**2*c) - 2*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a**2*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1106():
    f = exp(2*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = (2*a*x + 2)/(a*sqrt(-a**2*c*x**2 + c)) - atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1107():
    f = exp(2*atanh(a*x))/(x*sqrt(-a**2*c*x**2 + c))
    F = (2*a*x + 2)/sqrt(-a**2*c*x**2 + c) - atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/sqrt(c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1108():
    f = exp(2*atanh(a*x))/(x**2*sqrt(-a**2*c*x**2 + c))
    F = 2*a*(a*x + 1)/sqrt(-a**2*c*x**2 + c) - 2*a*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/sqrt(c) - sqrt(-a**2*c*x**2 + c)/(c*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1109():
    f = exp(2*atanh(a*x))/(x**3*sqrt(-a**2*c*x**2 + c))
    F = 2*a**2*(a*x + 1)/sqrt(-a**2*c*x**2 + c) - 5*a**2*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/(2*sqrt(c)) - 2*a*sqrt(-a**2*c*x**2 + c)/(c*x) - sqrt(-a**2*c*x**2 + c)/(2*c*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1110():
    f = exp(2*atanh(a*x))/(x**4*sqrt(-a**2*c*x**2 + c))
    F = 2*a**3*(a*x + 1)/sqrt(-a**2*c*x**2 + c) - 3*a**3*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/sqrt(c) - 8*a**2*sqrt(-a**2*c*x**2 + c)/(3*c*x) - a*sqrt(-a**2*c*x**2 + c)/(c*x**2) - sqrt(-a**2*c*x**2 + c)/(3*c*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1111():
    f = x**3*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = (a*x + 1)**2/(3*a**4*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - (8*a*x + 8)/(3*a**4*c*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*c*x**2 + c)/(a**4*c**2) + 2*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a**4*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1112():
    f = x**2*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = (a*x + 1)**2/(3*a**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - (5*a*x + 5)/(3*a**3*c*sqrt(-a**2*c*x**2 + c)) + atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a**3*c**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1113():
    f = x*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = (a*x + 1)**2/(3*a**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - (2*a*x + 2)/(3*a**2*c*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1114():
    f = exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x/(3*c*sqrt(-a**2*c*x**2 + c)) + (2*a*x + 2)/(3*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1115():
    f = exp(2*atanh(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = (2*a*x + 2)/(3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + (4*a*x + 3)/(3*c*sqrt(-a**2*c*x**2 + c)) - atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/c**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1116():
    f = exp(2*atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = 2*a*(a*x + 1)/(3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a*(7*a*x + 6)/(3*c*sqrt(-a**2*c*x**2 + c)) - 2*a*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/c**(sympy.S(3)/2) - sqrt(-a**2*c*x**2 + c)/(c**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1117():
    f = exp(2*atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = 2*a**2*(a*x + 1)/(3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + a**2*(10*a*x + 9)/(3*c*sqrt(-a**2*c*x**2 + c)) - 7*a**2*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/(2*c**(sympy.S(3)/2)) - 2*a*sqrt(-a**2*c*x**2 + c)/(c**2*x) - sqrt(-a**2*c*x**2 + c)/(2*c**2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1118():
    f = exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = x/(5*c*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + 2*x/(5*c**2*sqrt(-a**2*c*x**2 + c)) + (2*a*x + 2)/(5*a*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1119():
    f = exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = x/(7*c*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + 4*x/(21*c**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + 8*x/(21*c**3*sqrt(-a**2*c*x**2 + c)) + (2*a*x + 2)/(7*a*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1120():
    f = x**m*(-a**2*c*x**2 + c)**3*exp(2*atanh(a*x))
    F = a**6*c**3*x**(m + 7)/(m + 7) + 2*a**5*c**3*x**(m + 6)/(m + 6) - a**4*c**3*x**(m + 5)/(m + 5) - 4*a**3*c**3*x**(m + 4)/(m + 4) - a**2*c**3*x**(m + 3)/(m + 3) + 2*a*c**3*x**(m + 2)/(m + 2) + c**3*x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1121():
    f = x**m*(-a**2*c*x**2 + c)**2*exp(2*atanh(a*x))
    F = -a**4*c**2*x**(m + 5)/(m + 5) - 2*a**3*c**2*x**(m + 4)/(m + 4) + 2*a*c**2*x**(m + 2)/(m + 2) + c**2*x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1122():
    f = x**m*(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = a**2*c*x**(m + 3)/(m + 3) + 2*a*c*x**(m + 2)/(m + 2) + c*x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1123():
    f = x**m*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)
    F = x**(m + 1)*hyper((2, m + 1), (m + 2,), a*x)/(c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1124():
    f = x**m*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = x**(m + 1)*(2 - m)/(4*c**2*(-a*x + 1)) + x**(m + 1)/(4*c**2*(-a*x + 1)**2) + x**(m + 1)*(2*m**2 - 4*m + 1)*hyper((1, m + 1), (m + 2,), a*x)/(8*c**2*(m + 1)) + x**(m + 1)*hyper((1, m + 1), (m + 2,), -a*x)/(8*c**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1125():
    f = x**m*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = -x**(m + 1)*(2 - m)*(4 - m)/(24*c**3*(a*x + 1)) + x**(m + 1)*(2 - m)*(7 - 2*m)/(24*c**3*(-a*x + 1)*(a*x + 1)) + x**(m + 1)*(2 - m)*(2*m**2 - 8*m + 3)*hyper((1, m + 1), (m + 2,), a*x)/(48*c**3*(m + 1)) + x**(m + 1)*(2 - m)*hyper((1, m + 1), (m + 2,), -a*x)/(16*c**3*(m + 1)) + x**(m + 1)*(4 - m)/(12*c**3*(-a*x + 1)**2*(a*x + 1)) + x**(m + 1)/(6*c**3*(-a*x + 1)**3*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1126():
    f = x**m*(-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(2*atanh(a*x))
    F = 2*a*c**2*x**(m + 2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(-3)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*sqrt(-a**2*x**2 + 1)) + c**2*x**(m + 1)*(2*m + 7)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(-3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(m + 6)*sqrt(-a**2*x**2 + 1)) - x**(m + 1)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(m + 6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1127():
    f = x**m*(-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atanh(a*x))
    F = 2*a*c*x**(m + 2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(-1)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*sqrt(-a**2*x**2 + 1)) + c*x**(m + 1)*(2*m + 5)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(-1)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(m + 4)*sqrt(-a**2*x**2 + 1)) - x**(m + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(m + 4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1128():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(2*atanh(a*x))
    F = 2*a*c*x**(m + 2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*sqrt(-a**2*c*x**2 + c)) + c*x**(m + 1)*(2*m + 3)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(m + 2)*sqrt(-a**2*c*x**2 + c)) - x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(m + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1129():
    f = x**m*exp(2*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -2*a*x**(m + 2)*(m + 1)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*sqrt(-a**2*c*x**2 + c)) + 2*x**(m + 1)*(a*x + 1)/sqrt(-a**2*c*x**2 + c) - x**(m + 1)*(2*m + 1)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1130():
    f = x**m*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = 2*a*x**(m + 2)*(1 - m)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S(3)/2, m/2 + 1), (m/2 + 2,), a**2*x**2)/(3*c*(m + 2)*sqrt(-a**2*c*x**2 + c)) + 2*x**(m + 1)*(a*x + 1)/(3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + x**(m + 1)*(1 - 2*m)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/(3*c*(m + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1131():
    f = (-a**2*c*x**2 + c)**p*exp(2*atanh(a*x))
    F = -2**(p + 1)*(-a**2*c*x**2 + c)**p*hyper((p, -p - 1), (p + 1,), -a*x/2 + sympy.S.Half)/(a*p*(a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1132():
    f = x**3*(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -a*c*x**5*sqrt(-a**2*x**2 + 1)/6 - 3*c*x**4*sqrt(-a**2*x**2 + 1)/5 - 23*c*x**3*sqrt(-a**2*x**2 + 1)/(24*a) - 17*c*x**2*sqrt(-a**2*x**2 + 1)/(15*a**2) - c*(345*a*x + 544)*sqrt(-a**2*x**2 + 1)/(240*a**4) + 23*c*asin(a*x)/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1133():
    f = x**2*(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -a*c*x**4*sqrt(-a**2*x**2 + 1)/5 - 3*c*x**3*sqrt(-a**2*x**2 + 1)/4 - 19*c*x**2*sqrt(-a**2*x**2 + 1)/(15*a) - c*(195*a*x + 304)*sqrt(-a**2*x**2 + 1)/(120*a**3) + 13*c*asin(a*x)/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1134():
    f = x*(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -c*(a*x + 1)**3*sqrt(-a**2*x**2 + 1)/(4*a**2) - c*(a*x + 1)**2*sqrt(-a**2*x**2 + 1)/(4*a**2) - 5*c*(a*x + 1)*sqrt(-a**2*x**2 + 1)/(8*a**2) - 15*c*sqrt(-a**2*x**2 + 1)/(8*a**2) + 15*c*asin(a*x)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1135():
    f = (-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -c*(a*x + 1)**2*sqrt(-a**2*x**2 + 1)/(3*a) - 5*c*(a*x + 1)*sqrt(-a**2*x**2 + 1)/(6*a) - 5*c*sqrt(-a**2*x**2 + 1)/(2*a) + 5*c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1136():
    f = (-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x
    F = -a*c*x*sqrt(-a**2*x**2 + 1)/2 - 3*c*sqrt(-a**2*x**2 + 1) + 7*c*asin(a*x)/2 - c*atanh(sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1137():
    f = (-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**2
    F = -a*c*sqrt(-a**2*x**2 + 1) + 3*a*c*asin(a*x) - 3*a*c*atanh(sqrt(-a**2*x**2 + 1)) - c*sqrt(-a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1138():
    f = (-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**3
    F = a**2*c*asin(a*x) - 7*a**2*c*atanh(sqrt(-a**2*x**2 + 1))/2 - 3*a*c*sqrt(-a**2*x**2 + 1)/x - c*sqrt(-a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1139():
    f = (-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**4
    F = -5*a**3*c*atanh(sqrt(-a**2*x**2 + 1))/2 - 11*a**2*c*sqrt(-a**2*x**2 + 1)/(3*x) - 3*a*c*sqrt(-a**2*x**2 + 1)/(2*x**2) - c*sqrt(-a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1140():
    f = (-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**5
    F = -15*a**4*c*atanh(sqrt(-a**2*x**2 + 1))/8 - 3*a**3*c*sqrt(-a**2*x**2 + 1)/x - 15*a**2*c*sqrt(-a**2*x**2 + 1)/(8*x**2) - a*c*sqrt(-a**2*x**2 + 1)/x**3 - c*sqrt(-a**2*x**2 + 1)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1141():
    f = (-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**6
    F = -13*a**5*c*atanh(sqrt(-a**2*x**2 + 1))/8 - 38*a**4*c*sqrt(-a**2*x**2 + 1)/(15*x) - 13*a**3*c*sqrt(-a**2*x**2 + 1)/(8*x**2) - 19*a**2*c*sqrt(-a**2*x**2 + 1)/(15*x**3) - 3*a*c*sqrt(-a**2*x**2 + 1)/(4*x**4) - c*sqrt(-a**2*x**2 + 1)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1142():
    f = (-a**2*c*x**2 + c)**2*exp(3*atanh(a*x))
    F = 7*c**2*x*sqrt(-a**2*x**2 + 1)/8 - c**2*(a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a) - 7*c**2*(a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(20*a) - 7*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(12*a) + 7*c**2*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1143():
    f = (-a**2*c*x**2 + c)**3*exp(3*atanh(a*x))
    F = 3*c**3*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/8 + 9*c**3*x*sqrt(-a**2*x**2 + 1)/16 - c**3*(a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(7*a) - 3*c**3*(a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(14*a) - 3*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(10*a) + 9*c**3*asin(a*x)/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1144():
    f = (-a**2*c*x**2 + c)**4*exp(3*atanh(a*x))
    F = 11*c**4*x*(-a**2*x**2 + 1)**(sympy.S(5)/2)/48 + 55*c**4*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/192 + 55*c**4*x*sqrt(-a**2*x**2 + 1)/128 - c**4*(a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(9*a) - 11*c**4*(a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(72*a) - 11*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(56*a) + 55*c**4*asin(a*x)/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1145():
    f = x**2*exp(3*atanh(a*x))/(-a**2*c*x**2 + c)
    F = (a*x + 1)**3/(3*a**3*c*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 2*(a*x + 1)**2/(a**3*c*sqrt(-a**2*x**2 + 1)) - 3*sqrt(-a**2*x**2 + 1)/(a**3*c) + 3*asin(a*x)/(a**3*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1146():
    f = x*exp(3*atanh(a*x))/(-a**2*c*x**2 + c)
    F = (a*x + 1)**3/(3*a**2*c*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - (2*a*x + 2)/(a**2*c*sqrt(-a**2*x**2 + 1)) + asin(a*x)/(a**2*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1147():
    f = exp(3*atanh(a*x))/(-a**2*c*x**2 + c)
    F = exp(3*atanh(a*x))/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1148():
    f = exp(3*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = 2*sqrt(-a**2*x**2 + 1)/(15*a*c**2*(-a*x + 1)) + 2*sqrt(-a**2*x**2 + 1)/(15*a*c**2*(-a*x + 1)**2) + sqrt(-a**2*x**2 + 1)/(5*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1149():
    f = exp(3*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = 8*x/(35*c**3*sqrt(-a**2*x**2 + 1)) + 4/(35*a*c**3*(-a*x + 1)*sqrt(-a**2*x**2 + 1)) + 4/(35*a*c**3*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1)) + 1/(7*a*c**3*(-a*x + 1)**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1150():
    f = exp(3*atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = 16*x/(63*c**4*sqrt(-a**2*x**2 + 1)) + 8*x/(63*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 2/(21*a*c**4*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 2/(21*a*c**4*(-a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 1/(9*a*c**4*(-a*x + 1)**3*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1151():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -a*x**5*sqrt(-a**2*c*x**2 + c)/(5*sqrt(-a**2*x**2 + 1)) - 3*x**4*sqrt(-a**2*c*x**2 + c)/(4*sqrt(-a**2*x**2 + 1)) - 4*x**3*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1)) - 2*x**2*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(-a**2*x**2 + 1)) - 4*x*sqrt(-a**2*c*x**2 + c)/(a**3*sqrt(-a**2*x**2 + 1)) - 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**4*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1152():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -a*x**4*sqrt(-a**2*c*x**2 + c)/(4*sqrt(-a**2*x**2 + 1)) - x**3*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) - 2*x**2*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1)) - 4*x*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(-a**2*x**2 + 1)) - 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1153():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -a*x**3*sqrt(-a**2*c*x**2 + c)/(3*sqrt(-a**2*x**2 + 1)) - 3*x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1)) - 4*x*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1)) - 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1154():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -a*x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1)) - 3*x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) - 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/(a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1155():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x
    F = -a*x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) + sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1156():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**2
    F = 3*a*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*a*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1157():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**3
    F = 4*a**2*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**2*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - 3*a*sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*c*x**2 + c)/(2*x**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1158():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**4
    F = 4*a**3*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**3*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**2*sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1)) - 3*a*sqrt(-a**2*c*x**2 + c)/(2*x**2*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*c*x**2 + c)/(3*x**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1159():
    f = sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))/x**5
    F = 4*a**4*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**4*sqrt(-a**2*c*x**2 + c)*log(-a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**3*sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1)) - 2*a**2*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(-a**2*x**2 + 1)) - a*sqrt(-a**2*c*x**2 + c)/(x**3*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*c*x**2 + c)/(4*x**4*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1160():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(3*atanh(a*x))
    F = c*(a*x + 1)**4*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1161():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(3*atanh(a*x))
    F = -c**2*(a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(6*a*sqrt(-a**2*x**2 + 1)) + 2*c**2*(a*x + 1)**5*sqrt(-a**2*c*x**2 + c)/(5*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1162():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(3*atanh(a*x))
    F = c**3*(a*x + 1)**8*sqrt(-a**2*c*x**2 + c)/(8*a*sqrt(-a**2*x**2 + 1)) - 4*c**3*(a*x + 1)**7*sqrt(-a**2*c*x**2 + c)/(7*a*sqrt(-a**2*x**2 + 1)) + 2*c**3*(a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1163():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(3*atanh(a*x))
    F = -c**4*(a*x + 1)**10*sqrt(-a**2*c*x**2 + c)/(10*a*sqrt(-a**2*x**2 + 1)) + 2*c**4*(a*x + 1)**9*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1)) - 3*c**4*(a*x + 1)**8*sqrt(-a**2*c*x**2 + c)/(2*a*sqrt(-a**2*x**2 + 1)) + 8*c**4*(a*x + 1)**7*sqrt(-a**2*c*x**2 + c)/(7*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1164():
    f = exp(3*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(a*sqrt(-a**2*c*x**2 + c)) + 2*sqrt(-a**2*x**2 + 1)/(a*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1165():
    f = exp(3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = sqrt(-a**2*x**2 + 1)/(2*a*c*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1166():
    f = exp(3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = sqrt(-a**2*x**2 + 1)*atanh(a*x)/(8*a*c**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a*c**2*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(6*a*c**2*(-a*x + 1)**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1167():
    f = exp(3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = 5*sqrt(-a**2*x**2 + 1)*atanh(a*x)/(32*a*c**3*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(32*a*c**3*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a*c**3*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + 3*sqrt(-a**2*x**2 + 1)/(32*a*c**3*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(12*a*c**3*(-a*x + 1)**3*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(16*a*c**3*(-a*x + 1)**4*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1168():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(3*atanh(a*x))
    F = -a*x**(m + 2)*sqrt(-a**2*c*x**2 + c)/((m + 2)*sqrt(-a**2*x**2 + 1)) + 4*x**(m + 1)*sqrt(-a**2*c*x**2 + c)*hyper((1, m + 1), (m + 2,), a*x)/((m + 1)*sqrt(-a**2*x**2 + 1)) - 3*x**(m + 1)*sqrt(-a**2*c*x**2 + c)/((m + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1169():
    f = x**m*(-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))
    F = -a*x**(m + 2)*(-a**2*c*x**2 + c)**p/(sqrt(-a**2*x**2 + 1)*(m + 2*p + 1)) + a*x**(m + 2)*(-a**2*c*x**2 + c)**p*(4*m + 6*p + 5)*hyper((sympy.S(3)/2 - p, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*(-a**2*x**2 + 1)**p*(m + 2*p + 1)) - 3*x**(m + 1)*(-a**2*c*x**2 + c)**p/((m + 2*p)*sqrt(-a**2*x**2 + 1)) + x**(m + 1)*(-a**2*c*x**2 + c)**p*(4*m + 2*p + 3)*hyper((sympy.S(3)/2 - p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(m + 2*p)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1170():
    f = x**3*(-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))
    F = a*x**5*(6*p + 17)*(-a**2*c*x**2 + c)**p*hyper((sympy.S(5)/2, sympy.S(3)/2 - p), (sympy.S(7)/2,), a**2*x**2)/((10*p + 20)*(-a**2*x**2 + 1)**p) - a*x**5*(-a**2*c*x**2 + c)**p/((2*p + 4)*sqrt(-a**2*x**2 + 1)) - 3*(-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a**2*c*x**2 + c)**p/(a**4*(2*p + 3)) + 7*sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**4*(2*p + 1)) + 4*(-a**2*c*x**2 + c)**p/(a**4*(1 - 2*p)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1171():
    f = x**2*(-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))
    F = x**3*(2*p + 11)*(-a**2*c*x**2 + c)**p*hyper((sympy.S(3)/2, sympy.S(3)/2 - p), (sympy.S(5)/2,), a**2*x**2)/((6*p + 6)*(-a**2*x**2 + 1)**p) - 3*x**3*(-a**2*c*x**2 + c)**p/((2*p + 2)*sqrt(-a**2*x**2 + 1)) - (-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a**2*c*x**2 + c)**p/(a**3*(2*p + 3)) + 5*sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**3*(2*p + 1)) + 4*(-a**2*c*x**2 + c)**p/(a**3*(1 - 2*p)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1172():
    f = x*(-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))
    F = 3*2**(p + sympy.S(3)/2)*(-a*x + 1)**(p + sympy.S(-1)/2)*(-a**2*c*x**2 + c)**p*hyper((p + sympy.S(-1)/2, -p + sympy.S(-3)/2), (p + sympy.S.Half,), -a*x/2 + sympy.S.Half)/(a**2*(-a**2*x**2 + 1)**p*(-2*p**2 - p + 1)) - (a*x + 1)**3*(-a**2*c*x**2 + c)**p/(2*a**2*(p + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1173():
    f = (-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))
    F = 2**(p + sympy.S(5)/2)*(-a*x + 1)**(p + sympy.S(-1)/2)*(-a**2*c*x**2 + c)**p*hyper((p + sympy.S(-1)/2, -p + sympy.S(-3)/2), (p + sympy.S.Half,), -a*x/2 + sympy.S.Half)/(a*(1 - 2*p)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1174():
    f = (-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))/x
    F = a*x*(6*p + 1)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half, sympy.S(3)/2 - p), (sympy.S(3)/2,), a**2*x**2)/(2*p*(-a**2*x**2 + 1)**p) - a*x*(-a**2*c*x**2 + c)**p/(2*p*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) + 4*(-a**2*c*x**2 + c)**p/((1 - 2*p)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1175():
    f = (-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))/x**2
    F = a**2*x*(5 - 2*p)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half, sympy.S(3)/2 - p), (sympy.S(3)/2,), a**2*x**2)/(-a**2*x**2 + 1)**p - 3*a*sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) + 4*a*(-a**2*c*x**2 + c)**p/((1 - 2*p)*sqrt(-a**2*x**2 + 1)) - (-a**2*c*x**2 + c)**p/(x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1176():
    f = (-a**2*c*x**2 + c)**p*exp(3*atanh(a*x))/x**3
    F = a**3*x*(7 - 6*p)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half, sympy.S(3)/2 - p), (sympy.S(3)/2,), a**2*x**2)/(-a**2*x**2 + 1)**p + a**2*(9 - 2*p)*(-a**2*c*x**2 + c)**p*hyper((1, p + sympy.S(-1)/2), (p + sympy.S.Half,), -a**2*x**2 + 1)/((2 - 4*p)*sqrt(-a**2*x**2 + 1)) - 3*a*(-a**2*c*x**2 + c)**p/(x*sqrt(-a**2*x**2 + 1)) - (-a**2*c*x**2 + c)**p/(2*x**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1177():
    f = (-a**2*c*x**2 + c)**5*exp(4*atanh(a*x))
    F = -c**5*(a*x + 1)**11/(11*a) + 3*c**5*(a*x + 1)**10/(5*a) - 4*c**5*(a*x + 1)**9/(3*a) + c**5*(a*x + 1)**8/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1178():
    f = (-a**2*c*x**2 + c)**4*exp(4*atanh(a*x))
    F = c**4*(a*x + 1)**9/(9*a) - c**4*(a*x + 1)**8/(2*a) + 4*c**4*(a*x + 1)**7/(7*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1179():
    f = (-a**2*c*x**2 + c)**3*exp(4*atanh(a*x))
    F = -c**3*(a*x + 1)**7/(7*a) + c**3*(a*x + 1)**6/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1180():
    f = (-a**2*c*x**2 + c)**2*exp(4*atanh(a*x))
    F = c**2*(a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1181():
    f = exp(4*atanh(a*x))/(-a**2*c*x**2 + c)
    F = x/(c*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1182():
    f = exp(4*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = 1/(3*a*c**2*(-a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1183():
    f = exp(4*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = atanh(a*x)/(16*a*c**3) + 1/(16*a*c**3*(-a*x + 1)) + 1/(16*a*c**3*(-a*x + 1)**2) + 1/(12*a*c**3*(-a*x + 1)**3) + 1/(8*a*c**3*(-a*x + 1)**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1184():
    f = exp(4*atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = 3*atanh(a*x)/(32*a*c**4) - 1/(64*a*c**4*(a*x + 1)) + 5/(64*a*c**4*(-a*x + 1)) + 1/(16*a*c**4*(-a*x + 1)**2) + 1/(16*a*c**4*(-a*x + 1)**3) + 1/(16*a*c**4*(-a*x + 1)**4) + 1/(20*a*c**4*(-a*x + 1)**5)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1185():
    f = (-a**2*c*x**2 + c)**p*exp(4*atanh(a*x))
    F = 2**(p + 2)*c*(a*x + 1)**(1 - p)*(-a**2*c*x**2 + c)**(p - 1)*hyper((p - 1, -p - 2), (p,), -a*x/2 + sympy.S.Half)/(a*(1 - p))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1186():
    f = (-a**2*c*x**2 + c)**4*exp(-atanh(a*x))
    F = c**4*x*(-a**2*x**2 + 1)**(sympy.S(7)/2)/8 + 7*c**4*x*(-a**2*x**2 + 1)**(sympy.S(5)/2)/48 + 35*c**4*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/192 + 35*c**4*x*sqrt(-a**2*x**2 + 1)/128 + c**4*(-a**2*x**2 + 1)**(sympy.S(9)/2)/(9*a) + 35*c**4*asin(a*x)/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1187():
    f = (-a**2*c*x**2 + c)**3*exp(-atanh(a*x))
    F = c**3*x*(-a**2*x**2 + 1)**(sympy.S(5)/2)/6 + 5*c**3*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/24 + 5*c**3*x*sqrt(-a**2*x**2 + 1)/16 + c**3*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(7*a) + 5*c**3*asin(a*x)/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1188():
    f = (-a**2*c*x**2 + c)**2*exp(-atanh(a*x))
    F = c**2*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/4 + 3*c**2*x*sqrt(-a**2*x**2 + 1)/8 + c**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(5*a) + 3*c**2*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1189():
    f = (-a**2*c*x**2 + c)*exp(-atanh(a*x))
    F = c*x*sqrt(-a**2*x**2 + 1)/2 + c*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a) + c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1190():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)
    F = -exp(-atanh(a*x))/(a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1191():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = 2*x/(3*c**2*sqrt(-a**2*x**2 + 1)) - (-a*x + 1)/(3*a*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1192():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = 8*x/(15*c**3*sqrt(-a**2*x**2 + 1)) + 4*x/(15*c**3*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - (-a*x + 1)/(5*a*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1193():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = 16*x/(35*c**4*sqrt(-a**2*x**2 + 1)) + 8*x/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 6*x/(35*c**4*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - (-a*x + 1)/(7*a*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1194():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)**5
    F = 128*x/(315*c**5*sqrt(-a**2*x**2 + 1)) + 64*x/(315*c**5*(-a**2*x**2 + 1)**(sympy.S(3)/2)) + 16*x/(105*c**5*(-a**2*x**2 + 1)**(sympy.S(5)/2)) + 8*x/(63*c**5*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - (-a*x + 1)/(9*a*c**5*(-a**2*x**2 + 1)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1195():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(-atanh(a*x))
    F = -a*x**(m + 2)*sqrt(-a**2*c*x**2 + c)/((m + 2)*sqrt(-a**2*x**2 + 1)) + x**(m + 1)*sqrt(-a**2*c*x**2 + c)/((m + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1196():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(-atanh(a*x))
    F = -a*x**4*sqrt(-a**2*c*x**2 + c)/(4*sqrt(-a**2*x**2 + 1)) + x**3*sqrt(-a**2*c*x**2 + c)/(3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1197():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(-atanh(a*x))
    F = -a*x**3*sqrt(-a**2*c*x**2 + c)/(3*sqrt(-a**2*x**2 + 1)) + x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1198():
    f = sqrt(-a**2*c*x**2 + c)*exp(-atanh(a*x))
    F = -a*x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1)) + x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1199():
    f = sqrt(-a**2*c*x**2 + c)*exp(-atanh(a*x))/x
    F = -a*x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) + sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1200():
    f = sqrt(-a**2*c*x**2 + c)*exp(-atanh(a*x))/x**2
    F = -a*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1201():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-atanh(a*x))
    F = c*(-a*x + 1)**4*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(-a**2*x**2 + 1)) - 2*c*(-a*x + 1)**3*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1202():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(-atanh(a*x))
    F = -c**2*(-a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(6*a*sqrt(-a**2*x**2 + 1)) + 4*c**2*(-a*x + 1)**5*sqrt(-a**2*c*x**2 + c)/(5*a*sqrt(-a**2*x**2 + 1)) - c**2*(-a*x + 1)**4*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1203():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(-atanh(a*x))
    F = c**3*(-a*x + 1)**8*sqrt(-a**2*c*x**2 + c)/(8*a*sqrt(-a**2*x**2 + 1)) - 6*c**3*(-a*x + 1)**7*sqrt(-a**2*c*x**2 + c)/(7*a*sqrt(-a**2*x**2 + 1)) + 2*c**3*(-a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1)) - 8*c**3*(-a*x + 1)**5*sqrt(-a**2*c*x**2 + c)/(5*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1204():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(-atanh(a*x))
    F = -c**4*(-a*x + 1)**10*sqrt(-a**2*c*x**2 + c)/(10*a*sqrt(-a**2*x**2 + 1)) + 8*c**4*(-a*x + 1)**9*sqrt(-a**2*c*x**2 + c)/(9*a*sqrt(-a**2*x**2 + 1)) - 3*c**4*(-a*x + 1)**8*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1)) + 32*c**4*(-a*x + 1)**7*sqrt(-a**2*c*x**2 + c)/(7*a*sqrt(-a**2*x**2 + 1)) - 8*c**4*(-a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1205():
    f = exp(-atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(a*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1206():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = sqrt(-a**2*x**2 + 1)*atanh(a*x)/(2*a*c*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(2*a*c*(a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1207():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = 3*sqrt(-a**2*x**2 + 1)*atanh(a*x)/(8*a*c**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(4*a*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a*c**2*(a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a*c**2*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1208():
    f = exp(-atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = 5*sqrt(-a**2*x**2 + 1)*atanh(a*x)/(16*a*c**3*sqrt(-a**2*c*x**2 + c)) - 3*sqrt(-a**2*x**2 + 1)/(16*a*c**3*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - 3*sqrt(-a**2*x**2 + 1)/(32*a*c**3*(a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(24*a*c**3*(a*x + 1)**3*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(8*a*c**3*(-a*x + 1)*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(32*a*c**3*(-a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1209():
    f = x**m*(-a**2*c*x**2 + c)**p*exp(-atanh(a*x))
    F = -a*x**(m + 2)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half - p, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*(-a**2*x**2 + 1)**p) + x**(m + 1)*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half - p, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1210():
    f = x**3*(-a**2*x**2 + 1)**p*exp(-atanh(a*x))
    F = -a*x**5*hyper((sympy.S(5)/2, sympy.S.Half - p), (sympy.S(7)/2,), a**2*x**2)/5 + (-a**2*x**2 + 1)**(p + sympy.S(3)/2)/(a**4*(2*p + 3)) - (-a**2*x**2 + 1)**(p + sympy.S.Half)/(a**4*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1211():
    f = x**2*(-a**2*x**2 + 1)**p*exp(-atanh(a*x))
    F = x**3*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/3 - (-a**2*x**2 + 1)**(p + sympy.S(3)/2)/(a**3*(2*p + 3)) + (-a**2*x**2 + 1)**(p + sympy.S.Half)/(a**3*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1212():
    f = x*(-a**2*x**2 + 1)**p*exp(-atanh(a*x))
    F = -a*x**3*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/3 - (-a**2*x**2 + 1)**(p + sympy.S.Half)/(a**2*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1213():
    f = (-a**2*x**2 + 1)**p*exp(-atanh(a*x))
    F = -2**(p + sympy.S.Half)*(-a*x + 1)**(p + sympy.S(3)/2)*hyper((p + sympy.S(3)/2, sympy.S.Half - p), (p + sympy.S(5)/2,), -a*x/2 + sympy.S.Half)/(a*(2*p + 3))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1214():
    f = (-a**2*x**2 + 1)**p*exp(-atanh(a*x))/x
    F = -a*x*hyper((sympy.S.Half, sympy.S.Half - p), (sympy.S(3)/2,), a**2*x**2) - (-a**2*x**2 + 1)**(p + sympy.S.Half)*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1215():
    f = (-a**2*x**2 + 1)**p*exp(-atanh(a*x))/x**2
    F = a*(-a**2*x**2 + 1)**(p + sympy.S.Half)*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) - hyper((sympy.S(-1)/2, sympy.S.Half - p), (sympy.S.Half,), a**2*x**2)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1216():
    f = x**3*(-a**2*c*x**2 + c)**p*exp(-atanh(a*x))
    F = -a*x**5*(-a**2*c*x**2 + c)**p*hyper((sympy.S(5)/2, sympy.S.Half - p), (sympy.S(7)/2,), a**2*x**2)/(5*(-a**2*x**2 + 1)**p) + (-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a**2*c*x**2 + c)**p/(a**4*(2*p + 3)) - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**4*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1217():
    f = x**2*(-a**2*c*x**2 + c)**p*exp(-atanh(a*x))
    F = x**3*(-a**2*c*x**2 + c)**p*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/(3*(-a**2*x**2 + 1)**p) - (-a**2*x**2 + 1)**(sympy.S(3)/2)*(-a**2*c*x**2 + c)**p/(a**3*(2*p + 3)) + sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**3*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1218():
    f = x*(-a**2*c*x**2 + c)**p*exp(-atanh(a*x))
    F = -a*x**3*(-a**2*c*x**2 + c)**p*hyper((sympy.S(3)/2, sympy.S.Half - p), (sympy.S(5)/2,), a**2*x**2)/(3*(-a**2*x**2 + 1)**p) - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p/(a**2*(2*p + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1219():
    f = (-a**2*c*x**2 + c)**p*exp(-atanh(a*x))
    F = -2**(p + sympy.S.Half)*(-a*x + 1)**(p + sympy.S(3)/2)*(-a**2*c*x**2 + c)**p*hyper((p + sympy.S(3)/2, sympy.S.Half - p), (p + sympy.S(5)/2,), -a*x/2 + sympy.S.Half)/(a*(2*p + 3)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1220():
    f = (-a**2*c*x**2 + c)**p*exp(-atanh(a*x))/x
    F = -a*x*(-a**2*c*x**2 + c)**p*hyper((sympy.S.Half, sympy.S.Half - p), (sympy.S(3)/2,), a**2*x**2)/(-a**2*x**2 + 1)**p - sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1221():
    f = (-a**2*c*x**2 + c)**p*exp(-atanh(a*x))/x**2
    F = a*sqrt(-a**2*x**2 + 1)*(-a**2*c*x**2 + c)**p*hyper((1, p + sympy.S.Half), (p + sympy.S(3)/2,), -a**2*x**2 + 1)/(2*p + 1) - (-a**2*c*x**2 + c)**p*hyper((sympy.S(-1)/2, sympy.S.Half - p), (sympy.S.Half,), a**2*x**2)/(x*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1222():
    f = (-a**2*c*x**2 + c)**4*exp(-2*atanh(a*x))
    F = c**4*(-a*x + 1)**9/(9*a) - 3*c**4*(-a*x + 1)**8/(4*a) + 12*c**4*(-a*x + 1)**7/(7*a) - 4*c**4*(-a*x + 1)**6/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1223():
    f = (-a**2*c*x**2 + c)**3*exp(-2*atanh(a*x))
    F = -c**3*(-a*x + 1)**7/(7*a) + 2*c**3*(-a*x + 1)**6/(3*a) - 4*c**3*(-a*x + 1)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1224():
    f = (-a**2*c*x**2 + c)**2*exp(-2*atanh(a*x))
    F = c**2*(-a*x + 1)**5/(5*a) - c**2*(-a*x + 1)**4/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1225():
    f = (-a**2*c*x**2 + c)*exp(-2*atanh(a*x))
    F = -c*(-a*x + 1)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1226():
    f = exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)
    F = -1/(a*c*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1227():
    f = exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = atanh(a*x)/(4*a*c**2) - 1/(4*a*c**2*(a*x + 1)) - 1/(4*a*c**2*(a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1228():
    f = exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = atanh(a*x)/(4*a*c**3) - 3/(16*a*c**3*(a*x + 1)) - 1/(8*a*c**3*(a*x + 1)**2) - 1/(12*a*c**3*(a*x + 1)**3) + 1/(16*a*c**3*(-a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1229():
    f = exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = 15*atanh(a*x)/(64*a*c**4) - 5/(32*a*c**4*(a*x + 1)) - 3/(32*a*c**4*(a*x + 1)**2) - 1/(16*a*c**4*(a*x + 1)**3) - 1/(32*a*c**4*(a*x + 1)**4) + 5/(64*a*c**4*(-a*x + 1)) + 1/(64*a*c**4*(-a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1230():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))
    F = -x**4*sqrt(-a**2*c*x**2 + c)/5 + x**3*sqrt(-a**2*c*x**2 + c)/(2*a) - 3*x**2*sqrt(-a**2*c*x**2 + c)/(5*a**2) - 3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(4*a**4) - (-15*a*x + 24)*sqrt(-a**2*c*x**2 + c)/(20*a**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1231():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))
    F = -x**3*sqrt(-a**2*c*x**2 + c)/4 + 2*x**2*sqrt(-a**2*c*x**2 + c)/(3*a) + 7*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a**3) + (-21*a*x + 32)*sqrt(-a**2*c*x**2 + c)/(24*a**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1232():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))
    F = -x**2*sqrt(-a**2*c*x**2 + c)/3 - sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/a**2 - (-3*a*x + 5)*sqrt(-a**2*c*x**2 + c)/(3*a**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1233():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))
    F = 3*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(2*a) + (-a*x + 1)*sqrt(-a**2*c*x**2 + c)/(2*a) + 3*sqrt(-a**2*c*x**2 + c)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1234():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))/x
    F = -2*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) - sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - sqrt(-a**2*c*x**2 + c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1235():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))/x**2
    F = a*sqrt(c)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c)) + 2*a*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - sqrt(-a**2*c*x**2 + c)/x
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1236():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))/x**3
    F = -3*a**2*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/2 + 2*a*sqrt(-a**2*c*x**2 + c)/x - sqrt(-a**2*c*x**2 + c)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1237():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))/x**4
    F = a**3*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c)) - 5*a**2*sqrt(-a**2*c*x**2 + c)/(3*x) + a*sqrt(-a**2*c*x**2 + c)/x**2 - sqrt(-a**2*c*x**2 + c)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1238():
    f = sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))/x**5
    F = -7*a**4*sqrt(c)*atanh(sqrt(-a**2*c*x**2 + c)/sqrt(c))/8 + 4*a**3*sqrt(-a**2*c*x**2 + c)/(3*x) - 7*a**2*sqrt(-a**2*c*x**2 + c)/(8*x**2) + 2*a*sqrt(-a**2*c*x**2 + c)/(3*x**3) - sqrt(-a**2*c*x**2 + c)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1239():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-2*atanh(a*x))
    F = 5*c*x*sqrt(-a**2*c*x**2 + c)/8 + 5*c**(sympy.S(3)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(8*a) + (-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(4*a) + 5*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/(12*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1240():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(-2*atanh(a*x))
    F = 7*c**2*x*sqrt(-a**2*c*x**2 + c)/16 + 7*c*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/24 + 7*c**(sympy.S(5)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(16*a) + (-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(6*a) + 7*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/(30*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1241():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(-2*atanh(a*x))
    F = 45*c**3*x*sqrt(-a**2*c*x**2 + c)/128 + 15*c**2*x*(-a**2*c*x**2 + c)**(sympy.S(3)/2)/64 + 3*c*x*(-a**2*c*x**2 + c)**(sympy.S(5)/2)/16 + 45*c**(sympy.S(7)/2)*atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(128*a) + (-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(8*a) + 9*(-a**2*c*x**2 + c)**(sympy.S(7)/2)/(56*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1242():
    f = exp(-2*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -(-2*a*x + 2)/(a*sqrt(-a**2*c*x**2 + c)) - atan(a*sqrt(c)*x/sqrt(-a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1243():
    f = exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x/(3*c*sqrt(-a**2*c*x**2 + c)) - (-2*a*x + 2)/(3*a*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1244():
    f = exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = x/(5*c*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + 2*x/(5*c**2*sqrt(-a**2*c*x**2 + c)) - (-2*a*x + 2)/(5*a*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1245():
    f = exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = x/(7*c*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) + 4*x/(21*c**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + 8*x/(21*c**3*sqrt(-a**2*c*x**2 + c)) - (-2*a*x + 2)/(7*a*(-a**2*c*x**2 + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1246():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(-2*atanh(a*x))
    F = -2*a*c*x**(m + 2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), a**2*x**2)/((m + 2)*sqrt(-a**2*c*x**2 + c)) + c*x**(m + 1)*(2*m + 3)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), a**2*x**2)/((m + 1)*(m + 2)*sqrt(-a**2*c*x**2 + c)) - x**(m + 1)*sqrt(-a**2*c*x**2 + c)/(m + 2)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1247():
    f = (-a**2*c*x**2 + c)**p*exp(-2*atanh(a*x))
    F = 2**(p + 1)*(-a**2*c*x**2 + c)**p*hyper((p, -p - 1), (p + 1,), a*x/2 + sympy.S.Half)/(a*p*(-a*x + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1248():
    f = (-a**2*c*x**2 + c)**4*exp(-3*atanh(a*x))
    F = 11*c**4*x*(-a**2*x**2 + 1)**(sympy.S(5)/2)/48 + 55*c**4*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/192 + 55*c**4*x*sqrt(-a**2*x**2 + 1)/128 + c**4*(-a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(9*a) + 11*c**4*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(72*a) + 11*c**4*(-a**2*x**2 + 1)**(sympy.S(7)/2)/(56*a) + 55*c**4*asin(a*x)/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1249():
    f = (-a**2*c*x**2 + c)**3*exp(-3*atanh(a*x))
    F = 3*c**3*x*(-a**2*x**2 + 1)**(sympy.S(3)/2)/8 + 9*c**3*x*sqrt(-a**2*x**2 + 1)/16 + c**3*(-a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(7*a) + 3*c**3*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(14*a) + 3*c**3*(-a**2*x**2 + 1)**(sympy.S(5)/2)/(10*a) + 9*c**3*asin(a*x)/(16*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1250():
    f = (-a**2*c*x**2 + c)**2*exp(-3*atanh(a*x))
    F = 7*c**2*x*sqrt(-a**2*x**2 + 1)/8 + c**2*(-a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(5*a) + 7*c**2*(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(20*a) + 7*c**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)/(12*a) + 7*c**2*asin(a*x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1251():
    f = (-a**2*c*x**2 + c)*exp(-3*atanh(a*x))
    F = c*(-a*x + 1)**2*sqrt(-a**2*x**2 + 1)/(3*a) + 5*c*(-a*x + 1)*sqrt(-a**2*x**2 + 1)/(6*a) + 5*c*sqrt(-a**2*x**2 + 1)/(2*a) + 5*c*asin(a*x)/(2*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1252():
    f = exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)
    F = -exp(-3*atanh(a*x))/(3*a*c)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1253():
    f = exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -2*sqrt(-a**2*x**2 + 1)/(15*a*c**2*(a*x + 1)) - 2*sqrt(-a**2*x**2 + 1)/(15*a*c**2*(a*x + 1)**2) - sqrt(-a**2*x**2 + 1)/(5*a*c**2*(a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1254():
    f = exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = 8*x/(35*c**3*sqrt(-a**2*x**2 + 1)) - 4/(35*a*c**3*(a*x + 1)*sqrt(-a**2*x**2 + 1)) - 4/(35*a*c**3*(a*x + 1)**2*sqrt(-a**2*x**2 + 1)) - 1/(7*a*c**3*(a*x + 1)**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1255():
    f = exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = 16*x/(63*c**4*sqrt(-a**2*x**2 + 1)) + 8*x/(63*c**4*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 2/(21*a*c**4*(a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 2/(21*a*c**4*(a*x + 1)**2*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 1/(9*a*c**4*(a*x + 1)**3*(-a**2*x**2 + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1256():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))
    F = a*x**5*sqrt(-a**2*c*x**2 + c)/(5*sqrt(-a**2*x**2 + 1)) - 3*x**4*sqrt(-a**2*c*x**2 + c)/(4*sqrt(-a**2*x**2 + 1)) + 4*x**3*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1)) - 2*x**2*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(-a**2*x**2 + 1)) + 4*x*sqrt(-a**2*c*x**2 + c)/(a**3*sqrt(-a**2*x**2 + 1)) - 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**4*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1257():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))
    F = a*x**4*sqrt(-a**2*c*x**2 + c)/(4*sqrt(-a**2*x**2 + 1)) - x**3*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) + 2*x**2*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1)) - 4*x*sqrt(-a**2*c*x**2 + c)/(a**2*sqrt(-a**2*x**2 + 1)) + 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1258():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))
    F = a*x**3*sqrt(-a**2*c*x**2 + c)/(3*sqrt(-a**2*x**2 + 1)) - 3*x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1)) + 4*x*sqrt(-a**2*c*x**2 + c)/(a*sqrt(-a**2*x**2 + 1)) - 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1259():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))
    F = a*x**2*sqrt(-a**2*c*x**2 + c)/(2*sqrt(-a**2*x**2 + 1)) - 3*x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) + 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/(a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1260():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))/x
    F = a*x*sqrt(-a**2*c*x**2 + c)/sqrt(-a**2*x**2 + 1) + sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/sqrt(-a**2*x**2 + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1261():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))/x**2
    F = -3*a*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) + 4*a*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/sqrt(-a**2*x**2 + 1) - sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1262():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))/x**3
    F = 4*a**2*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**2*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/sqrt(-a**2*x**2 + 1) + 3*a*sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*c*x**2 + c)/(2*x**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1263():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))/x**4
    F = -4*a**3*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) + 4*a**3*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/sqrt(-a**2*x**2 + 1) - 4*a**2*sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1)) + 3*a*sqrt(-a**2*c*x**2 + c)/(2*x**2*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*c*x**2 + c)/(3*x**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1264():
    f = sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))/x**5
    F = 4*a**4*sqrt(-a**2*c*x**2 + c)*log(x)/sqrt(-a**2*x**2 + 1) - 4*a**4*sqrt(-a**2*c*x**2 + c)*log(a*x + 1)/sqrt(-a**2*x**2 + 1) + 4*a**3*sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1)) - 2*a**2*sqrt(-a**2*c*x**2 + c)/(x**2*sqrt(-a**2*x**2 + 1)) + a*sqrt(-a**2*c*x**2 + c)/(x**3*sqrt(-a**2*x**2 + 1)) - sqrt(-a**2*c*x**2 + c)/(4*x**4*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1265():
    f = (-a**2*c*x**2 + c)**(sympy.S(9)/2)*exp(-3*atanh(a*x))
    F = c**4*(-a*x + 1)**10*sqrt(-a**2*c*x**2 + c)/(10*a*sqrt(-a**2*x**2 + 1)) - 2*c**4*(-a*x + 1)**9*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1)) + 3*c**4*(-a*x + 1)**8*sqrt(-a**2*c*x**2 + c)/(2*a*sqrt(-a**2*x**2 + 1)) - 8*c**4*(-a*x + 1)**7*sqrt(-a**2*c*x**2 + c)/(7*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1266():
    f = (-a**2*c*x**2 + c)**(sympy.S(7)/2)*exp(-3*atanh(a*x))
    F = -c**3*(-a*x + 1)**8*sqrt(-a**2*c*x**2 + c)/(8*a*sqrt(-a**2*x**2 + 1)) + 4*c**3*(-a*x + 1)**7*sqrt(-a**2*c*x**2 + c)/(7*a*sqrt(-a**2*x**2 + 1)) - 2*c**3*(-a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(3*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1267():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(-3*atanh(a*x))
    F = c**2*(-a*x + 1)**6*sqrt(-a**2*c*x**2 + c)/(6*a*sqrt(-a**2*x**2 + 1)) - 2*c**2*(-a*x + 1)**5*sqrt(-a**2*c*x**2 + c)/(5*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1268():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-3*atanh(a*x))
    F = -c*(-a*x + 1)**4*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1269():
    f = exp(-3*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(a*sqrt(-a**2*c*x**2 + c)) - 2*sqrt(-a**2*x**2 + 1)/(a*(a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1270():
    f = exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -sqrt(-a**2*x**2 + 1)/(2*a*c*(a*x + 1)**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1271():
    f = exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = sqrt(-a**2*x**2 + 1)*atanh(a*x)/(8*a*c**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a*c**2*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a*c**2*(a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(6*a*c**2*(a*x + 1)**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1272():
    f = exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = 5*sqrt(-a**2*x**2 + 1)*atanh(a*x)/(32*a*c**3*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(8*a*c**3*(a*x + 1)*sqrt(-a**2*c*x**2 + c)) - 3*sqrt(-a**2*x**2 + 1)/(32*a*c**3*(a*x + 1)**2*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(12*a*c**3*(a*x + 1)**3*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(16*a*c**3*(a*x + 1)**4*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(32*a*c**3*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1273():
    f = x**m*sqrt(-a**2*c*x**2 + c)*exp(-3*atanh(a*x))
    F = a*x**(m + 2)*sqrt(-a**2*c*x**2 + c)/((m + 2)*sqrt(-a**2*x**2 + 1)) + 4*x**(m + 1)*sqrt(-a**2*c*x**2 + c)*hyper((1, m + 1), (m + 2,), -a*x)/((m + 1)*sqrt(-a**2*x**2 + 1)) - 3*x**(m + 1)*sqrt(-a**2*c*x**2 + c)/((m + 1)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1274():
    f = (-a**2*c*x**2 + c)**p*exp(-3*atanh(a*x))
    F = -2**(p + sympy.S(-1)/2)*(-a*x + 1)**(p + sympy.S(5)/2)*(-a**2*c*x**2 + c)**p*hyper((p + sympy.S(5)/2, sympy.S(3)/2 - p), (p + sympy.S(7)/2,), -a*x/2 + sympy.S.Half)/(a*(2*p + 5)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1275():
    f = (-a**2*x**2 + 1)**(sympy.S(5)/2)*exp(atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(13)/4)*(a*x + 1)**(sympy.S(11)/4)/(6*a) - 11*(-a*x + 1)**(sympy.S(13)/4)*(a*x + 1)**(sympy.S(7)/4)/(60*a) - 77*(-a*x + 1)**(sympy.S(13)/4)*(a*x + 1)**(sympy.S(3)/4)/(480*a) + 77*(-a*x + 1)**(sympy.S(9)/4)*(a*x + 1)**(sympy.S(3)/4)/(960*a) + 231*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(1280*a) + 231*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(512*a) + 231*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2048*a) - 231*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2048*a) - 231*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(1024*a) - 231*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(1024*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1276():
    f = (-a**2*x**2 + 1)**(sympy.S(3)/2)*exp(atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(9)/4)*(a*x + 1)**(sympy.S(7)/4)/(4*a) - 7*(-a*x + 1)**(sympy.S(9)/4)*(a*x + 1)**(sympy.S(3)/4)/(24*a) + 7*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(32*a) + 35*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(64*a) + 35*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a) - 35*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a) - 35*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(128*a) - 35*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(128*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1277():
    f = sqrt(-a**2*x**2 + 1)*exp(atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)/(2*a) + 3*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)/(4*a) + 3*sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a) - 3*sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a) - 3*sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1278():
    f = exp(atanh(a*x)/2)/sqrt(-a**2*x**2 + 1)
    F = sqrt(2)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2*a) - sqrt(2)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2*a) - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/a - sqrt(2)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/a
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1279():
    f = exp(atanh(a*x)/2)/(-a**2*x**2 + 1)**(sympy.S(3)/2)
    F = -2*(-2*a*x + 1)*exp(atanh(a*x)/2)/(3*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1280():
    f = exp(atanh(a*x)/2)/(-a**2*x**2 + 1)**(sympy.S(5)/2)
    F = -2*(-6*a*x + 1)*exp(atanh(a*x)/2)/(35*a*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 16*(-2*a*x + 1)*exp(atanh(a*x)/2)/(35*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1281():
    f = exp(atanh(a*x)/2)/(-a**2*x**2 + 1)**(sympy.S(7)/2)
    F = -2*(-10*a*x + 1)*exp(atanh(a*x)/2)/(99*a*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 32*(-6*a*x + 1)*exp(atanh(a*x)/2)/(693*a*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 256*(-2*a*x + 1)*exp(atanh(a*x)/2)/(693*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1282():
    f = exp(atanh(a*x)/2)/(-a**2*x**2 + 1)**(sympy.S(9)/2)
    F = -2*(-14*a*x + 1)*exp(atanh(a*x)/2)/(195*a*(-a**2*x**2 + 1)**(sympy.S(7)/2)) - 112*(-10*a*x + 1)*exp(atanh(a*x)/2)/(6435*a*(-a**2*x**2 + 1)**(sympy.S(5)/2)) - 256*(-6*a*x + 1)*exp(atanh(a*x)/2)/(6435*a*(-a**2*x**2 + 1)**(sympy.S(3)/2)) - 2048*(-2*a*x + 1)*exp(atanh(a*x)/2)/(6435*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1283():
    f = (-a**2*c*x**2 + c)**(sympy.S(5)/2)*exp(atanh(a*x)/2)
    F = -c**2*(-a*x + 1)**(sympy.S(13)/4)*(a*x + 1)**(sympy.S(11)/4)*sqrt(-a**2*c*x**2 + c)/(6*a*sqrt(-a**2*x**2 + 1)) - 11*c**2*(-a*x + 1)**(sympy.S(13)/4)*(a*x + 1)**(sympy.S(7)/4)*sqrt(-a**2*c*x**2 + c)/(60*a*sqrt(-a**2*x**2 + 1)) - 77*c**2*(-a*x + 1)**(sympy.S(13)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(480*a*sqrt(-a**2*x**2 + 1)) + 77*c**2*(-a*x + 1)**(sympy.S(9)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(960*a*sqrt(-a**2*x**2 + 1)) + 231*c**2*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(1280*a*sqrt(-a**2*x**2 + 1)) + 231*c**2*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(512*a*sqrt(-a**2*x**2 + 1)) + 231*sqrt(2)*c**2*sqrt(-a**2*c*x**2 + c)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2048*a*sqrt(-a**2*x**2 + 1)) - 231*sqrt(2)*c**2*sqrt(-a**2*c*x**2 + c)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2048*a*sqrt(-a**2*x**2 + 1)) - 231*sqrt(2)*c**2*sqrt(-a**2*c*x**2 + c)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(1024*a*sqrt(-a**2*x**2 + 1)) - 231*sqrt(2)*c**2*sqrt(-a**2*c*x**2 + c)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(1024*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1284():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(atanh(a*x)/2)
    F = -c*(-a*x + 1)**(sympy.S(9)/4)*(a*x + 1)**(sympy.S(7)/4)*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(-a**2*x**2 + 1)) - 7*c*(-a*x + 1)**(sympy.S(9)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(24*a*sqrt(-a**2*x**2 + 1)) + 7*c*(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(32*a*sqrt(-a**2*x**2 + 1)) + 35*c*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(64*a*sqrt(-a**2*x**2 + 1)) + 35*sqrt(2)*c*sqrt(-a**2*c*x**2 + c)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a*sqrt(-a**2*x**2 + 1)) - 35*sqrt(2)*c*sqrt(-a**2*c*x**2 + c)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(256*a*sqrt(-a**2*x**2 + 1)) - 35*sqrt(2)*c*sqrt(-a**2*c*x**2 + c)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(128*a*sqrt(-a**2*x**2 + 1)) - 35*sqrt(2)*c*sqrt(-a**2*c*x**2 + c)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(128*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1285():
    f = sqrt(-a**2*c*x**2 + c)*exp(atanh(a*x)/2)
    F = -(-a*x + 1)**(sympy.S(5)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(2*a*sqrt(-a**2*x**2 + 1)) + 3*(-a*x + 1)**(sympy.S(1)/4)*(a*x + 1)**(sympy.S(3)/4)*sqrt(-a**2*c*x**2 + c)/(4*a*sqrt(-a**2*x**2 + 1)) + 3*sqrt(2)*sqrt(-a**2*c*x**2 + c)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a*sqrt(-a**2*x**2 + 1)) - 3*sqrt(2)*sqrt(-a**2*c*x**2 + c)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(16*a*sqrt(-a**2*x**2 + 1)) - 3*sqrt(2)*sqrt(-a**2*c*x**2 + c)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(8*a*sqrt(-a**2*x**2 + 1)) - 3*sqrt(2)*sqrt(-a**2*c*x**2 + c)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(8*a*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1286():
    f = exp(atanh(a*x)/2)/sqrt(-a**2*c*x**2 + c)
    F = sqrt(2)*sqrt(-a**2*x**2 + 1)*log(-sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2*a*sqrt(-a**2*c*x**2 + c)) - sqrt(2)*sqrt(-a**2*x**2 + 1)*log(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + sqrt(-a*x + 1)/sqrt(a*x + 1) + 1)/(2*a*sqrt(-a**2*c*x**2 + c)) - sqrt(2)*sqrt(-a**2*x**2 + 1)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) - 1)/(a*sqrt(-a**2*c*x**2 + c)) - sqrt(2)*sqrt(-a**2*x**2 + 1)*atan(sqrt(2)*(-a*x + 1)**(sympy.S(1)/4)/(a*x + 1)**(sympy.S(1)/4) + 1)/(a*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1287():
    f = exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -2*(-2*a*x + 1)*exp(atanh(a*x)/2)/(3*a*c*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1288():
    f = exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -2*(-6*a*x + 1)*exp(atanh(a*x)/2)/(35*a*c*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 16*(-2*a*x + 1)*exp(atanh(a*x)/2)/(35*a*c**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1289():
    f = exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -2*(-10*a*x + 1)*exp(atanh(a*x)/2)/(99*a*c*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 32*(-6*a*x + 1)*exp(atanh(a*x)/2)/(693*a*c**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 256*(-2*a*x + 1)*exp(atanh(a*x)/2)/(693*a*c**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1290():
    f = exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(9)/2)
    F = -2*(-14*a*x + 1)*exp(atanh(a*x)/2)/(195*a*c*(-a**2*c*x**2 + c)**(sympy.S(7)/2)) - 112*(-10*a*x + 1)*exp(atanh(a*x)/2)/(6435*a*c**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 256*(-6*a*x + 1)*exp(atanh(a*x)/2)/(6435*a*c**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 2048*(-2*a*x + 1)*exp(atanh(a*x)/2)/(6435*a*c**4*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1291():
    f = x**3*exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(5)/4)
    F = -2*(-a*x + 1)**(sympy.S(3)/2)*(-a**2*x**2 + 1)**(sympy.S(1)/4)/(3*a**4*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + 2*sqrt(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(1)/4)/(a**4*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + sqrt(2)*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(2)*sqrt(-a*x + 1)/2)/(2*a**4*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + (-a**2*x**2 + 1)**(sympy.S(1)/4)/(a**4*c*sqrt(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1292():
    f = x**2*exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(5)/4)
    F = 2*sqrt(-a*x + 1)*(-a**2*x**2 + 1)**(sympy.S(1)/4)/(a**3*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) - sqrt(2)*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(2)*sqrt(-a*x + 1)/2)/(2*a**3*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + (-a**2*x**2 + 1)**(sympy.S(1)/4)/(a**3*c*sqrt(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1293():
    f = x*exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(5)/4)
    F = sqrt(2)*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(2)*sqrt(-a*x + 1)/2)/(2*a**2*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + (-a**2*x**2 + 1)**(sympy.S(1)/4)/(a**2*c*sqrt(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1294():
    f = exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(5)/4)
    F = -sqrt(2)*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(2)*sqrt(-a*x + 1)/2)/(2*a*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + (-a**2*x**2 + 1)**(sympy.S(1)/4)/(a*c*sqrt(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1295():
    f = exp(atanh(a*x)/2)/(x*(-a**2*c*x**2 + c)**(sympy.S(5)/4))
    F = sqrt(2)*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(2)*sqrt(-a*x + 1)/2)/(2*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) - 2*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(-a*x + 1))/(c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + (-a**2*x**2 + 1)**(sympy.S(1)/4)/(c*sqrt(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1296():
    f = exp(atanh(a*x)/2)/(x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/4))
    F = -sqrt(2)*a*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(2)*sqrt(-a*x + 1)/2)/(2*c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) - a*(-a**2*x**2 + 1)**(sympy.S(1)/4)*atanh(sqrt(-a*x + 1))/(c*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) + 2*a*(-a**2*x**2 + 1)**(sympy.S(1)/4)/(c*sqrt(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(1)/4)) - (-a**2*x**2 + 1)**(sympy.S(1)/4)/(c*x*sqrt(-a*x + 1)*(-a**2*c*x**2 + c)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1297():
    f = x**3*exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(9)/8)
    F = -4*x**2*(a*x + 1)**(sympy.S(1)/8)*(-a**2*x**2 + 1)**(sympy.S(1)/8)/(7*a**2*c*(-a*x + 1)**(sympy.S(3)/8)*(-a**2*c*x**2 + c)**(sympy.S(1)/8)) + (-8*a*x + 48)*(a*x + 1)**(sympy.S(1)/8)*(-a**2*x**2 + 1)**(sympy.S(1)/8)/(21*a**4*c*(-a*x + 1)**(sympy.S(3)/8)*(-a**2*c*x**2 + c)**(sympy.S(1)/8)) + 64*2**(sympy.S(1)/8)*(-a*x + 1)**(sympy.S(5)/8)*(-a**2*x**2 + 1)**(sympy.S(1)/8)*hyper((sympy.S(5)/8, sympy.S(7)/8), (sympy.S(13)/8,), -a*x/2 + sympy.S.Half)/(105*a**4*c*(-a**2*c*x**2 + c)**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1298():
    f = x**2*exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(9)/8)
    F = 4*(-a*x + 2)*exp(atanh(a*x)/2)/(3*a**3*c*(-a**2*c*x**2 + c)**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1299():
    f = x*exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(9)/8)
    F = 8*2**(sympy.S(1)/8)*(-a*x + 1)**(sympy.S(5)/8)*(-a**2*x**2 + 1)**(sympy.S(1)/8)*hyper((sympy.S(5)/8, sympy.S(7)/8), (sympy.S(13)/8,), -a*x/2 + sympy.S.Half)/(15*a**2*c*(-a**2*c*x**2 + c)**(sympy.S(1)/8)) + 4*(a*x + 1)**(sympy.S(1)/8)*(-a**2*x**2 + 1)**(sympy.S(1)/8)/(3*a**2*c*(-a*x + 1)**(sympy.S(3)/8)*(-a**2*c*x**2 + c)**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1300():
    f = exp(atanh(a*x)/2)/(-a**2*c*x**2 + c)**(sympy.S(9)/8)
    F = 4*2**(sympy.S(1)/8)*(-a**2*x**2 + 1)**(sympy.S(1)/8)*hyper((sympy.S(-3)/8, sympy.S(7)/8), (sympy.S(5)/8,), -a*x/2 + sympy.S.Half)/(3*a*c*(-a*x + 1)**(sympy.S(3)/8)*(-a**2*c*x**2 + c)**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1301():
    f = exp(atanh(a*x)/2)/(x*(-a**2*c*x**2 + c)**(sympy.S(9)/8))
    F = -2*2**(sympy.S(5)/8)*(a*x + 1)**(sympy.S(1)/8)*(-a**2*x**2 + 1)**(sympy.S(1)/8)*appellf1(sympy.S(1)/8, 1, sympy.S(11)/8, sympy.S(9)/8, a*x + 1, a*x/2 + sympy.S.Half)/(c*(-a**2*c*x**2 + c)**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1302():
    f = (-a**2*c*x**2 + c)*exp(n*atanh(a*x))
    F = -2**(n/2 + 2)*c*(-a*x + 1)**(2 - n/2)*hyper((2 - n/2, -n/2 - 1), (3 - n/2,), -a*x/2 + sympy.S.Half)/(a*(4 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1303():
    f = (-a**2*c*x**2 + c)**2*exp(n*atanh(a*x))
    F = -2**(n/2 + 3)*c**2*(-a*x + 1)**(3 - n/2)*hyper((3 - n/2, -n/2 - 2), (4 - n/2,), -a*x/2 + sympy.S.Half)/(a*(6 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1304():
    f = (-a**2*c*x**2 + c)**3*exp(n*atanh(a*x))
    F = -2**(n/2 + 4)*c**3*(-a*x + 1)**(4 - n/2)*hyper((4 - n/2, -n/2 - 3), (5 - n/2,), -a*x/2 + sympy.S.Half)/(a*(8 - n))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1305():
    f = x**4*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)
    F = 2**(n/2 - 1)*n*(n**2 + 8)*(-a*x + 1)**(1 - n/2)*hyper((1 - n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(3*a**5*c*(2 - n)) - x**3*(a*x + 1)**(n/2)/(3*a**2*c*(-a*x + 1)**(n/2)) - n*x**2*(a*x + 1)**(n/2)/(6*a**3*c*(-a*x + 1)**(n/2)) + (a*x + 1)**(n/2)*(-a*n*x*(n**2 + 6) + n**3 + n**2 + 8*n + 6)/(6*a**5*c*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1306():
    f = x**3*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)
    F = 2**(n/2 - 1)*(n**2 + 2)*(-a*x + 1)**(1 - n/2)*hyper((1 - n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(a**4*c*(2 - n)) - x**2*(a*x + 1)**(n/2)/(2*a**2*c*(-a*x + 1)**(n/2)) + (a*x + 1)**(n/2)*(-a*n**2*x + n**2 + n + 2)/(2*a**4*c*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1307():
    f = x*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)
    F = 2**(n/2 + 1)*hyper((-n/2, -n/2), (1 - n/2,), -a*x/2 + sympy.S.Half)/(a**2*c*n*(-a*x + 1)**(n/2)) - (a*x + 1)**(n/2)/(a**2*c*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1308():
    f = exp(n*atanh(a*x))/(-a**2*c*x**2 + c)
    F = exp(n*atanh(a*x))/(a*c*n)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1309():
    f = exp(n*atanh(a*x))/(x*(-a**2*c*x**2 + c))
    F = -2*(a*x + 1)**(n/2)*hyper((1, n/2), (n/2 + 1,), (a*x + 1)/(-a*x + 1))/(c*n*(-a*x + 1)**(n/2)) + (a*x + 1)**(n/2)/(c*n*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1310():
    f = exp(n*atanh(a*x))/(x**2*(-a**2*c*x**2 + c))
    F = -2*a*(a*x + 1)**(n/2)*hyper((1, n/2), (n/2 + 1,), (a*x + 1)/(-a*x + 1))/(c*(-a*x + 1)**(n/2)) + a*(n + 1)*(a*x + 1)**(n/2)/(c*n*(-a*x + 1)**(n/2)) - (a*x + 1)**(n/2)/(c*x*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1311():
    f = x**4*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -2**(n/2)*n*(-a*x + 1)**(1 - n/2)*hyper((1 - n/2, 1 - n/2), (2 - n/2,), -a*x/2 + sympy.S.Half)/(a**5*c**2*(2 - n)) - x**3*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(a**2*c**2) + x*(n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(a**4*c**2) + (1 - n)*(n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(a**5*c**2*(2 - n)) - (2 - n**2)*(n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2)/(a**5*c**2*(4 - n**2)) - (a*x + 1)**(n/2 - 1)/(a**5*c**2*(-a*x + 1)**(n/2)) + (-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)/(a**5*c**2*(2 - n)) - (2 - n**2)*(n + 3)*(a*x + 1)**(n/2)/(a**5*c**2*n*(4 - n**2)*(-a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1312():
    f = x**3*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = 2**(n/2 + 2)*(-a*x + 1)**(-n/2 - 1)*hyper((-n/2 - 1, -n/2 - 1), (-n/2,), -a*x/2 + sympy.S.Half)/(a**4*c**2*(n + 2)) + 3*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2)/(a**4*c**2*(n + 2)) - (-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(a**4*c**2*(n + 2)) - 3*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 + 1)/(a**4*c**2*(n + 2)) + 3*(a*x + 1)**(n/2)/(a**4*c**2*n*(n + 2)*(-a*x + 1)**(n/2)) - 2*(a*x + 1)**(n/2 - 1)/(a**4*c**2*n*(n + 2)*(-a*x + 1)**(n/2)) + 2*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)/(a**4*c**2*n*(4 - n**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1313():
    f = x**2*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -(-2*a*x + n)*exp(n*atanh(a*x))/(a**3*c**2*(4 - n**2)*(-a**2*x**2 + 1)) - (2 - n**2)*exp(n*atanh(a*x))/(a**3*c**2*n*(4 - n**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1314():
    f = exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = -(-2*a*x + n)*exp(n*atanh(a*x))/(a*c**2*(4 - n**2)*(-a**2*x**2 + 1)) + 2*exp(n*atanh(a*x))/(a*c**2*n*(4 - n**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1315():
    f = exp(n*atanh(a*x))/(x*(-a**2*c*x**2 + c)**2)
    F = (-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(c**2*(n + 2)) - 2*(a*x + 1)**(n/2)*hyper((1, n/2), (n/2 + 1,), (a*x + 1)/(-a*x + 1))/(c**2*n*(-a*x + 1)**(n/2)) + (n + 4)*(a*x + 1)**(n/2 - 1)/(c**2*n*(n + 2)*(-a*x + 1)**(n/2)) - (-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)*(-n**2 - n + 4)/(c**2*n*(4 - n**2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1316():
    f = exp(n*atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**2)
    F = -2*a*(a*x + 1)**(n/2)*hyper((1, n/2), (n/2 + 1,), (a*x + 1)/(-a*x + 1))/(c**2*(-a*x + 1)**(n/2)) + a*(n + 3)*(-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(c**2*(n + 2)) + a*(a*x + 1)**(n/2 - 1)*(n**2 + 4*n + 6)/(c**2*n*(n + 2)*(-a*x + 1)**(n/2)) - a*(-a*x + 1)**(1 - n/2)*(a*x + 1)**(n/2 - 1)*(-n**3 - n**2 + 4*n + 6)/(c**2*n*(4 - n**2)) - (-a*x + 1)**(-n/2 - 1)*(a*x + 1)**(n/2 - 1)/(c**2*x)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1317():
    f = exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = -(-4*a*x + n)*exp(n*atanh(a*x))/(a*c**3*(16 - n**2)*(-a**2*x**2 + 1)**2) - 12*(-2*a*x + n)*exp(n*atanh(a*x))/(a*c**3*(4 - n**2)*(16 - n**2)*(-a**2*x**2 + 1)) + 24*exp(n*atanh(a*x))/(a*c**3*n*(n**4 - 20*n**2 + 64))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1318():
    f = exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**4
    F = -(-6*a*x + n)*exp(n*atanh(a*x))/(a*c**4*(36 - n**2)*(-a**2*x**2 + 1)**3) - 30*(-4*a*x + n)*exp(n*atanh(a*x))/(a*c**4*(16 - n**2)*(36 - n**2)*(-a**2*x**2 + 1)**2) - 360*(-2*a*x + n)*exp(n*atanh(a*x))/(a*c**4*(4 - n**2)*(16 - n**2)*(36 - n**2)*(-a**2*x**2 + 1)) + 720*exp(n*atanh(a*x))/(a*c**4*n*(36 - n**2)*(n**4 - 20*n**2 + 64))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1319():
    f = x**3*sqrt(-a**2*c*x**2 + c)*exp(n*atanh(a*x))
    F = -2**(n/2 + sympy.S(-1)/2)*n*(n**2 + 11)*(-a*x + 1)**(sympy.S(3)/2 - n/2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(3)/2 - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(15*a**4*(3 - n)*sqrt(-a**2*x**2 + 1)) - x**2*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(3)/2)*sqrt(-a**2*c*x**2 + c)/(5*a**2*sqrt(-a**2*x**2 + 1)) - (-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(3)/2)*sqrt(-a**2*c*x**2 + c)*(3*a*n*x + n**2 + 8)/(60*a**4*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1320():
    f = x**2*sqrt(-a**2*c*x**2 + c)*exp(n*atanh(a*x))
    F = -2**(n/2 + sympy.S(-1)/2)*(n**2 + 3)*(-a*x + 1)**(sympy.S(3)/2 - n/2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(3)/2 - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(3*a**3*(3 - n)*sqrt(-a**2*x**2 + 1)) - x*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(3)/2)*sqrt(-a**2*c*x**2 + c)/(4*a**2*sqrt(-a**2*x**2 + 1)) - n*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(3)/2)*sqrt(-a**2*c*x**2 + c)/(12*a**3*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1321():
    f = x*sqrt(-a**2*c*x**2 + c)*exp(n*atanh(a*x))
    F = -2**(n/2 + sympy.S(3)/2)*n*(-a*x + 1)**(sympy.S(3)/2 - n/2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(3)/2 - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(3*a**2*(3 - n)*sqrt(-a**2*x**2 + 1)) - (-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(3)/2)*sqrt(-a**2*c*x**2 + c)/(3*a**2*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1322():
    f = sqrt(-a**2*c*x**2 + c)*exp(n*atanh(a*x))
    F = -2**(n/2 + sympy.S(3)/2)*(-a*x + 1)**(sympy.S(3)/2 - n/2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(3)/2 - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*(3 - n)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1323():
    f = sqrt(-a**2*c*x**2 + c)*exp(n*atanh(a*x))/x
    F = 2**(n/2 + sympy.S.Half)*n*(-a*x + 1)**(sympy.S(3)/2 - n/2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S.Half - n/2, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(sqrt(-a**2*x**2 + 1)*(n**2 - 4*n + 3)) + 2*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*c*x**2 + c)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/((1 - n)*sqrt(-a**2*x**2 + 1)) - (-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*c*x**2 + c)/((1 - n)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1324():
    f = sqrt(-a**2*c*x**2 + c)*exp(n*atanh(a*x))/x**2
    F = 2**(n/2 + sympy.S.Half)*a*(-a*x + 1)**(sympy.S.Half - n/2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S.Half - n/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), -a*x/2 + sympy.S.Half)/((1 - n)*sqrt(-a**2*x**2 + 1)) - 2*a*n*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*c*x**2 + c)*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (-a*x + 1)/(a*x + 1))/((1 - n)*sqrt(-a**2*x**2 + 1)) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*c*x**2 + c)/(x*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1325():
    f = (-a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(n*atanh(a*x))
    F = -2**(n/2 + sympy.S(5)/2)*c*(-a*x + 1)**(sympy.S(5)/2 - n/2)*sqrt(-a**2*c*x**2 + c)*hyper((sympy.S(5)/2 - n/2, -n/2 + sympy.S(-3)/2), (sympy.S(7)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*(5 - n)*sqrt(-a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1326():
    f = x**3*exp(n*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -2**(n/2 + sympy.S(-1)/2)*n*(n**2 + 5)*(-a*x + 1)**(sympy.S(3)/2 - n/2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half - n/2, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(3*a**4*(1 - n)*(3 - n)*sqrt(-a**2*c*x**2 + c)) - x**2*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(3*a**2*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)*(a*n*x*(1 - n) + n**2 + n + 4)/(6*a**4*(1 - n)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1327():
    f = x**2*exp(n*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -2**(n/2 + sympy.S.Half)*(n**2 + 1)*(-a*x + 1)**(sympy.S.Half - n/2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(3)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a**3*(1 - n**2)*sqrt(-a**2*c*x**2 + c)) - x*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(2*a**2*sqrt(-a**2*c*x**2 + c)) + (1 - n)*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(2*a**3*(n + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1328():
    f = x*exp(n*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -2**(n/2 + sympy.S(3)/2)*n*(-a*x + 1)**(sympy.S.Half - n/2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(3)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a**2*(1 - n**2)*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(a**2*(n + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1329():
    f = exp(n*atanh(a*x))/sqrt(-a**2*c*x**2 + c)
    F = -2**(n/2 + sympy.S.Half)*(-a*x + 1)**(sympy.S.Half - n/2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half - n/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a*(1 - n)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1330():
    f = exp(n*atanh(a*x))/(x*sqrt(-a**2*c*x**2 + c))
    F = -2*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (-a*x + 1)/(a*x + 1))/((1 - n)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1331():
    f = exp(n*atanh(a*x))/(x**2*sqrt(-a**2*c*x**2 + c))
    F = -2*a*n*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (-a*x + 1)/(a*x + 1))/((1 - n)*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(x*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1332():
    f = exp(n*atanh(a*x))/(x**3*sqrt(-a**2*c*x**2 + c))
    F = -a**2*(n**2 + 1)*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), (-a*x + 1)/(a*x + 1))/((1 - n)*sqrt(-a**2*c*x**2 + c)) - a*n*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(2*x*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S.Half)*sqrt(-a**2*x**2 + 1)/(2*x**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1333():
    f = x**3*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -2**(n/2 + sympy.S(-1)/2)*n*(-a*x + 1)**(sympy.S(3)/2 - n/2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S(3)/2 - n/2, sympy.S(3)/2 - n/2), (sympy.S(5)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a**4*c*(3 - n)*sqrt(-a**2*c*x**2 + c)) - x**2*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)/(a**2*c*sqrt(-a**2*c*x**2 + c)) + (-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*(-a*n*x*(2*n + 3) + n**2 + 2*n + 2)/(a**4*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1334():
    f = x**2*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = 2**(n/2 + sympy.S.Half)*(-a*x + 1)**(sympy.S.Half - n/2)*sqrt(-a**2*x**2 + 1)*hyper((sympy.S.Half - n/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), -a*x/2 + sympy.S.Half)/(a**3*c*(1 - n)*sqrt(-a**2*c*x**2 + c)) - (-a*x + n)*exp(n*atanh(a*x))/(a**3*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1335():
    f = x*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = (-a*n*x + 1)*exp(n*atanh(a*x))/(a**2*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1336():
    f = exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -(-a*x + n)*exp(n*atanh(a*x))/(a*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1337():
    f = exp(n*atanh(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = (-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)/(c*(n + 1)*sqrt(-a**2*c*x**2 + c)) - (n + 2)*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)/(c*(1 - n**2)*sqrt(-a**2*c*x**2 + c)) + 2*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/(c*(1 - n)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1338():
    f = exp(n*atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = 2*a*n*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/(c*(1 - n)*sqrt(-a**2*c*x**2 + c)) + a*(n + 2)*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)/(c*(n + 1)*sqrt(-a**2*c*x**2 + c)) - a*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*(n**2 + 2*n + 2)/(c*(1 - n**2)*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)/(c*x*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1339():
    f = exp(n*atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**(sympy.S(3)/2))
    F = a**2*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*(n**2 + 2*n + 3)/(2*c*(n + 1)*sqrt(-a**2*c*x**2 + c)) - a**2*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*(n**3 + 2*n**2 + 5*n + 6)/(2*c*(1 - n**2)*sqrt(-a**2*c*x**2 + c)) + a**2*(n**2 + 3)*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/(c*(1 - n)*sqrt(-a**2*c*x**2 + c)) - a*n*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)/(2*c*x*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)/(2*c*x**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1340():
    f = x**3*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = x**3*(-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(a*c**2*(n + 3)*sqrt(-a**2*c*x**2 + c)) - 3*x*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(a**3*c**2*(n + 3)*sqrt(-a**2*c*x**2 + c)) - (6 - 3*n)*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(a**4*c**2*(9 - n**2)*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*(-3*n**2 + 6*n + 3)/(a**4*c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9)) + (-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*(-3*n**2 + 6*n + 3)/(a**4*c**2*(3 - n)*(n + 1)*(n + 3)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1341():
    f = x**2*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -(-3*a*x + n)*exp(n*atanh(a*x))/(a**3*c*(9 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + (3 - n**2)*(-a*x + n)*exp(n*atanh(a*x))/(a**3*c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1342():
    f = x*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = n*(-3*a*x + n)*exp(n*atanh(a*x))/(3*a**2*c*(9 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + exp(n*atanh(a*x))/(3*a**2*c*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) + 2*n*(-a*x + n)*exp(n*atanh(a*x))/(a**2*c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1343():
    f = exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -(-3*a*x + n)*exp(n*atanh(a*x))/(a*c*(9 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 6*(-a*x + n)*exp(n*atanh(a*x))/(a*c**2*(1 - n**2)*(9 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1344():
    f = exp(n*atanh(a*x))/(x*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = (-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(-n**3 - 2*n**2 + 7*n + 18)/(c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9)) + (-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(c**2*(n + 3)*sqrt(-a**2*c*x**2 + c)) + (n + 6)*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(c**2*(n + 1)*(n + 3)*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(n**2 + 6*n + 15)/(c**2*(1 - n**2)*(n + 3)*sqrt(-a**2*c*x**2 + c)) + 2*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/(c**2*(1 - n)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1345():
    f = exp(n*atanh(a*x))/(x**2*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = 2*a*n*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/(c**2*(1 - n)*sqrt(-a**2*c*x**2 + c)) + a*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(-n**4 - 2*n**3 + 7*n**2 + 18*n + 24)/(c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9)) + a*(n + 4)*(-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(c**2*(n + 3)*sqrt(-a**2*c*x**2 + c)) + a*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(n**2 + 6*n + 12)/(c**2*(n + 1)*(n + 3)*sqrt(-a**2*c*x**2 + c)) - a*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(n**3 + 6*n**2 + 15*n + 24)/(c**2*(1 - n**2)*(n + 3)*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(c**2*x*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1346():
    f = exp(n*atanh(a*x))/(x**3*(-a**2*c*x**2 + c)**(sympy.S(5)/2))
    F = a**2*(-a*x + 1)**(sympy.S(3)/2 - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(-n**5 - 2*n**4 + 2*n**3 + 8*n**2 + 59*n + 90)/(2*c**2*sqrt(-a**2*c*x**2 + c)*(n**4 - 10*n**2 + 9)) + a**2*(-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(n**2 + 4*n + 5)/(2*c**2*(n + 3)*sqrt(-a**2*c*x**2 + c)) + a**2*(-a*x + 1)**(-n/2 + sympy.S(-1)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(n**3 + 6*n**2 + 17*n + 30)/(2*c**2*(n + 1)*(n + 3)*sqrt(-a**2*c*x**2 + c)) - a**2*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)*(n**4 + 6*n**3 + 20*n**2 + 54*n + 75)/(2*c**2*(1 - n**2)*(n + 3)*sqrt(-a**2*c*x**2 + c)) + a**2*(n**2 + 5)*(-a*x + 1)**(sympy.S.Half - n/2)*(a*x + 1)**(n/2 + sympy.S(-1)/2)*sqrt(-a**2*x**2 + 1)*hyper((1, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), (a*x + 1)/(-a*x + 1))/(c**2*(1 - n)*sqrt(-a**2*c*x**2 + c)) - a*n*(-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(2*c**2*x*sqrt(-a**2*c*x**2 + c)) - (-a*x + 1)**(-n/2 + sympy.S(-3)/2)*(a*x + 1)**(n/2 + sympy.S(-3)/2)*sqrt(-a**2*x**2 + 1)/(2*c**2*x**2*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1347():
    f = exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -(-5*a*x + n)*exp(n*atanh(a*x))/(a*c*(25 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(5)/2)) - 20*(-3*a*x + n)*exp(n*atanh(a*x))/(a*c**2*(9 - n**2)*(25 - n**2)*(-a**2*c*x**2 + c)**(sympy.S(3)/2)) - 120*(-a*x + n)*exp(n*atanh(a*x))/(a*c**3*(1 - n**2)*(9 - n**2)*(25 - n**2)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1348():
    f = x**m*(-a**2*c*x**2 + c)**2*exp(n*atanh(a*x))
    F = c**2*x**(m + 1)*appellf1(m + 1, -n/2 - 2, n/2 - 2, m + 2, -a*x, a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1349():
    f = x**m*(-a**2*c*x**2 + c)*exp(n*atanh(a*x))
    F = c*x**(m + 1)*appellf1(m + 1, -n/2 - 1, n/2 - 1, m + 2, -a*x, a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1350():
    f = x**m*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)
    F = x**(m + 1)*appellf1(m + 1, 1 - n/2, n/2 + 1, m + 2, -a*x, a*x)/(c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1351():
    f = x**m*exp(n*atanh(a*x))/(-a**2*c*x**2 + c)**2
    F = x**(m + 1)*appellf1(m + 1, 2 - n/2, n/2 + 2, m + 2, -a*x, a*x)/(c**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1352():
    f = x**m*(-a**2*c*x**2 + c)**p*exp(n*atanh(a*x))
    F = x**(m + 1)*(-a**2*c*x**2 + c)**p*appellf1(m + 1, -n/2 - p, n/2 - p, m + 2, -a*x, a*x)/((m + 1)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1353():
    f = x*(-a**2*c*x**2 + c)**p*exp(n*atanh(a*x))
    F = -2**(n/2 + p)*n*(-a*x + 1)**(-n/2 + p + 1)*(-a**2*c*x**2 + c)**p*hyper((-n/2 + p + 1, -n/2 - p), (-n/2 + p + 2,), -a*x/2 + sympy.S.Half)/(a**2*(p + 1)*(-a**2*x**2 + 1)**p*(-n + 2*p + 2)) - (-a*x + 1)**(-n/2 + p + 1)*(a*x + 1)**(n/2 + p + 1)*(-a**2*c*x**2 + c)**p/(2*a**2*(p + 1)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1354():
    f = (-a**2*c*x**2 + c)**p*exp(n*atanh(a*x))
    F = -2**(n/2 + p + 1)*(-a*x + 1)**(-n/2 + p + 1)*(-a**2*c*x**2 + c)**p*hyper((-n/2 + p + 1, -n/2 - p), (-n/2 + p + 2,), -a*x/2 + sympy.S.Half)/(a*(-a**2*x**2 + 1)**p*(-n + 2*p + 2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1355():
    f = exp((2*p + 2)*atanh(a*x))/(-a**2*x**2 + 1)**p
    F = (-a*x + 1)**(1 - 2*p)/(a*(1 - 2*p)) + 1/(a*p*(-a*x + 1)**(2*p))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1356():
    f = exp((2*p + 2)*atanh(a*x))/(-a**2*c*x**2 + c)**p
    F = (-a*x + 1)**(1 - 2*p)*(-a**2*x**2 + 1)**p/(a*(1 - 2*p)*(-a**2*c*x**2 + c)**p) + (-a**2*x**2 + 1)**p/(a*p*(-a*x + 1)**(2*p)*(-a**2*c*x**2 + c)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1357():
    f = (-a**2*c*x**2 + c)**p*exp(2*p*atanh(a*x))
    F = (a*x + 1)**(2*p + 1)*(-a**2*c*x**2 + c)**p/(a*(2*p + 1)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1358():
    f = (-a**2*c*x**2 + c)**p*exp(-2*p*atanh(a*x))
    F = -(-a*x + 1)**(2*p + 1)*(-a**2*c*x**2 + c)**p/(a*(2*p + 1)*(-a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1359():
    f = x**2*(-a**2*c*x**2 + c)**(-n**2/2 - 1)*exp(n*atanh(a*x))
    F = (-a*n*x + 1)*exp(n*atanh(a*x))/(a**3*c*n*(1 - n**2)*(-a**2*c*x**2 + c)**(n**2/2))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1360():
    f = x**2*exp(6*atanh(a*x))/(-a**2*c*x**2 + c)**19
    F = -(-6*a*x + 1)/(210*a**3*c**19*(-a*x + 1)**21*(a*x + 1)**15)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1361():
    f = x**2*exp(4*atanh(a*x))/(-a**2*c*x**2 + c)**9
    F = -(-4*a*x + 1)/(60*a**3*c**9*(-a*x + 1)**10*(a*x + 1)**6)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1362():
    f = x**2*exp(2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = -(-2*a*x + 1)/(6*a**3*c**3*(-a*x + 1)**3*(a*x + 1))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1363():
    f = x**2*exp(-2*atanh(a*x))/(-a**2*c*x**2 + c)**3
    F = (2*a*x + 1)/(6*a**3*c**3*(-a*x + 1)*(a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1364():
    f = x**2*exp(-4*atanh(a*x))/(-a**2*c*x**2 + c)**9
    F = (4*a*x + 1)/(60*a**3*c**9*(-a*x + 1)**6*(a*x + 1)**10)
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1365():
    f = x**2*exp(5*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(27)/2)
    F = -(-5*a*x + 1)*sqrt(-a**2*x**2 + 1)/(120*a**3*c**13*(-a*x + 1)**15*(a*x + 1)**10*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1366():
    f = x**2*exp(3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(11)/2)
    F = -(-3*a*x + 1)*sqrt(-a**2*x**2 + 1)/(24*a**3*c**5*(-a*x + 1)**6*(a*x + 1)**3*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1367():
    f = x**2*exp(atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = 3*sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*a**3*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*a**3*c*sqrt(-a**2*c*x**2 + c)) + sqrt(-a**2*x**2 + 1)/(2*a**3*c*(-a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1368():
    f = x**2*exp(-atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -sqrt(-a**2*x**2 + 1)*log(-a*x + 1)/(4*a**3*c*sqrt(-a**2*c*x**2 + c)) - 3*sqrt(-a**2*x**2 + 1)*log(a*x + 1)/(4*a**3*c*sqrt(-a**2*c*x**2 + c)) - sqrt(-a**2*x**2 + 1)/(2*a**3*c*(a*x + 1)*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1369():
    f = x**2*exp(-3*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(11)/2)
    F = (3*a*x + 1)*sqrt(-a**2*x**2 + 1)/(24*a**3*c**5*(-a*x + 1)**3*(a*x + 1)**6*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_7_Inverse_hyperbolic_functions_7_3_Inverse_hyperbolic_tangent_7_3_6_Exponentials_of_inverse_hyperbolic_tangent_functions_1370():
    f = x**2*exp(-5*atanh(a*x))/(-a**2*c*x**2 + c)**(sympy.S(27)/2)
    F = (5*a*x + 1)*sqrt(-a**2*x**2 + 1)/(120*a**3*c**13*(-a*x + 1)**10*(a*x + 1)**15*sqrt(-a**2*c*x**2 + c))
    assert integrate(f, x) == F

