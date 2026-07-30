"""Generated from MathematicaSyntaxTestSuite.

Source: 5 Inverse trig functions/5.3 Inverse tangent/5.3.6 Exponentials of inverse tangent.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, m, n, p = symbols('a b c m n p')

def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_1():
    f = x**4*exp(I*atan(a*x))
    F = I*x**4*sqrt(a**2*x**2 + 1)/(5*a) + x**3*sqrt(a**2*x**2 + 1)/(4*a**2) - 4*I*x**2*sqrt(a**2*x**2 + 1)/(15*a**3) + (-45*a*x + 64*I)*sqrt(a**2*x**2 + 1)/(120*a**5) + 3*asinh(a*x)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_2():
    f = x**3*exp(I*atan(a*x))
    F = I*x**3*sqrt(a**2*x**2 + 1)/(4*a) + x**2*sqrt(a**2*x**2 + 1)/(3*a**2) - sqrt(a**2*x**2 + 1)*(9*I*a*x + 16)/(24*a**4) + 3*I*asinh(a*x)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_3():
    f = x**2*exp(I*atan(a*x))
    F = x*sqrt(a**2*x**2 + 1)/(2*a**2) + I*(a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**3) - I*sqrt(a**2*x**2 + 1)/a**3 - asinh(a*x)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_4():
    f = x*exp(I*atan(a*x))
    F = sqrt(a**2*x**2 + 1)*(I*a*x + 2)/(2*a**2) - I*asinh(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_5():
    f = exp(I*atan(a*x))
    F = I*sqrt(a**2*x**2 + 1)/a + asinh(a*x)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_6():
    f = exp(I*atan(a*x))/x
    F = I*asinh(a*x) - atanh(sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_7():
    f = exp(I*atan(a*x))/x**2
    F = -I*a*atanh(sqrt(a**2*x**2 + 1)) - sqrt(a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_8():
    f = exp(I*atan(a*x))/x**3
    F = a**2*atanh(sqrt(a**2*x**2 + 1))/2 - I*a*sqrt(a**2*x**2 + 1)/x - sqrt(a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_9():
    f = exp(I*atan(a*x))/x**4
    F = I*a**3*atanh(sqrt(a**2*x**2 + 1))/2 + 2*a**2*sqrt(a**2*x**2 + 1)/(3*x) - I*a*sqrt(a**2*x**2 + 1)/(2*x**2) - sqrt(a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_10():
    f = exp(I*atan(a*x))/x**5
    F = -3*a**4*atanh(sqrt(a**2*x**2 + 1))/8 + 2*I*a**3*sqrt(a**2*x**2 + 1)/(3*x) + 3*a**2*sqrt(a**2*x**2 + 1)/(8*x**2) - I*a*sqrt(a**2*x**2 + 1)/(3*x**3) - sqrt(a**2*x**2 + 1)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_11():
    f = x**3*exp(2*I*atan(a*x))
    F = -x**4/4 + 2*I*x**3/(3*a) + x**2/a**2 - 2*I*x/a**3 - 2*log(a*x + I)/a**4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_12():
    f = x**2*exp(2*I*atan(a*x))
    F = -x**3/3 + I*x**2/a + 2*x/a**2 - 2*I*log(a*x + I)/a**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_13():
    f = x*exp(2*I*atan(a*x))
    F = -x**2/2 + 2*I*x/a + 2*log(a*x + I)/a**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_14():
    f = exp(2*I*atan(a*x))
    F = -x + 2*I*log(a*x + I)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_15():
    f = exp(2*I*atan(a*x))/x
    F = log(x) - 2*log(a*x + I)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_16():
    f = exp(2*I*atan(a*x))/x**2
    F = 2*I*a*log(x) - 2*I*a*log(a*x + I) - 1/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_17():
    f = exp(2*I*atan(a*x))/x**3
    F = -2*a**2*log(x) + 2*a**2*log(a*x + I) - 2*I*a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_18():
    f = exp(2*I*atan(a*x))/x**4
    F = -2*I*a**3*log(x) + 2*I*a**3*log(a*x + I) + 2*a**2/x - I*a/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_19():
    f = x**3*exp(3*I*atan(a*x))
    F = -I*x**3*sqrt(a**2*x**2 + 1)/(4*a) - x**2*sqrt(a**2*x**2 + 1)/a**2 - 9*I*(-3*a*x + 2*I)*sqrt(a**2*x**2 + 1)/(8*a**4) + 27*sqrt(a**2*x**2 + 1)/(4*a**4) - 51*I*asinh(a*x)/(8*a**4) + (I*a*x + 1)**3/(a**4*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_20():
    f = x**2*exp(3*I*atan(a*x))
    F = (-3*a*x + 28*I)*sqrt(a**2*x**2 + 1)/(6*a**3) + I*sqrt(a**2*x**2 + 1)*(I*a*x + 3)**2/(3*a**3) + 11*asinh(a*x)/(2*a**3) + I*(I*a*x + 1)**3/(a**3*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_21():
    f = x*exp(3*I*atan(a*x))
    F = -(a**2*x**2 + 1)**(sympy.S(5)/2)/(a**2*(-I*a*x + 1)**3) - 3*(a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**2*(-I*a*x + 1)) - 9*sqrt(a**2*x**2 + 1)/(2*a**2) + 9*I*asinh(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_22():
    f = exp(3*I*atan(a*x))
    F = -3*I*sqrt(a**2*x**2 + 1)/a - 3*asinh(a*x)/a - 2*I*(I*a*x + 1)**2/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_23():
    f = exp(3*I*atan(a*x))/x
    F = -I*asinh(a*x) - atanh(sqrt(a**2*x**2 + 1)) + 4*I*sqrt(a**2*x**2 + 1)/(a*x + I)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_24():
    f = exp(3*I*atan(a*x))/x**2
    F = -3*I*a*atanh(sqrt(a**2*x**2 + 1)) - 4*a*sqrt(a**2*x**2 + 1)/(a*x + I) - sqrt(a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_25():
    f = exp(3*I*atan(a*x))/x**3
    F = 9*a**2*atanh(sqrt(a**2*x**2 + 1))/2 - 4*I*a**2*sqrt(a**2*x**2 + 1)/(a*x + I) - 3*I*a*sqrt(a**2*x**2 + 1)/x - sqrt(a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_26():
    f = exp(3*I*atan(a*x))/x**4
    F = 11*I*a**3*atanh(sqrt(a**2*x**2 + 1))/2 + 4*a**3*sqrt(a**2*x**2 + 1)/(a*x + I) + 14*a**2*sqrt(a**2*x**2 + 1)/(3*x) - 3*I*a*sqrt(a**2*x**2 + 1)/(2*x**2) - sqrt(a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_27():
    f = x**3*exp(4*I*atan(a*x))
    F = x**4/4 - 4*I*x**3/(3*a) - 4*x**2/a**2 + 12*I*x/a**3 + 16*log(a*x + I)/a**4 + 4*I/(a**4*(a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_28():
    f = x**2*exp(4*I*atan(a*x))
    F = x**3/3 - 2*I*x**2/a - 8*x/a**2 + 12*I*log(a*x + I)/a**3 - 4/(a**3*(a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_29():
    f = x*exp(4*I*atan(a*x))
    F = x**2/2 - 4*I*x/a - 8*log(a*x + I)/a**2 - 4*I/(a**2*(a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_30():
    f = exp(4*I*atan(a*x))
    F = x - 4*I*log(a*x + I)/a + 4/(a*(a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_31():
    f = exp(4*I*atan(a*x))/x
    F = log(x) + 4*I/(a*x + I)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_32():
    f = exp(4*I*atan(a*x))/x**2
    F = 4*I*a*log(x) - 4*I*a*log(a*x + I) - 4*a/(a*x + I) - 1/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_33():
    f = exp(4*I*atan(a*x))/x**3
    F = -8*a**2*log(x) + 8*a**2*log(a*x + I) - 4*I*a**2/(a*x + I) - 4*I*a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_34():
    f = exp(4*I*atan(a*x))/x**4
    F = -12*I*a**3*log(x) + 12*I*a**3*log(a*x + I) + 4*a**3/(a*x + I) + 8*a**2/x - 2*I*a/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_35():
    f = x**3*exp(-I*atan(a*x))
    F = -I*x**3*sqrt(a**2*x**2 + 1)/(4*a) + x**2*sqrt(a**2*x**2 + 1)/(3*a**2) - sqrt(a**2*x**2 + 1)*(-9*I*a*x + 16)/(24*a**4) - 3*I*asinh(a*x)/(8*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_36():
    f = x**2*exp(-I*atan(a*x))
    F = x*sqrt(a**2*x**2 + 1)/(2*a**2) - I*(a**2*x**2 + 1)**(sympy.S(3)/2)/(3*a**3) + I*sqrt(a**2*x**2 + 1)/a**3 - asinh(a*x)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_37():
    f = x*exp(-I*atan(a*x))
    F = sqrt(a**2*x**2 + 1)*(-I*a*x + 2)/(2*a**2) + I*asinh(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_38():
    f = exp(-I*atan(a*x))
    F = -I*sqrt(a**2*x**2 + 1)/a + asinh(a*x)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_39():
    f = exp(-I*atan(a*x))/x
    F = -I*asinh(a*x) - atanh(sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_40():
    f = exp(-I*atan(a*x))/x**2
    F = I*a*atanh(sqrt(a**2*x**2 + 1)) - sqrt(a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_41():
    f = exp(-I*atan(a*x))/x**3
    F = a**2*atanh(sqrt(a**2*x**2 + 1))/2 + I*a*sqrt(a**2*x**2 + 1)/x - sqrt(a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_42():
    f = exp(-I*atan(a*x))/x**4
    F = -I*a**3*atanh(sqrt(a**2*x**2 + 1))/2 + 2*a**2*sqrt(a**2*x**2 + 1)/(3*x) + I*a*sqrt(a**2*x**2 + 1)/(2*x**2) - sqrt(a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_43():
    f = exp(-I*atan(a*x))/x**5
    F = -3*a**4*atanh(sqrt(a**2*x**2 + 1))/8 - 2*I*a**3*sqrt(a**2*x**2 + 1)/(3*x) + 3*a**2*sqrt(a**2*x**2 + 1)/(8*x**2) + I*a*sqrt(a**2*x**2 + 1)/(3*x**3) - sqrt(a**2*x**2 + 1)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_44():
    f = x**3*exp(-2*I*atan(a*x))
    F = -x**4/4 - 2*I*x**3/(3*a) + x**2/a**2 + 2*I*x/a**3 - 2*log(-a*x + I)/a**4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_45():
    f = x**2*exp(-2*I*atan(a*x))
    F = -x**3/3 - I*x**2/a + 2*x/a**2 + 2*I*log(-a*x + I)/a**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_46():
    f = x*exp(-2*I*atan(a*x))
    F = -x**2/2 - 2*I*x/a + 2*log(-a*x + I)/a**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_47():
    f = exp(-2*I*atan(a*x))
    F = -x - 2*I*log(-a*x + I)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_48():
    f = exp(-2*I*atan(a*x))/x
    F = log(x) - 2*log(-a*x + I)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_49():
    f = exp(-2*I*atan(a*x))/x**2
    F = -2*I*a*log(x) + 2*I*a*log(-a*x + I) - 1/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_50():
    f = exp(-2*I*atan(a*x))/x**3
    F = -2*a**2*log(x) + 2*a**2*log(-a*x + I) + 2*I*a/x - 1/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_51():
    f = exp(-2*I*atan(a*x))/x**4
    F = 2*I*a**3*log(x) - 2*I*a**3*log(-a*x + I) + 2*a**2/x + I*a/x**2 - 1/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_52():
    f = x**3*exp(-3*I*atan(a*x))
    F = I*x**3*sqrt(a**2*x**2 + 1)/(4*a) - x**2*sqrt(a**2*x**2 + 1)/a**2 - 9*I*(3*a*x + 2*I)*sqrt(a**2*x**2 + 1)/(8*a**4) + 27*sqrt(a**2*x**2 + 1)/(4*a**4) + 51*I*asinh(a*x)/(8*a**4) + (-I*a*x + 1)**3/(a**4*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_53():
    f = x**2*exp(-3*I*atan(a*x))
    F = -(3*a*x + 28*I)*sqrt(a**2*x**2 + 1)/(6*a**3) - I*sqrt(a**2*x**2 + 1)*(-I*a*x + 3)**2/(3*a**3) + 11*asinh(a*x)/(2*a**3) - I*(-I*a*x + 1)**3/(a**3*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_54():
    f = x*exp(-3*I*atan(a*x))
    F = -(a**2*x**2 + 1)**(sympy.S(5)/2)/(a**2*(I*a*x + 1)**3) - 3*(a**2*x**2 + 1)**(sympy.S(3)/2)/(2*a**2*(I*a*x + 1)) - 9*sqrt(a**2*x**2 + 1)/(2*a**2) - 9*I*asinh(a*x)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_55():
    f = exp(-3*I*atan(a*x))
    F = 3*I*sqrt(a**2*x**2 + 1)/a - 3*asinh(a*x)/a + 2*I*(-I*a*x + 1)**2/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_56():
    f = exp(-3*I*atan(a*x))/x
    F = I*asinh(a*x) - atanh(sqrt(a**2*x**2 + 1)) + 4*I*sqrt(a**2*x**2 + 1)/(-a*x + I)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_57():
    f = exp(-3*I*atan(a*x))/x**2
    F = 3*I*a*atanh(sqrt(a**2*x**2 + 1)) + 4*a*sqrt(a**2*x**2 + 1)/(-a*x + I) - sqrt(a**2*x**2 + 1)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_58():
    f = exp(-3*I*atan(a*x))/x**3
    F = 9*a**2*atanh(sqrt(a**2*x**2 + 1))/2 - 4*I*a**2*sqrt(a**2*x**2 + 1)/(-a*x + I) + 3*I*a*sqrt(a**2*x**2 + 1)/x - sqrt(a**2*x**2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_59():
    f = exp(-3*I*atan(a*x))/x**4
    F = -11*I*a**3*atanh(sqrt(a**2*x**2 + 1))/2 - 4*a**3*sqrt(a**2*x**2 + 1)/(-a*x + I) + 14*a**2*sqrt(a**2*x**2 + 1)/(3*x) + 3*I*a*sqrt(a**2*x**2 + 1)/(2*x**2) - sqrt(a**2*x**2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_60():
    f = exp(-3*I*atan(a*x))/x**5
    F = -51*a**4*atanh(sqrt(a**2*x**2 + 1))/8 + 4*I*a**4*sqrt(a**2*x**2 + 1)/(-a*x + I) - 6*I*a**3*sqrt(a**2*x**2 + 1)/x + 19*a**2*sqrt(a**2*x**2 + 1)/(8*x**2) + I*a*sqrt(a**2*x**2 + 1)/x**3 - sqrt(a**2*x**2 + 1)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_61():
    f = x**2*exp(I*atan(a*x)/2)
    F = x*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(3*a**2) - I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(12*a**3) - 3*I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(8*a**3) - 3*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) + 3*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) - 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_62():
    f = x*exp(I*atan(a*x)/2)
    F = (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(2*a**2) + (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*a**2) + sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) - sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_63():
    f = exp(I*atan(a*x)/2)
    F = I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/a + sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) - sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_64():
    f = exp(I*atan(a*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 + sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 - 2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_65():
    f = exp(I*atan(a*x)/2)/x**2
    F = -I*a*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - I*a*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_66():
    f = exp(I*atan(a*x)/2)/x**3
    F = a**2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + a**2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 - I*a*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*x) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_67():
    f = exp(I*atan(a*x)/2)/x**4
    F = 3*I*a**3*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 3*I*a**3*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 11*a**2*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(24*x) - 5*I*a*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(12*x**2) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_68():
    f = exp(I*atan(a*x)/2)/x**5
    F = -11*a**4*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 11*a**4*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 + 83*I*a**3*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(192*x) + 29*a**2*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(96*x**2) - 7*I*a*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(24*x**3) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_69():
    f = exp(I*atan(a*x)/2)/x**6
    F = -31*I*a**5*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/128 - 31*I*a**5*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/128 - 611*a**4*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(1920*x) + 269*I*a**3*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(960*x**2) + 11*a**2*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(48*x**3) - 9*I*a*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(40*x**4) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(5*x**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_70():
    f = x**3*exp(3*I*atan(a*x)/2)
    F = x**2*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(7)/4)/(4*a**2) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(7)/4)*(4*I*a*x + 11)/(32*a**4) - 41*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(64*a**4) + 123*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) - 123*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) - 123*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) - 123*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_71():
    f = x**2*exp(3*I*atan(a*x)/2)
    F = x*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(7)/4)/(3*a**2) - I*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(7)/4)/(4*a**3) - 17*I*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(24*a**3) + 17*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) - 17*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) - 17*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 17*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_72():
    f = x*exp(3*I*atan(a*x)/2)
    F = (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(7)/4)/(2*a**2) + 3*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*a**2) - 9*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) + 9*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) + 9*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) + 9*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_73():
    f = exp(3*I*atan(a*x)/2)
    F = I*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/a - 3*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) + 3*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_74():
    f = exp(3*I*atan(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 + 2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_75():
    f = exp(3*I*atan(a*x)/2)/x**2
    F = 3*I*a*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - 3*I*a*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_76():
    f = exp(3*I*atan(a*x)/2)/x**3
    F = -9*a**2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + 9*a**2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 - 3*I*a*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*x) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(7)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_77():
    f = exp(3*I*atan(a*x)/2)/x**4
    F = -17*I*a**3*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 17*I*a**3*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 23*a**2*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(24*x) - 7*I*a*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(12*x**2) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_78():
    f = exp(3*I*atan(a*x)/2)/x**5
    F = 123*a**4*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 123*a**4*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 + 63*I*a**3*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(64*x) + 15*a**2*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(32*x**2) - 3*I*a*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(8*x**3) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_79():
    f = x**3*exp(5*I*atan(a*x)/2)
    F = -4*I*x**3*(I*a*x + 1)**(sympy.S(5)/4)/(a*(-I*a*x + 1)**(sympy.S(1)/4)) - 17*x**2*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(4*a**2) - I*(-452*a*x + 521*I)*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(96*a**4) + 475*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(64*a**4) + 475*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) - 475*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) + 475*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) + 475*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_80():
    f = x**2*exp(5*I*atan(a*x)/2)
    F = I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(9)/4)/(3*a**3) + 11*I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(4*a**3) + 55*I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(8*a**3) + 55*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) - 55*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) + 55*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) + 55*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3) + 2*I*(I*a*x + 1)**(sympy.S(9)/4)/(a**3*(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_81():
    f = x*exp(5*I*atan(a*x)/2)
    F = -5*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(5)/4)/(2*a**2) - 25*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*a**2) - 25*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) + 25*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) - 25*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) - 25*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2) - 2*(I*a*x + 1)**(sympy.S(9)/4)/(a**2*(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_82():
    f = exp(5*I*atan(a*x)/2)
    F = -5*I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/a - 5*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) + 5*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) - 5*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) - 5*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(2*a) - 4*I*(I*a*x + 1)**(sympy.S(5)/4)/(a*(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_83():
    f = exp(5*I*atan(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 - 2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + 8*(I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_84():
    f = exp(5*I*atan(a*x)/2)/x**2
    F = -5*I*a*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - 5*I*a*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + 10*I*a*(I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4) - (I*a*x + 1)**(sympy.S(5)/4)/(x*(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_85():
    f = exp(5*I*atan(a*x)/2)/x**3
    F = 25*a**2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + 25*a**2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 - 25*a**2*(I*a*x + 1)**(sympy.S(1)/4)/(2*(-I*a*x + 1)**(sympy.S(1)/4)) - 5*I*a*(I*a*x + 1)**(sympy.S(5)/4)/(4*x*(-I*a*x + 1)**(sympy.S(1)/4)) - (I*a*x + 1)**(sympy.S(9)/4)/(2*x**2*(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_86():
    f = exp(5*I*atan(a*x)/2)/x**4
    F = 55*I*a**3*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 55*I*a**3*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 - 287*I*a**3*(I*a*x + 1)**(sympy.S(1)/4)/(24*(-I*a*x + 1)**(sympy.S(1)/4)) + 61*a**2*(I*a*x + 1)**(sympy.S(1)/4)/(24*x*(-I*a*x + 1)**(sympy.S(1)/4)) - 13*I*a*(I*a*x + 1)**(sympy.S(1)/4)/(12*x**2*(-I*a*x + 1)**(sympy.S(1)/4)) - (I*a*x + 1)**(sympy.S(1)/4)/(3*x**3*(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_87():
    f = exp(5*I*atan(a*x)/2)/x**5
    F = -475*a**4*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 475*a**4*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 + 2467*a**4*(I*a*x + 1)**(sympy.S(1)/4)/(192*(-I*a*x + 1)**(sympy.S(1)/4)) + 521*I*a**3*(I*a*x + 1)**(sympy.S(1)/4)/(192*x*(-I*a*x + 1)**(sympy.S(1)/4)) + 113*a**2*(I*a*x + 1)**(sympy.S(1)/4)/(96*x**2*(-I*a*x + 1)**(sympy.S(1)/4)) - 17*I*a*(I*a*x + 1)**(sympy.S(1)/4)/(24*x**3*(-I*a*x + 1)**(sympy.S(1)/4)) - (I*a*x + 1)**(sympy.S(1)/4)/(4*x**4*(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_88():
    f = x**3*exp(-I*atan(a*x)/2)
    F = x**2*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*a**2) - (-4*I*a*x + 25)*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(96*a**4) - 11*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(64*a**4) - 11*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) + 11*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) + 11*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) + 11*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_89():
    f = x**2*exp(-I*atan(a*x)/2)
    F = x*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(3*a**2) + I*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(12*a**3) + 3*I*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(8*a**3) + 3*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) - 3*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) - 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_90():
    f = x*exp(-I*atan(a*x)/2)
    F = (-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(2*a**2) + (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*a**2) + sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) - sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_91():
    f = exp(-I*atan(a*x)/2)
    F = -I*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/a - sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) + sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) + sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_92():
    f = exp(-I*atan(a*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 + sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 + 2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_93():
    f = exp(-I*atan(a*x)/2)/x**2
    F = -I*a*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + I*a*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_94():
    f = exp(-I*atan(a*x)/2)/x**3
    F = -a**2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + a**2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + I*a*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*x) - (-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_95():
    f = exp(-I*atan(a*x)/2)/x**4
    F = 3*I*a**3*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 - 3*I*a**3*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 11*a**2*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(24*x) + 5*I*a*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(12*x**2) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_96():
    f = exp(-I*atan(a*x)/2)/x**5
    F = 11*a**4*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 11*a**4*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 83*I*a**3*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(192*x) + 29*a**2*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(96*x**2) + 7*I*a*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(24*x**3) - (-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_97():
    f = x**3*exp(-3*I*atan(a*x)/2)
    F = x**2*(-I*a*x + 1)**(sympy.S(7)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*a**2) - (-4*I*a*x + 11)*(-I*a*x + 1)**(sympy.S(7)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(32*a**4) - 41*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(64*a**4) + 123*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) - 123*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) + 123*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) + 123*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_98():
    f = x**2*exp(-3*I*atan(a*x)/2)
    F = x*(-I*a*x + 1)**(sympy.S(7)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(3*a**2) + I*(-I*a*x + 1)**(sympy.S(7)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*a**3) + 17*I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(24*a**3) - 17*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) + 17*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) - 17*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) - 17*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_99():
    f = x*exp(-3*I*atan(a*x)/2)
    F = (-I*a*x + 1)**(sympy.S(7)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(2*a**2) + 3*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*a**2) - 9*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) + 9*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) - 9*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) - 9*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_100():
    f = exp(-3*I*atan(a*x)/2)
    F = -I*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/a + 3*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) - 3*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_101():
    f = exp(-3*I*atan(a*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 - 2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_102():
    f = exp(-3*I*atan(a*x)/2)/x**2
    F = 3*I*a*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + 3*I*a*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_103():
    f = exp(-3*I*atan(a*x)/2)/x**3
    F = 9*a**2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + 9*a**2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + 3*I*a*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*x) - (-I*a*x + 1)**(sympy.S(7)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_104():
    f = exp(-3*I*atan(a*x)/2)/x**4
    F = -17*I*a**3*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 - 17*I*a**3*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 23*a**2*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(24*x) + 7*I*a*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(12*x**2) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_105():
    f = exp(-3*I*atan(a*x)/2)/x**5
    F = -123*a**4*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 123*a**4*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 63*I*a**3*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(64*x) + 15*a**2*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(32*x**2) + 3*I*a*(-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(8*x**3) - (-I*a*x + 1)**(sympy.S(3)/4)*(I*a*x + 1)**(sympy.S(1)/4)/(4*x**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_106():
    f = x**3*exp(-5*I*atan(a*x)/2)
    F = 4*I*x**3*(-I*a*x + 1)**(sympy.S(5)/4)/(a*(I*a*x + 1)**(sympy.S(1)/4)) - 17*x**2*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*a**2) - I*(452*a*x + 521*I)*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(96*a**4) + 475*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(64*a**4) + 475*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) - 475*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(256*a**4) - 475*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(128*a**4) - 475*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(128*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_107():
    f = x**2*exp(-5*I*atan(a*x)/2)
    F = -I*(-I*a*x + 1)**(sympy.S(9)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(3*a**3) - 2*I*(-I*a*x + 1)**(sympy.S(9)/4)/(a**3*(I*a*x + 1)**(sympy.S(1)/4)) - 11*I*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*a**3) - 55*I*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(8*a**3) - 55*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) + 55*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(32*a**3) + 55*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(16*a**3) + 55*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(16*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_108():
    f = x*exp(-5*I*atan(a*x)/2)
    F = -2*(-I*a*x + 1)**(sympy.S(9)/4)/(a**2*(I*a*x + 1)**(sympy.S(1)/4)) - 5*(-I*a*x + 1)**(sympy.S(5)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(2*a**2) - 25*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/(4*a**2) - 25*sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) + 25*sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(16*a**2) + 25*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(8*a**2) + 25*sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_109():
    f = exp(-5*I*atan(a*x)/2)
    F = 4*I*(-I*a*x + 1)**(sympy.S(5)/4)/(a*(I*a*x + 1)**(sympy.S(1)/4)) + 5*I*(-I*a*x + 1)**(sympy.S(1)/4)*(I*a*x + 1)**(sympy.S(3)/4)/a + 5*sqrt(2)*I*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) - 5*sqrt(2)*I*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/(4*a) - 5*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1)/(2*a) - 5*sqrt(2)*I*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(2*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_110():
    f = exp(-5*I*atan(a*x)/2)/x
    F = 8*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(2)*log(-sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + sqrt(-I*a*x + 1)/sqrt(I*a*x + 1) + 1)/2 + 2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1) - 2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_111():
    f = exp(-5*I*atan(a*x)/2)/x**2
    F = -10*I*a*(-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) - 5*I*a*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) + 5*I*a*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(5)/4)/(x*(I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_112():
    f = exp(-5*I*atan(a*x)/2)/x**3
    F = -25*a**2*(-I*a*x + 1)**(sympy.S(1)/4)/(2*(I*a*x + 1)**(sympy.S(1)/4)) - 25*a**2*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + 25*a**2*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/4 + 5*I*a*(-I*a*x + 1)**(sympy.S(5)/4)/(4*x*(I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(9)/4)/(2*x**2*(I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_113():
    f = exp(-5*I*atan(a*x)/2)/x**4
    F = 287*I*a**3*(-I*a*x + 1)**(sympy.S(1)/4)/(24*(I*a*x + 1)**(sympy.S(1)/4)) + 55*I*a**3*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 - 55*I*a**3*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/8 + 61*a**2*(-I*a*x + 1)**(sympy.S(1)/4)/(24*x*(I*a*x + 1)**(sympy.S(1)/4)) + 13*I*a*(-I*a*x + 1)**(sympy.S(1)/4)/(12*x**2*(I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(1)/4)/(3*x**3*(I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_114():
    f = exp(-5*I*atan(a*x)/2)/x**5
    F = 2467*a**4*(-I*a*x + 1)**(sympy.S(1)/4)/(192*(I*a*x + 1)**(sympy.S(1)/4)) + 475*a**4*atan((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 475*a**4*atanh((I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4))/64 - 521*I*a**3*(-I*a*x + 1)**(sympy.S(1)/4)/(192*x*(I*a*x + 1)**(sympy.S(1)/4)) + 113*a**2*(-I*a*x + 1)**(sympy.S(1)/4)/(96*x**2*(I*a*x + 1)**(sympy.S(1)/4)) + 17*I*a*(-I*a*x + 1)**(sympy.S(1)/4)/(24*x**3*(I*a*x + 1)**(sympy.S(1)/4)) - (-I*a*x + 1)**(sympy.S(1)/4)/(4*x**4*(I*a*x + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_115():
    f = x**2*exp(I*atan(x)/3)
    F = x*(-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(7)/6)/3 - I*(-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(7)/6)/18 - 19*I*(-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6)/54 - 19*sqrt(3)*I*log(-sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/324 + 19*sqrt(3)*I*log(sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/324 - 19*I*atan((-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6))/81 - 19*I*atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) - sqrt(3))/162 - 19*I*atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + sqrt(3))/162
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_116():
    f = x*exp(I*atan(x)/3)
    F = (-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(7)/6)/2 + (-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6)/6 + sqrt(3)*log(-sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/36 - sqrt(3)*log(sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/36 + atan((-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6))/9 + atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) - sqrt(3))/18 + atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + sqrt(3))/18
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_117():
    f = exp(I*atan(x)/3)
    F = I*(-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6) + sqrt(3)*I*log(-sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/6 - sqrt(3)*I*log(sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/6 + 2*I*atan((-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6))/3 + I*atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) - sqrt(3))/3 + I*atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + sqrt(3))/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_118():
    f = exp(I*atan(x)/3)/x
    F = log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) - (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/2 - log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) + (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/2 - sqrt(3)*log(-sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/2 + sqrt(3)*log(sqrt(3)*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + (-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/2 + sqrt(3)*atan(sqrt(3)*(1 - 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3) - sqrt(3)*atan(sqrt(3)*(1 + 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3) - 2*atan((-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6)) - atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) - sqrt(3)) - atan(2*(-I*x + 1)**(sympy.S(1)/6)/(I*x + 1)**(sympy.S(1)/6) + sqrt(3)) - 2*atanh((I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_119():
    f = exp(I*atan(x)/3)/x**2
    F = I*log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) - (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/6 - I*log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) + (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/6 + sqrt(3)*I*atan(sqrt(3)*(1 - 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3)/3 - sqrt(3)*I*atan(sqrt(3)*(1 + 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3)/3 - 2*I*atanh((I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3 - (-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_120():
    f = exp(I*atan(x)/3)/x**3
    F = -log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) - (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/36 + log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) + (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/36 - sqrt(3)*atan(sqrt(3)*(1 - 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3)/18 + sqrt(3)*atan(sqrt(3)*(1 + 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3)/18 + atanh((I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/9 - I*(-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6)/(6*x) - (-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(7)/6)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_121():
    f = exp(I*atan(x)/3)/x**4
    F = -19*I*log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) - (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/324 + 19*I*log(1 + (I*x + 1)**(sympy.S(1)/3)/(-I*x + 1)**(sympy.S(1)/3) + (I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/324 - 19*sqrt(3)*I*atan(sqrt(3)*(1 - 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3)/162 + 19*sqrt(3)*I*atan(sqrt(3)*(1 + 2*(I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/3)/162 + 19*I*atanh((I*x + 1)**(sympy.S(1)/6)/(-I*x + 1)**(sympy.S(1)/6))/81 + 11*(-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6)/(27*x) - 7*I*(-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6)/(18*x**2) - (-I*x + 1)**(sympy.S(5)/6)*(I*x + 1)**(sympy.S(1)/6)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_122():
    f = x**2*exp(2*I*atan(x)/3)
    F = x*(-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(4)/3)/3 - I*(-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(4)/3)/9 - 11*I*(-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(1)/3)/27 + 11*I*log(I*x + 1)/81 + 11*I*log((-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/27 - 22*sqrt(3)*I*atan(2*sqrt(3)*(-I*x + 1)**(sympy.S(1)/3)/(3*(I*x + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/81
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_123():
    f = x*exp(2*I*atan(x)/3)
    F = (-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(4)/3)/2 + (-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(1)/3)/3 - log(I*x + 1)/9 - log((-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/3 + 2*sqrt(3)*atan(2*sqrt(3)*(-I*x + 1)**(sympy.S(1)/3)/(3*(I*x + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/9
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_124():
    f = exp(2*I*atan(x)/3)
    F = I*(-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(1)/3) - I*log(I*x + 1)/3 - I*log((-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1) + 2*sqrt(3)*I*atan(2*sqrt(3)*(-I*x + 1)**(sympy.S(1)/3)/(3*(I*x + 1)**(sympy.S(1)/3)) - sqrt(3)/3)/3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_125():
    f = exp(2*I*atan(x)/3)/x
    F = -log(x)/2 + log(I*x + 1)/2 + 3*log((-I*x + 1)**(sympy.S(1)/3)/(I*x + 1)**(sympy.S(1)/3) + 1)/2 + 3*log((-I*x + 1)**(sympy.S(1)/3) - (I*x + 1)**(sympy.S(1)/3))/2 - sqrt(3)*atan(2*sqrt(3)*(-I*x + 1)**(sympy.S(1)/3)/(3*(I*x + 1)**(sympy.S(1)/3)) - sqrt(3)/3) + sqrt(3)*atan(2*sqrt(3)*(-I*x + 1)**(sympy.S(1)/3)/(3*(I*x + 1)**(sympy.S(1)/3)) + sqrt(3)/3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_126():
    f = exp(2*I*atan(x)/3)/x**2
    F = -I*log(x)/3 + I*log((-I*x + 1)**(sympy.S(1)/3) - (I*x + 1)**(sympy.S(1)/3)) + 2*sqrt(3)*I*atan(2*sqrt(3)*(-I*x + 1)**(sympy.S(1)/3)/(3*(I*x + 1)**(sympy.S(1)/3)) + sqrt(3)/3)/3 - (-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(1)/3)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_127():
    f = exp(2*I*atan(x)/3)/x**3
    F = log(x)/9 - log((-I*x + 1)**(sympy.S(1)/3) - (I*x + 1)**(sympy.S(1)/3))/3 - 2*sqrt(3)*atan(2*sqrt(3)*(-I*x + 1)**(sympy.S(1)/3)/(3*(I*x + 1)**(sympy.S(1)/3)) + sqrt(3)/3)/9 - I*(-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(1)/3)/(3*x) - (-I*x + 1)**(sympy.S(2)/3)*(I*x + 1)**(sympy.S(4)/3)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_128():
    f = x**2*exp(I*atan(a*x)/4)
    F = x*(-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(9)/8)/(3*a**2) - I*(-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(9)/8)/(24*a**3) - 11*I*(-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(1)/8)/(32*a**3) - 11*I*sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) + 11*I*sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) - 11*I*sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) + 11*I*sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(256*a**3) + 11*I*sqrt(2 - sqrt(2))*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(128*a**3) - 11*I*sqrt(2 - sqrt(2))*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(128*a**3) + 11*I*sqrt(sqrt(2) + 2)*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(128*a**3) - 11*I*sqrt(sqrt(2) + 2)*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(128*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_129():
    f = x*exp(I*atan(a*x)/4)
    F = (-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(9)/8)/(2*a**2) + (-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(1)/8)/(8*a**2) + sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) - sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) + sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) - sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(64*a**2) - sqrt(2 - sqrt(2))*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(32*a**2) + sqrt(2 - sqrt(2))*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(32*a**2) - sqrt(sqrt(2) + 2)*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(32*a**2) + sqrt(sqrt(2) + 2)*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(32*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_130():
    f = exp(I*atan(a*x)/4)
    F = I*(-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(1)/8)/a + I*sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) - I*sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) + I*sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) - I*sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/(8*a) - I*sqrt(2 - sqrt(2))*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*a) + I*sqrt(2 - sqrt(2))*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2)))/(4*a) - I*sqrt(sqrt(2) + 2)*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*a) + I*sqrt(sqrt(2) + 2)*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2))/(4*a)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_131():
    f = exp(I*atan(a*x)/4)/x
    F = sqrt(2)*log(1 + (I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4) - sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/2 - sqrt(2)*log(1 + (I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4) + sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/2 - sqrt(2 - sqrt(2))*log(-sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/2 + sqrt(2 - sqrt(2))*log(sqrt(2 - sqrt(2))*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/2 - sqrt(sqrt(2) + 2)*log(-sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/2 + sqrt(sqrt(2) + 2)*log(sqrt(sqrt(2) + 2)*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + (-I*a*x + 1)**(sympy.S(1)/4)/(I*a*x + 1)**(sympy.S(1)/4) + 1)/2 + sqrt(2 - sqrt(2))*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) - sqrt(2 - sqrt(2))*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(sqrt(2) + 2))/sqrt(2 - sqrt(2))) + sqrt(sqrt(2) + 2)*atan((-2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2)) - sqrt(sqrt(2) + 2)*atan((2*(-I*a*x + 1)**(sympy.S(1)/8)/(I*a*x + 1)**(sympy.S(1)/8) + sqrt(2 - sqrt(2)))/sqrt(sqrt(2) + 2)) - 2*atan((I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8)) + sqrt(2)*atan(1 - sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8)) - sqrt(2)*atan(1 + sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8)) - 2*atanh((I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_132():
    f = exp(I*atan(a*x)/4)/x**2
    F = sqrt(2)*I*a*log(1 + (I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4) - sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/8 - sqrt(2)*I*a*log(1 + (I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4) + sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/8 - I*a*atan((I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/2 + sqrt(2)*I*a*atan(1 - sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/4 - sqrt(2)*I*a*atan(1 + sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/4 - I*a*atanh((I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/2 - (-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(1)/8)/x
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_133():
    f = exp(I*atan(a*x)/4)/x**3
    F = -sqrt(2)*a**2*log(1 + (I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4) - sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/64 + sqrt(2)*a**2*log(1 + (I*a*x + 1)**(sympy.S(1)/4)/(-I*a*x + 1)**(sympy.S(1)/4) + sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/64 + a**2*atan((I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/16 - sqrt(2)*a**2*atan(1 - sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/32 + sqrt(2)*a**2*atan(1 + sqrt(2)*(I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/32 + a**2*atanh((I*a*x + 1)**(sympy.S(1)/8)/(-I*a*x + 1)**(sympy.S(1)/8))/16 - I*a*(-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(1)/8)/(8*x) - (-I*a*x + 1)**(sympy.S(7)/8)*(I*a*x + 1)**(sympy.S(9)/8)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_134():
    f = x**m*exp(6*I*atan(a*x))
    F = x**(m + 1)*(4*m**2 + 8*m + 6)*hyper((1, m + 1), (m + 2,), I*a*x)/(m + 1) - x**(m + 1)*(I*a*x + 1)**2/((m + 1)*(-I*a*x + 1)**2) + 4*I*x**(m + 1)*(a*x*(m**2 + 3*m + 3) + I*(m + 1)**2)/((m + 1)*(-I*a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_135():
    f = x**m*exp(4*I*atan(a*x))
    F = -4*x**(m + 1)*hyper((1, m + 1), (m + 2,), I*a*x) + 4*x**(m + 1)/(-I*a*x + 1) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_136():
    f = x**m*exp(2*I*atan(a*x))
    F = 2*x**(m + 1)*hyper((1, m + 1), (m + 2,), I*a*x)/(m + 1) - x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_137():
    f = x**m*exp(-2*I*atan(a*x))
    F = 2*x**(m + 1)*hyper((1, m + 1), (m + 2,), -I*a*x)/(m + 1) - x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_138():
    f = x**m*exp(-4*I*atan(a*x))
    F = -4*x**(m + 1)*hyper((1, m + 1), (m + 2,), -I*a*x) + 4*x**(m + 1)/(I*a*x + 1) + x**(m + 1)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_139():
    f = x**m*exp(-6*I*atan(a*x))
    F = -x**(m + 1)*(-I*a*x + 1)**2/((m + 1)*(I*a*x + 1)**2) + x**(m + 1)*(4*m**2 + 8*m + 6)*hyper((1, m + 1), (m + 2,), -I*a*x)/(m + 1) + 4*I*x**(m + 1)*(-a*x*(m**2 + 3*m + 3) + I*(m + 1)**2)/((m + 1)*(I*a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_140():
    f = x**m*exp(3*I*atan(a*x))
    F = -I*a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), -a**2*x**2)/(m + 2) + 4*I*a*x**(m + 2)*hyper((sympy.S(3)/2, m/2 + 1), (m/2 + 2,), -a**2*x**2)/(m + 2) - 3*x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -a**2*x**2)/(m + 1) + 4*x**(m + 1)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_141():
    f = x**m*exp(I*atan(a*x))
    F = I*a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), -a**2*x**2)/(m + 2) + x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_142():
    f = x**m*exp(-I*atan(a*x))
    F = -I*a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), -a**2*x**2)/(m + 2) + x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_143():
    f = x**m*exp(-3*I*atan(a*x))
    F = I*a*x**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), -a**2*x**2)/(m + 2) - 4*I*a*x**(m + 2)*hyper((sympy.S(3)/2, m/2 + 1), (m/2 + 2,), -a**2*x**2)/(m + 2) - 3*x**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -a**2*x**2)/(m + 1) + 4*x**(m + 1)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -a**2*x**2)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_144():
    f = x**m*exp(5*I*atan(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-5)/4, sympy.S(5)/4, m + 2, -I*a*x, I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_145():
    f = x**m*exp(3*I*atan(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-3)/4, sympy.S(3)/4, m + 2, -I*a*x, I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_146():
    f = x**m*exp(I*atan(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/4, sympy.S(1)/4, m + 2, -I*a*x, I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_147():
    f = x**m*exp(-I*atan(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/4, sympy.S(1)/4, m + 2, I*a*x, -I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_148():
    f = x**m*exp(-3*I*atan(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-3)/4, sympy.S(3)/4, m + 2, I*a*x, -I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_149():
    f = x**m*exp(-5*I*atan(a*x)/2)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-5)/4, sympy.S(5)/4, m + 2, I*a*x, -I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_150():
    f = x**m*exp(2*atan(x)/3)
    F = x**(m + 1)*appellf1(m + 1, -I/3, I/3, m + 2, I*x, -I*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_151():
    f = x**m*exp(atan(x)/3)
    F = x**(m + 1)*appellf1(m + 1, -I/6, I/6, m + 2, I*x, -I*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_152():
    f = x**m*exp(I*atan(a*x)/4)
    F = x**(m + 1)*appellf1(m + 1, sympy.S(-1)/8, sympy.S(1)/8, m + 2, -I*a*x, I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_153():
    f = x**m*exp(I*n*atan(a*x))
    F = x**(m + 1)*appellf1(m + 1, -n/2, n/2, m + 2, -I*a*x, I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_154():
    f = x**3*exp(I*n*atan(a*x))
    F = -2**(n/2 - 2)*n*(n**2 + 8)*(-I*a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -I*a*x/2 + sympy.S.Half)/(3*a**4*(2 - n)) + x**2*(-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)/(4*a**2) - (-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)*(2*I*a*n*x + n**2 + 6)/(24*a**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_155():
    f = x**2*exp(I*n*atan(a*x))
    F = -2**(n/2)*I*(n**2 + 2)*(-I*a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -I*a*x/2 + sympy.S.Half)/(3*a**3*(2 - n)) + x*(-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)/(3*a**2) - I*n*(-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)/(6*a**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_156():
    f = x*exp(I*n*atan(a*x))
    F = 2**(n/2)*n*(-I*a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -I*a*x/2 + sympy.S.Half)/(a**2*(2 - n)) + (-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)/(2*a**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_157():
    f = exp(I*n*atan(a*x))
    F = 2**(n/2 + 1)*I*(-I*a*x + 1)**(1 - n/2)*hyper((-n/2, 1 - n/2), (2 - n/2,), -I*a*x/2 + sympy.S.Half)/(a*(2 - n))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_158():
    f = exp(I*n*atan(a*x))/x
    F = -2**(n/2 + 1)*hyper((-n/2, -n/2), (1 - n/2,), -I*a*x/2 + sympy.S.Half)/(n*(-I*a*x + 1)**(n/2)) + 2*(I*a*x + 1)**(n/2)*hyper((1, -n/2), (1 - n/2,), (-I*a*x + 1)/(I*a*x + 1))/(n*(-I*a*x + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_159():
    f = exp(I*n*atan(a*x))/x**2
    F = -4*I*a*(-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (-I*a*x + 1)/(I*a*x + 1))/(2 - n)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_160():
    f = exp(I*n*atan(a*x))/x**3
    F = 2*a**2*n*(-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (-I*a*x + 1)/(I*a*x + 1))/(2 - n) - (-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)/(2*x**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_161():
    f = exp(I*n*atan(a*x))/x**4
    F = 2*I*a**3*(n**2 + 2)*(-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 - 1)*hyper((2, 1 - n/2), (2 - n/2,), (-I*a*x + 1)/(I*a*x + 1))/(6 - 3*n) - I*a*n*(-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)/(6*x**2) - (-I*a*x + 1)**(1 - n/2)*(I*a*x + 1)**(n/2 + 1)/(3*x**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_162():
    f = x**4*exp(I*atan(a + b*x))
    F = x**3*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(5*b**2) - x**2*(8*a + I)*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(20*b**3) + sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)*(-96*a**3 - 86*I*a**2 + 114*a - b*x*(-72*a**2 - 28*I*a + 26) + 19*I)/(120*b**5) + sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(8*I*a**4 - 16*a**3 - 24*I*a**2 + 12*a + 3*I)/(8*b**5) + (8*a**4 + 16*I*a**3 - 24*a**2 - 12*I*a + 3)*asinh(a + b*x)/(8*b**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_163():
    f = x**3*exp(I*atan(a + b*x))
    F = x**2*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(4*b**2) - sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)*(-18*a**2 - 10*I*a + b*x*(12*a + 2*I) + 7)/(24*b**4) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(8*I*a**3 - 12*a**2 - 12*I*a + 3)/(8*b**4) + (-8*a**3 - 12*I*a**2 + 12*a + 3*I)*asinh(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_164():
    f = x**2*exp(I*atan(a + b*x))
    F = x*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(3*b**2) - (4*a + I)*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(6*b**3) - (-2*a**2 - 2*I*a + 1)*asinh(a + b*x)/(2*b**3) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(-2*I*a**2 + 2*a + I)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_165():
    f = x*exp(I*atan(a + b*x))
    F = -(2*a + I)*asinh(a + b*x)/(2*b**2) + (-2*I*a + 1)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(2*b**2) + sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_166():
    f = exp(I*atan(a + b*x))
    F = I*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/b + asinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_167():
    f = exp(I*atan(a + b*x))/x
    F = -2*sqrt(-a + I)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/sqrt(a + I) + I*asinh(a + b*x)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_168():
    f = exp(I*atan(a + b*x))/x**2
    F = 2*I*b*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/(sqrt(-a + I)*(a + I)**(sympy.S(3)/2)) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x*(-I*a + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_169():
    f = exp(I*atan(a + b*x))/x**3
    F = b**2*(2*I*a + 1)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(3)/2)*(a + I)**(sympy.S(5)/2)) - b*(2*I*a + 1)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x*(-2*a + 2*I)*(a + I)**2) - sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(x**2*(2*a**2 + 2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_170():
    f = exp(I*atan(a + b*x))/x**4
    F = b**3*(2*a - I*(1 - 2*a**2))*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(5)/2)*(a + I)**(sympy.S(7)/2)) + b**2*(-2*a**2 + 9*I*a + 4)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x*(a**2 + 1)**2*(-6*I*a + 6)) - b*(-2*a + 3*I)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x**2*(a**2 + 1)*(-6*I*a + 6)) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x**3*(-3*I*a + 3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_171():
    f = x**4*exp(2*I*atan(a + b*x))
    F = -x**5/5 + I*x**4/(2*b) + x**3*(-2*I*a + 2)/(3*b**2) + I*x**2*(a + I)**2/b**3 - 2*x*(-I*a + 1)**3/b**4 + 2*I*(a + I)**4*log(a + b*x + I)/b**5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_172():
    f = x**3*exp(2*I*atan(a + b*x))
    F = -x**4/4 + 2*I*x**3/(3*b) + x**2*(-I*a + 1)/b**2 + 2*I*x*(a + I)**2/b**3 - 2*(-I*a + 1)**3*log(a + b*x + I)/b**4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_173():
    f = x**2*exp(2*I*atan(a + b*x))
    F = -x**3/3 + I*x**2/b + x*(-2*I*a + 2)/b**2 + 2*I*(a + I)**2*log(a + b*x + I)/b**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_174():
    f = x*exp(2*I*atan(a + b*x))
    F = -x**2/2 + 2*I*x/b + (-2*I*a + 2)*log(a + b*x + I)/b**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_175():
    f = exp(2*I*atan(a + b*x))
    F = -x + 2*I*log(a + b*x + I)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_176():
    f = exp(2*I*atan(a + b*x))/x
    F = (-a + I)*log(x)/(a + I) - 2*log(a + b*x + I)/(-I*a + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_177():
    f = exp(2*I*atan(a + b*x))/x**2
    F = -2*I*b*log(x)/(a + I)**2 + 2*I*b*log(a + b*x + I)/(a + I)**2 - (-a + I)/(x*(a + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_178():
    f = exp(2*I*atan(a + b*x))/x**3
    F = -2*b**2*log(x)/(-I*a + 1)**3 + 2*b**2*log(a + b*x + I)/(-I*a + 1)**3 + 2*I*b/(x*(a + I)**2) - (-a + I)/(x**2*(2*a + 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_179():
    f = exp(2*I*atan(a + b*x))/x**4
    F = -2*I*b**3*log(x)/(a + I)**4 + 2*I*b**3*log(a + b*x + I)/(a + I)**4 + 2*b**2/(x*(-I*a + 1)**3) + I*b/(x**2*(a + I)**2) - (-a + I)/(x**3*(3*a + 3*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_180():
    f = x**4*exp(3*I*atan(a + b*x))
    F = -2*I*x**4*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(b*sqrt(-I*a - I*b*x + 1)) - 11*x**3*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(5*b**2) + x**2*(48*a + 51*I)*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(20*b**3) - I*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)*(112*I*a**3 - 422*a**2 - 458*I*a + b*x*(-104*I*a**2 + 236*a + 122*I) + 163)/(40*b**5) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(24*I*a**4 - 144*a**3 - 264*I*a**2 + 204*a + 57*I)/(8*b**5) - (24*a**4 + 144*I*a**3 - 264*a**2 - 204*I*a + 57)*asinh(a + b*x)/(8*b**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_181():
    f = x**3*exp(3*I*atan(a + b*x))
    F = -2*I*x**3*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(b*sqrt(-I*a - I*b*x + 1)) - 9*x**2*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(4*b**2) - I*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)*(-22*I*a**2 + 54*a - b*x*(-20*I*a + 22) + 29*I)/(8*b**4) + sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(24*I*a**3 - 108*a**2 - 132*I*a + 51)/(8*b**4) - (-24*a**3 - 108*I*a**2 + 132*a + 51*I)*asinh(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_182():
    f = x**2*exp(3*I*atan(a + b*x))
    F = -I*(a + I)**2*(I*a + I*b*x + 1)**(sympy.S(5)/2)/(b**3*sqrt(-I*a - I*b*x + 1)) + (-6*a**2 - 18*I*a + 11)*asinh(a + b*x)/(2*b**3) + I*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(5)/2)/(3*b**3) + sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)*(-6*I*a**2 + 18*a + 11*I)/(6*b**3) + sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(-6*I*a**2 + 18*a + 11*I)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_183():
    f = x*exp(3*I*atan(a + b*x))
    F = (6*a + 9*I)*asinh(a + b*x)/(2*b**2) - (-6*I*a + 9)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(2*b**2) - (-2*I*a + 3)*sqrt(-I*a - I*b*x + 1)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(2*b**2) - (-I*a + 1)*(I*a + I*b*x + 1)**(sympy.S(5)/2)/(b**2*sqrt(-I*a - I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_184():
    f = exp(3*I*atan(a + b*x))
    F = -3*I*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/b - 3*asinh(a + b*x)/b - 2*I*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(b*sqrt(-I*a - I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_185():
    f = exp(3*I*atan(a + b*x))/x
    F = -2*(-a + I)**(sympy.S(3)/2)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/(a + I)**(sympy.S(3)/2) - I*asinh(a + b*x) + 4*sqrt(I*a + I*b*x + 1)/((-I*a + 1)*sqrt(-I*a - I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_186():
    f = exp(3*I*atan(a + b*x))/x**2
    F = 6*I*b*sqrt(-a + I)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/(a + I)**(sympy.S(5)/2) - 6*I*b*sqrt(I*a + I*b*x + 1)/((a + I)**2*sqrt(-I*a - I*b*x + 1)) - (I*a + I*b*x + 1)**(sympy.S(3)/2)/(x*(-I*a + 1)*sqrt(-I*a - I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_187():
    f = exp(3*I*atan(a + b*x))/x**3
    F = b**2*(-6*a + 9*I)*sqrt(I*a + I*b*x + 1)/((a + I)**3*(I*a + 1)*sqrt(-I*a - I*b*x + 1)) + b**2*(6*I*a + 9)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/(sqrt(-a + I)*(a + I)**(sympy.S(7)/2)) + b*(-2*a + 3*I)*(I*a + I*b*x + 1)**(sympy.S(3)/2)/(x*(a + I)**2*(2*I*a + 2)*sqrt(-I*a - I*b*x + 1)) - (I*a + I*b*x + 1)**(sympy.S(5)/2)/(x**2*(2*a**2 + 2)*sqrt(-I*a - I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_188():
    f = exp(3*I*atan(a + b*x))/x**4
    F = -b**3*(-6*I*a**2 - 18*a + 11*I)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(3)/2)*(a + I)**(sympy.S(9)/2)) + b**3*(-2*a**2 + 51*I*a + 52)*sqrt(I*a + I*b*x + 1)/((-6*a + 6*I)*(a + I)**4*sqrt(-I*a - I*b*x + 1)) + b**2*(16*I*a + 19)*sqrt(I*a + I*b*x + 1)/(x*(-6*a + 6*I)*(a + I)**3*sqrt(-I*a - I*b*x + 1)) + 7*I*b*sqrt(I*a + I*b*x + 1)/(6*x**2*(a + I)**2*sqrt(-I*a - I*b*x + 1)) - (-a + I)*sqrt(I*a + I*b*x + 1)/(x**3*(3*a + 3*I)*sqrt(-I*a - I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_189():
    f = x**4*exp(-I*atan(a + b*x))
    F = x**3*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(5*b**2) + x**2*(-8*a + I)*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(20*b**3) - (-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)*(96*a**3 - 86*I*a**2 - 114*a + b*x*(-72*a**2 + 28*I*a + 26) + 19*I)/(120*b**5) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(8*I*a**4 + 16*a**3 - 24*I*a**2 - 12*a + 3*I)/(8*b**5) + (8*a**4 - 16*I*a**3 - 24*a**2 + 12*I*a + 3)*asinh(a + b*x)/(8*b**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_190():
    f = x**3*exp(-I*atan(a + b*x))
    F = x**2*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(4*b**2) - (-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)*(-18*a**2 + 10*I*a - b*x*(-12*a + 2*I) + 7)/(24*b**4) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(-8*I*a**3 - 12*a**2 + 12*I*a + 3)/(8*b**4) - (8*a**3 - 12*I*a**2 - 12*a + 3*I)*asinh(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_191():
    f = x**2*exp(-I*atan(a + b*x))
    F = x*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(3*b**2) + (-4*a + I)*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(6*b**3) - (-2*a**2 + 2*I*a + 1)*asinh(a + b*x)/(2*b**3) + sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(-2*I*a**2 - 2*a + I)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_192():
    f = x*exp(-I*atan(a + b*x))
    F = (-2*a + I)*asinh(a + b*x)/(2*b**2) + (2*I*a + 1)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(2*b**2) + (-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_193():
    f = exp(-I*atan(a + b*x))
    F = -I*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/b + asinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_194():
    f = exp(-I*atan(a + b*x))/x
    F = -I*asinh(a + b*x) - 2*sqrt(a + I)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/sqrt(-a + I)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_195():
    f = exp(-I*atan(a + b*x))/x**2
    F = -2*I*b*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(3)/2)*sqrt(a + I)) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x*(I*a + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_196():
    f = exp(-I*atan(a + b*x))/x**3
    F = b**2*(-2*I*a + 1)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(5)/2)*(a + I)**(sympy.S(3)/2)) + b*(-2*I*a + 1)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(2*x*(-a + I)**2*(a + I)) - (-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(x**2*(2*a**2 + 2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_197():
    f = exp(-I*atan(a + b*x))/x**4
    F = b**3*(2*a + I*(1 - 2*a**2))*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(7)/2)*(a + I)**(sympy.S(5)/2)) + b**2*(-2*a**2 - 9*I*a + 4)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x*(a**2 + 1)**2*(6*I*a + 6)) + b*(-2*I*a + 3)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(6*x**2*(-a + I)**2*(a + I)) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(x**3*(3*I*a + 3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_198():
    f = x**4*exp(-2*I*atan(a + b*x))
    F = -x**5/5 - I*x**4/(2*b) + x**3*(2*I*a + 2)/(3*b**2) - I*x**2*(-a + I)**2/b**3 - 2*x*(I*a + 1)**3/b**4 - 2*I*(-a + I)**4*log(-a - b*x + I)/b**5
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_199():
    f = x**3*exp(-2*I*atan(a + b*x))
    F = -x**4/4 - 2*I*x**3/(3*b) + x**2*(I*a + 1)/b**2 - 2*I*x*(-a + I)**2/b**3 - 2*(I*a + 1)**3*log(-a - b*x + I)/b**4
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_200():
    f = x**2*exp(-2*I*atan(a + b*x))
    F = -x**3/3 - I*x**2/b + x*(2*I*a + 2)/b**2 - 2*I*(-a + I)**2*log(-a - b*x + I)/b**3
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_201():
    f = x*exp(-2*I*atan(a + b*x))
    F = -x**2/2 - 2*I*x/b + (2*I*a + 2)*log(-a - b*x + I)/b**2
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_202():
    f = exp(-2*I*atan(a + b*x))
    F = -x - 2*I*log(-a - b*x + I)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_203():
    f = exp(-2*I*atan(a + b*x))/x
    F = -2*log(-a - b*x + I)/(I*a + 1) + (a + I)*log(x)/(-a + I)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_204():
    f = exp(-2*I*atan(a + b*x))/x**2
    F = 2*I*b*log(x)/(-a + I)**2 - 2*I*b*log(-a - b*x + I)/(-a + I)**2 - (a + I)/(x*(-a + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_205():
    f = exp(-2*I*atan(a + b*x))/x**3
    F = -2*b**2*log(x)/(I*a + 1)**3 + 2*b**2*log(-a - b*x + I)/(I*a + 1)**3 - 2*I*b/(x*(-a + I)**2) + (-a - I)/(x**2*(-2*a + 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_206():
    f = exp(-2*I*atan(a + b*x))/x**4
    F = 2*I*b**3*log(x)/(-a + I)**4 - 2*I*b**3*log(-a - b*x + I)/(-a + I)**4 + 2*b**2/(x*(I*a + 1)**3) - I*b/(x**2*(-a + I)**2) + (-a - I)/(x**3*(-3*a + 3*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_207():
    f = x**4*exp(-3*I*atan(a + b*x))
    F = 2*I*x**4*(-I*a - I*b*x + 1)**(sympy.S(3)/2)/(b*sqrt(I*a + I*b*x + 1)) - 11*x**3*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(5*b**2) - x**2*(-48*a + 51*I)*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(20*b**3) + I*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)*(-112*I*a**3 - 422*a**2 + 458*I*a - b*x*(-104*I*a**2 - 236*a + 122*I) + 163)/(40*b**5) + sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(24*I*a**4 + 144*a**3 - 264*I*a**2 - 204*a + 57*I)/(8*b**5) - (24*a**4 - 144*I*a**3 - 264*a**2 + 204*I*a + 57)*asinh(a + b*x)/(8*b**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_208():
    f = x**3*exp(-3*I*atan(a + b*x))
    F = 2*I*x**3*(-I*a - I*b*x + 1)**(sympy.S(3)/2)/(b*sqrt(I*a + I*b*x + 1)) - 9*x**2*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(4*b**2) - I*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)*(-22*I*a**2 - 54*a + b*x*(20*I*a + 22) + 29*I)/(8*b**4) + sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(-24*I*a**3 - 108*a**2 + 132*I*a + 51)/(8*b**4) + (24*a**3 - 108*I*a**2 - 132*a + 51*I)*asinh(a + b*x)/(8*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_209():
    f = x**2*exp(-3*I*atan(a + b*x))
    F = I*(-a + I)**2*(-I*a - I*b*x + 1)**(sympy.S(5)/2)/(b**3*sqrt(I*a + I*b*x + 1)) + (-6*a**2 + 18*I*a + 11)*asinh(a + b*x)/(2*b**3) - I*(-I*a - I*b*x + 1)**(sympy.S(5)/2)*sqrt(I*a + I*b*x + 1)/(3*b**3) - (-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)*(-6*I*a**2 - 18*a + 11*I)/(6*b**3) - sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)*(-6*I*a**2 - 18*a + 11*I)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_210():
    f = x*exp(-3*I*atan(a + b*x))
    F = -(-6*a + 9*I)*asinh(a + b*x)/(2*b**2) - (I*a + 1)*(-I*a - I*b*x + 1)**(sympy.S(5)/2)/(b**2*sqrt(I*a + I*b*x + 1)) - (2*I*a + 3)*(-I*a - I*b*x + 1)**(sympy.S(3)/2)*sqrt(I*a + I*b*x + 1)/(2*b**2) - (6*I*a + 9)*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_211():
    f = exp(-3*I*atan(a + b*x))
    F = 2*I*(-I*a - I*b*x + 1)**(sympy.S(3)/2)/(b*sqrt(I*a + I*b*x + 1)) + 3*I*sqrt(-I*a - I*b*x + 1)*sqrt(I*a + I*b*x + 1)/b - 3*asinh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_212():
    f = exp(-3*I*atan(a + b*x))/x
    F = I*asinh(a + b*x) + 4*sqrt(-I*a - I*b*x + 1)/((I*a + 1)*sqrt(I*a + I*b*x + 1)) - 2*(a + I)**(sympy.S(3)/2)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/(-a + I)**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_213():
    f = exp(-3*I*atan(a + b*x))/x**2
    F = 6*I*b*sqrt(-I*a - I*b*x + 1)/((-a + I)**2*sqrt(I*a + I*b*x + 1)) - 6*I*b*sqrt(a + I)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/(-a + I)**(sympy.S(5)/2) - (-I*a - I*b*x + 1)**(sympy.S(3)/2)/(x*(I*a + 1)*sqrt(I*a + I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_214():
    f = exp(-3*I*atan(a + b*x))/x**3
    F = -b**2*(6*a + 9*I)*sqrt(-I*a - I*b*x + 1)/((a + I)*(I*a + 1)**3*sqrt(I*a + I*b*x + 1)) + b**2*(-6*I*a + 9)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(7)/2)*sqrt(a + I)) + b*(-2*I*a + 3)*(-I*a - I*b*x + 1)**(sympy.S(3)/2)/(2*x*(-a + I)**2*(a + I)*sqrt(I*a + I*b*x + 1)) - (-I*a - I*b*x + 1)**(sympy.S(5)/2)/(x**2*(2*a**2 + 2)*sqrt(I*a + I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_215():
    f = exp(-3*I*atan(a + b*x))/x**4
    F = -b**3*(-2*a**2 - 51*I*a + 52)*sqrt(-I*a - I*b*x + 1)/(6*(-a + I)**4*(a + I)*sqrt(I*a + I*b*x + 1)) + b**3*(-6*I*a**2 + 18*a + 11*I)*atanh(sqrt(a + I)*sqrt(I*a + I*b*x + 1)/(sqrt(-a + I)*sqrt(-I*a - I*b*x + 1)))/((-a + I)**(sympy.S(9)/2)*(a + I)**(sympy.S(3)/2)) + b**2*(-16*I*a + 19)*sqrt(-I*a - I*b*x + 1)/(6*x*(-a + I)**3*(a + I)*sqrt(I*a + I*b*x + 1)) - 7*I*b*sqrt(-I*a - I*b*x + 1)/(6*x**2*(-a + I)**2*sqrt(I*a + I*b*x + 1)) - (a + I)*sqrt(-I*a - I*b*x + 1)/(x**3*(-3*a + 3*I)*sqrt(I*a + I*b*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_216():
    f = x**2*exp(I*atan(a + b*x)/2)
    F = x*(-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(5)/4)/(3*b**2) - (8*a + I)*(-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(5)/4)/(12*b**3) - (-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)*(-8*I*a**2 + 4*a + 3*I)/(8*b**3) - sqrt(2)*(-8*I*a**2 + 4*a + 3*I)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) + sqrt(2)*(-8*I*a**2 + 4*a + 3*I)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) - sqrt(2)*(-8*I*a**2 + 4*a + 3*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(16*b**3) - sqrt(2)*(-8*I*a**2 + 4*a + 3*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(16*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_217():
    f = x*exp(I*atan(a + b*x)/2)
    F = (-4*I*a + 1)*(-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/(4*b**2) + sqrt(2)*(-4*I*a + 1)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) - sqrt(2)*(-4*I*a + 1)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) + sqrt(2)*(-4*I*a + 1)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(8*b**2) + sqrt(2)*(-4*I*a + 1)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(8*b**2) + (-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(5)/4)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_218():
    f = exp(I*atan(a + b*x)/2)
    F = I*(-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/b + sqrt(2)*I*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) - sqrt(2)*I*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) + sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(2*b) + sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_219():
    f = exp(I*atan(a + b*x)/2)/x
    F = -2*(-a + I)**(sympy.S(1)/4)*atan((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)))/(a + I)**(sympy.S(1)/4) - 2*(-a + I)**(sympy.S(1)/4)*atanh((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)))/(a + I)**(sympy.S(1)/4) - sqrt(2)*log(1 + sqrt(I*(a + b*x) + 1)/sqrt(-I*(a + b*x) + 1) - sqrt(2)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/(-I*(a + b*x) + 1)**(sympy.S(1)/4))/2 + sqrt(2)*log(1 + sqrt(I*(a + b*x) + 1)/sqrt(-I*(a + b*x) + 1) + sqrt(2)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/(-I*(a + b*x) + 1)**(sympy.S(1)/4))/2 - sqrt(2)*atan(1 - sqrt(2)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/(-I*(a + b*x) + 1)**(sympy.S(1)/4)) + sqrt(2)*atan(1 + sqrt(2)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/(-I*(a + b*x) + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_220():
    f = exp(I*atan(a + b*x)/2)/x**2
    F = I*b*atan((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(3)/4)*(a + I)**(sympy.S(5)/4)) + I*b*atanh((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(3)/4)*(a + I)**(sympy.S(5)/4)) - (I*(a + b*x) + 1)**(sympy.S(1)/4)*(a + b*x + I)/(x*(a + I)*(-I*(a + b*x) + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_221():
    f = x**2*exp(3*I*atan(a + b*x)/2)
    F = x*(-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(7)/4)/(3*b**2) - (8*a + 3*I)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(7)/4)/(12*b**3) - (-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)*(-24*I*a**2 + 36*a + 17*I)/(24*b**3) + sqrt(2)*(-24*I*a**2 + 36*a + 17*I)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) - sqrt(2)*(-24*I*a**2 + 36*a + 17*I)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) - sqrt(2)*(-24*I*a**2 + 36*a + 17*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(16*b**3) - sqrt(2)*(-24*I*a**2 + 36*a + 17*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(16*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_222():
    f = x*exp(3*I*atan(a + b*x)/2)
    F = -sqrt(2)*(-12*I*a + 9)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) + sqrt(2)*(-12*I*a + 9)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) + sqrt(2)*(-12*I*a + 9)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(8*b**2) + sqrt(2)*(-12*I*a + 9)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(8*b**2) + (-4*I*a + 3)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/(4*b**2) + (-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(7)/4)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_223():
    f = exp(3*I*atan(a + b*x)/2)
    F = I*(-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/b - 3*sqrt(2)*I*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) + 3*sqrt(2)*I*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(2*b) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_224():
    f = exp(3*I*atan(a + b*x)/2)/x
    F = 2*(-a + I)**(sympy.S(3)/4)*atan((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/(a + I)**(sympy.S(3)/4) - 2*(-a + I)**(sympy.S(3)/4)*atanh((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/(a + I)**(sympy.S(3)/4) + sqrt(2)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/2 - sqrt(2)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1) - sqrt(2)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_225():
    f = exp(3*I*atan(a + b*x)/2)/x**2
    F = -3*I*b*atan((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(1)/4)*(a + I)**(sympy.S(7)/4)) + 3*I*b*atanh((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(1)/4)*(a + I)**(sympy.S(7)/4)) - (-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/(x*(-I*a + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_226():
    f = x**2*exp(-I*atan(a + b*x)/2)
    F = x*(-I*a - I*b*x + 1)**(sympy.S(5)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/(3*b**2) + (-8*a + I)*(-I*a - I*b*x + 1)**(sympy.S(5)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/(12*b**3) + (-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)*(-8*I*a**2 - 4*a + 3*I)/(8*b**3) + sqrt(2)*(-8*I*a**2 - 4*a + 3*I)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) - sqrt(2)*(-8*I*a**2 - 4*a + 3*I)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) - sqrt(2)*(-8*I*a**2 - 4*a + 3*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(16*b**3) - sqrt(2)*(-8*I*a**2 - 4*a + 3*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(16*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_227():
    f = x*exp(-I*atan(a + b*x)/2)
    F = (4*I*a + 1)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/(4*b**2) + sqrt(2)*(4*I*a + 1)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) - sqrt(2)*(4*I*a + 1)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) - sqrt(2)*(4*I*a + 1)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(8*b**2) - sqrt(2)*(4*I*a + 1)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(8*b**2) + (-I*a - I*b*x + 1)**(sympy.S(5)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_228():
    f = exp(-I*atan(a + b*x)/2)
    F = -I*(-I*a - I*b*x + 1)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(3)/4)/b - sqrt(2)*I*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) + sqrt(2)*I*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) + sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(2*b) + sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_229():
    f = exp(-I*atan(a + b*x)/2)/x
    F = -sqrt(2)*log(-sqrt(2)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/(I*(a + b*x) + 1)**(sympy.S(1)/4) + sqrt(-I*(a + b*x) + 1)/sqrt(I*(a + b*x) + 1) + 1)/2 + sqrt(2)*log(sqrt(2)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/(I*(a + b*x) + 1)**(sympy.S(1)/4) + sqrt(-I*(a + b*x) + 1)/sqrt(I*(a + b*x) + 1) + 1)/2 + sqrt(2)*atan(sqrt(2)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/(I*(a + b*x) + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/(I*(a + b*x) + 1)**(sympy.S(1)/4) + 1) - 2*(a + I)**(sympy.S(1)/4)*atan((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)))/(-a + I)**(sympy.S(1)/4) - 2*(a + I)**(sympy.S(1)/4)*atanh((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)))/(-a + I)**(sympy.S(1)/4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_230():
    f = exp(-I*atan(a + b*x)/2)/x**2
    F = -I*b*atan((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(5)/4)*(a + I)**(sympy.S(3)/4)) - I*b*atanh((-a + I)**(sympy.S(1)/4)*(-I*(a + b*x) + 1)**(sympy.S(1)/4)/((a + I)**(sympy.S(1)/4)*(I*(a + b*x) + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(5)/4)*(a + I)**(sympy.S(3)/4)) - (-I*(a + b*x) + 1)**(sympy.S(1)/4)*(-a - b*x + I)/(x*(-a + I)*(I*(a + b*x) + 1)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_231():
    f = x**2*exp(-3*I*atan(a + b*x)/2)
    F = x*(-I*a - I*b*x + 1)**(sympy.S(7)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/(3*b**2) + (-8*a + 3*I)*(-I*a - I*b*x + 1)**(sympy.S(7)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/(12*b**3) + (-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)*(-24*I*a**2 - 36*a + 17*I)/(24*b**3) - sqrt(2)*(-24*I*a**2 - 36*a + 17*I)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) + sqrt(2)*(-24*I*a**2 - 36*a + 17*I)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(32*b**3) - sqrt(2)*(-24*I*a**2 - 36*a + 17*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(16*b**3) - sqrt(2)*(-24*I*a**2 - 36*a + 17*I)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(16*b**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_232():
    f = x*exp(-3*I*atan(a + b*x)/2)
    F = (4*I*a + 3)*(-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/(4*b**2) - sqrt(2)*(12*I*a + 9)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) + sqrt(2)*(12*I*a + 9)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(16*b**2) - sqrt(2)*(12*I*a + 9)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(8*b**2) - sqrt(2)*(12*I*a + 9)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(8*b**2) + (-I*a - I*b*x + 1)**(sympy.S(7)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/(2*b**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_233():
    f = exp(-3*I*atan(a + b*x)/2)
    F = -I*(-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/b + 3*sqrt(2)*I*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) - 3*sqrt(2)*I*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/(4*b) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1)/(2*b) + 3*sqrt(2)*I*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1)/(2*b)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_234():
    f = exp(-3*I*atan(a + b*x)/2)/x
    F = sqrt(2)*log(-sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/2 - sqrt(2)*log(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + sqrt(-I*a - I*b*x + 1)/sqrt(I*a + I*b*x + 1) + 1)/2 + sqrt(2)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) - 1) + sqrt(2)*atan(sqrt(2)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)/(I*a + I*b*x + 1)**(sympy.S(1)/4) + 1) - 2*(a + I)**(sympy.S(3)/4)*atan((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/(-a + I)**(sympy.S(3)/4) - 2*(a + I)**(sympy.S(3)/4)*atanh((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/(-a + I)**(sympy.S(3)/4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_235():
    f = exp(-3*I*atan(a + b*x)/2)/x**2
    F = -3*I*b*atan((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(7)/4)*(a + I)**(sympy.S(1)/4)) - 3*I*b*atanh((a + I)**(sympy.S(1)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/((-a + I)**(sympy.S(1)/4)*(-I*a - I*b*x + 1)**(sympy.S(1)/4)))/((-a + I)**(sympy.S(7)/4)*(a + I)**(sympy.S(1)/4)) - (-I*a - I*b*x + 1)**(sympy.S(3)/4)*(I*a + I*b*x + 1)**(sympy.S(1)/4)/(x*(I*a + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_236():
    f = x**m*exp(n*atan(a + b*x))
    F = x**(m + 1)*(-b*x/(-a + I) + 1)**(I*n/2)*(-I*a - I*b*x + 1)**(I*n/2)*appellf1(m + 1, -I*n/2, I*n/2, m + 2, -b*x/(a + I), b*x/(-a + I))/((m + 1)*(b*x/(a + I) + 1)**(I*n/2)*(I*a + I*b*x + 1)**(I*n/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_237():
    f = x**3*exp(n*atan(a + b*x))
    F = 2**(-I*n/2 - 2)*(-I*a - I*b*x + 1)**(I*n/2 + 1)*(24*a**3 + 36*a**2*n - 12*a*(2 - n**2) - n*(8 - n**2))*hyper((I*n/2, I*n/2 + 1), (I*n/2 + 2,), -I*a/2 - I*b*x/2 + sympy.S.Half)/(3*b**4*(-n + 2*I)) + x**2*(-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 + 1)/(4*b**2) - (-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 + 1)*(-18*a**2 - 10*a*n + 2*b*x*(6*a + n) - n**2 + 6)/(24*b**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_238():
    f = x**2*exp(n*atan(a + b*x))
    F = x*(-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 + 1)/(3*b**2) - (4*a + n)*(-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 + 1)/(6*b**3) + (-I*a - I*b*x + 1)**(I*n/2 + 1)*(-6*a**2 - 6*a*n - n**2 + 2)*hyper((I*n/2, I*n/2 + 1), (I*n/2 + 2,), -I*a/2 - I*b*x/2 + sympy.S.Half)/(3*2**(I*n/2)*b**3*(-n + 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_239():
    f = x*exp(n*atan(a + b*x))
    F = (-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 + 1)/(2*b**2) + (2*a + n)*(-I*a - I*b*x + 1)**(I*n/2 + 1)*hyper((I*n/2, I*n/2 + 1), (I*n/2 + 2,), -I*a/2 - I*b*x/2 + sympy.S.Half)/(2**(I*n/2)*b**2*(-n + 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_240():
    f = exp(n*atan(a + b*x))
    F = -2**(-I*n/2 + 1)*(-I*a - I*b*x + 1)**(I*n/2 + 1)*hyper((I*n/2, I*n/2 + 1), (I*n/2 + 2,), -I*a/2 - I*b*x/2 + sympy.S.Half)/(b*(-n + 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_241():
    f = exp(n*atan(a + b*x))/x
    F = -2**(-I*n/2 + 1)*I*(-I*a - I*b*x + 1)**(I*n/2)*hyper((I*n/2, I*n/2), (I*n/2 + 1,), -I*a/2 - I*b*x/2 + sympy.S.Half)/n + 2*I*(-I*a - I*b*x + 1)**(I*n/2)*hyper((1, I*n/2), (I*n/2 + 1,), (-a + I)*(-I*a - I*b*x + 1)/((a + I)*(I*a + I*b*x + 1)))/(n*(I*a + I*b*x + 1)**(I*n/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_242():
    f = exp(n*atan(a + b*x))/x**2
    F = -4*b*(-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 - 1)*hyper((2, I*n/2 + 1), (I*n/2 + 2,), (-a + I)*(-I*a - I*b*x + 1)/((a + I)*(I*a + I*b*x + 1)))/((a + I)**2*(-n + 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_243():
    f = exp(n*atan(a + b*x))/x**3
    F = -2*b**2*(2*a - n)*(-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 - 1)*hyper((2, I*n/2 + 1), (I*n/2 + 2,), (-a + I)*(-I*a - I*b*x + 1)/((a + I)*(I*a + I*b*x + 1)))/((-a + I)*(a + I)**3*(-n + 2*I)) - (-I*a - I*b*x + 1)**(I*n/2 + 1)*(I*a + I*b*x + 1)**(-I*n/2 + 1)/(x**2*(2*a**2 + 2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_244():
    f = (a**2*c*x**2 + c)**p*exp(atan(a*x))
    F = 2**(p + 1 - I/2)*I*(-I*a*x + 1)**(p + 1 + I/2)*(a**2*c*x**2 + c)**p*hyper((p + 1 + I/2, -p + I/2), (p + 2 + I/2,), -I*a*x/2 + sympy.S.Half)/(a*(a**2*x**2 + 1)**p*(2*p + 2 + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_245():
    f = (a**2*c*x**2 + c)**2*exp(atan(a*x))
    F = 2**(3 - I/2)*c**2*(sympy.S(1)/37 + 6*I/37)*(-I*a*x + 1)**(3 + I/2)*hyper((-2 + I/2, 3 + I/2), (4 + I/2,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_246():
    f = (a**2*c*x**2 + c)*exp(atan(a*x))
    F = 2**(2 - I/2)*c*(sympy.S(1)/17 + 4*I/17)*(-I*a*x + 1)**(2 + I/2)*hyper((-1 + I/2, 2 + I/2), (3 + I/2,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_247():
    f = exp(atan(a*x))
    F = 2**(1 - I/2)*(sympy.S(1)/5 + 2*I/5)*(-I*a*x + 1)**(1 + I/2)*hyper((I/2, 1 + I/2), (2 + I/2,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_248():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)
    F = exp(atan(a*x))/(a*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_249():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)**2
    F = (2*a*x + 1)*exp(atan(a*x))/(5*a*c**2*(a**2*x**2 + 1)) + 2*exp(atan(a*x))/(5*a*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_250():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)**3
    F = 12*(2*a*x + 1)*exp(atan(a*x))/(85*a*c**3*(a**2*x**2 + 1)) + (4*a*x + 1)*exp(atan(a*x))/(17*a*c**3*(a**2*x**2 + 1)**2) + 24*exp(atan(a*x))/(85*a*c**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_251():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)**4
    F = 72*(2*a*x + 1)*exp(atan(a*x))/(629*a*c**4*(a**2*x**2 + 1)) + 30*(4*a*x + 1)*exp(atan(a*x))/(629*a*c**4*(a**2*x**2 + 1)**2) + (6*a*x + 1)*exp(atan(a*x))/(37*a*c**4*(a**2*x**2 + 1)**3) + 144*exp(atan(a*x))/(629*a*c**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_252():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)**5
    F = 4032*(2*a*x + 1)*exp(atan(a*x))/(40885*a*c**5*(a**2*x**2 + 1)) + 336*(4*a*x + 1)*exp(atan(a*x))/(8177*a*c**5*(a**2*x**2 + 1)**2) + 56*(6*a*x + 1)*exp(atan(a*x))/(2405*a*c**5*(a**2*x**2 + 1)**3) + (8*a*x + 1)*exp(atan(a*x))/(65*a*c**5*(a**2*x**2 + 1)**4) + 8064*exp(atan(a*x))/(40885*a*c**5)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_253():
    f = (a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(atan(a*x))
    F = 2**(sympy.S(3)/2 - I/2)*c*(sympy.S(1)/13 + 5*I/13)*(-I*a*x + 1)**(sympy.S(5)/2 + I/2)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-3)/2 + I/2, sympy.S(5)/2 + I/2), (sympy.S(7)/2 + I/2,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_254():
    f = sqrt(a**2*c*x**2 + c)*exp(atan(a*x))
    F = 2**(sympy.S.Half - I/2)*(sympy.S(1)/5 + 3*I/5)*(-I*a*x + 1)**(sympy.S(3)/2 + I/2)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-1)/2 + I/2, sympy.S(3)/2 + I/2), (sympy.S(5)/2 + I/2,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_255():
    f = exp(atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = 2**(sympy.S(-1)/2 - I/2)*(1 + I)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(sympy.S.Half + I/2)*hyper((sympy.S.Half + I/2, sympy.S.Half + I/2), (sympy.S(3)/2 + I/2,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_256():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = (a*x + 1)*exp(atan(a*x))/(2*a*c*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_257():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = (3*a*x + 1)*exp(atan(a*x))/(10*a*c*(a**2*c*x**2 + c)**(sympy.S(3)/2)) + 3*(a*x + 1)*exp(atan(a*x))/(10*a*c**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_258():
    f = exp(atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = (5*a*x + 1)*exp(atan(a*x))/(26*a*c*(a**2*c*x**2 + c)**(sympy.S(5)/2)) + (3*a*x + 1)*exp(atan(a*x))/(13*a*c**2*(a**2*c*x**2 + c)**(sympy.S(3)/2)) + 3*(a*x + 1)*exp(atan(a*x))/(13*a*c**3*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_259():
    f = (a**2*c*x**2 + c)**p*exp(2*atan(a*x))
    F = 2**(p - I)*I*(-I*a*x + 1)**(p + 1 + I)*(a**2*c*x**2 + c)**p*hyper((p + 1 + I, -p + I), (p + 2 + I,), -I*a*x/2 + sympy.S.Half)/(a*(a**2*x**2 + 1)**p*(p + 1 + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_260():
    f = (a**2*c*x**2 + c)**2*exp(2*atan(a*x))
    F = 2**(1 - I)*c**2*(sympy.S(1)/5 + 3*I/5)*(-I*a*x + 1)**(3 + I)*hyper((-2 + I, 3 + I), (4 + I,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_261():
    f = (a**2*c*x**2 + c)*exp(2*atan(a*x))
    F = 2**(1 - I)*c*(sympy.S(1)/5 + 2*I/5)*(-I*a*x + 1)**(2 + I)*hyper((-1 + I, 2 + I), (3 + I,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_262():
    f = exp(2*atan(a*x))
    F = 2**(-1 - I)*(1 + I)*(-I*a*x + 1)**(1 + I)*hyper((I, 1 + I), (2 + I,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_263():
    f = exp(2*atan(a*x))/(a**2*c*x**2 + c)
    F = exp(2*atan(a*x))/(2*a*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_264():
    f = exp(2*atan(a*x))/(a**2*c*x**2 + c)**2
    F = (a*x + 1)*exp(2*atan(a*x))/(4*a*c**2*(a**2*x**2 + 1)) + exp(2*atan(a*x))/(8*a*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_265():
    f = exp(2*atan(a*x))/(a**2*c*x**2 + c)**3
    F = 3*(a*x + 1)*exp(2*atan(a*x))/(20*a*c**3*(a**2*x**2 + 1)) + (2*a*x + 1)*exp(2*atan(a*x))/(10*a*c**3*(a**2*x**2 + 1)**2) + 3*exp(2*atan(a*x))/(40*a*c**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_266():
    f = exp(2*atan(a*x))/(a**2*c*x**2 + c)**4
    F = 9*(a*x + 1)*exp(2*atan(a*x))/(80*a*c**4*(a**2*x**2 + 1)) + 3*(2*a*x + 1)*exp(2*atan(a*x))/(40*a*c**4*(a**2*x**2 + 1)**2) + (3*a*x + 1)*exp(2*atan(a*x))/(20*a*c**4*(a**2*x**2 + 1)**3) + 9*exp(2*atan(a*x))/(160*a*c**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_267():
    f = (a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(2*atan(a*x))
    F = 2**(sympy.S(5)/2 - I)*c*(sympy.S(2)/29 + 5*I/29)*(-I*a*x + 1)**(sympy.S(5)/2 + I)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-3)/2 + I, sympy.S(5)/2 + I), (sympy.S(7)/2 + I,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_268():
    f = sqrt(a**2*c*x**2 + c)*exp(2*atan(a*x))
    F = 2**(sympy.S(3)/2 - I)*(sympy.S(2)/13 + 3*I/13)*(-I*a*x + 1)**(sympy.S(3)/2 + I)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-1)/2 + I, sympy.S(3)/2 + I), (sympy.S(5)/2 + I,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_269():
    f = exp(2*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = 2**(sympy.S.Half - I)*(sympy.S(2)/5 + I/5)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(sympy.S.Half + I)*hyper((sympy.S.Half + I, sympy.S.Half + I), (sympy.S(3)/2 + I,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_270():
    f = exp(2*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = (a*x + 2)*exp(2*atan(a*x))/(5*a*c*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_271():
    f = exp(2*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = (3*a*x + 2)*exp(2*atan(a*x))/(13*a*c*(a**2*c*x**2 + c)**(sympy.S(3)/2)) + 6*(a*x + 2)*exp(2*atan(a*x))/(65*a*c**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_272():
    f = exp(2*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = (5*a*x + 2)*exp(2*atan(a*x))/(29*a*c*(a**2*c*x**2 + c)**(sympy.S(5)/2)) + 20*(3*a*x + 2)*exp(2*atan(a*x))/(377*a*c**2*(a**2*c*x**2 + c)**(sympy.S(3)/2)) + 24*(a*x + 2)*exp(2*atan(a*x))/(377*a*c**3*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_273():
    f = (a**2*c*x**2 + c)**p*exp(-atan(a*x))
    F = 2**(p + 1 + I/2)*(-I*a*x + 1)**(p + 1 - I/2)*(a**2*c*x**2 + c)**p*hyper((p + 1 - I/2, -p - I/2), (p + 2 - I/2,), -I*a*x/2 + sympy.S.Half)/(a*(a**2*x**2 + 1)**p*(-2*I*p - 1 - 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_274():
    f = (a**2*c*x**2 + c)**2*exp(-atan(a*x))
    F = -2**(3 + I/2)*c**2*(sympy.S(1)/37 - 6*I/37)*(-I*a*x + 1)**(3 - I/2)*hyper((-2 - I/2, 3 - I/2), (4 - I/2,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_275():
    f = (a**2*c*x**2 + c)*exp(-atan(a*x))
    F = -2**(2 + I/2)*c*(sympy.S(1)/17 - 4*I/17)*(-I*a*x + 1)**(2 - I/2)*hyper((-1 - I/2, 2 - I/2), (3 - I/2,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_276():
    f = exp(-atan(a*x))
    F = -2**(1 + I/2)*(sympy.S(1)/5 - 2*I/5)*(-I*a*x + 1)**(1 - I/2)*hyper((-I/2, 1 - I/2), (2 - I/2,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_277():
    f = exp(-atan(a*x))/(a**2*c*x**2 + c)
    F = -exp(-atan(a*x))/(a*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_278():
    f = exp(-atan(a*x))/(a**2*c*x**2 + c)**2
    F = -(-2*a*x + 1)*exp(-atan(a*x))/(5*a*c**2*(a**2*x**2 + 1)) - 2*exp(-atan(a*x))/(5*a*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_279():
    f = exp(-atan(a*x))/(a**2*c*x**2 + c)**3
    F = -(-24*a*x + 12)*exp(-atan(a*x))/(85*a*c**3*(a**2*x**2 + 1)) - (-4*a*x + 1)*exp(-atan(a*x))/(17*a*c**3*(a**2*x**2 + 1)**2) - 24*exp(-atan(a*x))/(85*a*c**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_280():
    f = exp(-atan(a*x))/(a**2*c*x**2 + c)**4
    F = -(-144*a*x + 72)*exp(-atan(a*x))/(629*a*c**4*(a**2*x**2 + 1)) - (-120*a*x + 30)*exp(-atan(a*x))/(629*a*c**4*(a**2*x**2 + 1)**2) - (-6*a*x + 1)*exp(-atan(a*x))/(37*a*c**4*(a**2*x**2 + 1)**3) - 144*exp(-atan(a*x))/(629*a*c**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_281():
    f = (a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-atan(a*x))
    F = -2**(sympy.S(3)/2 + I/2)*c*(sympy.S(1)/13 - 5*I/13)*(-I*a*x + 1)**(sympy.S(5)/2 - I/2)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-3)/2 - I/2, sympy.S(5)/2 - I/2), (sympy.S(7)/2 - I/2,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_282():
    f = sqrt(a**2*c*x**2 + c)*exp(-atan(a*x))
    F = -2**(sympy.S.Half + I/2)*(sympy.S(1)/5 - 3*I/5)*(-I*a*x + 1)**(sympy.S(3)/2 - I/2)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-1)/2 - I/2, sympy.S(3)/2 - I/2), (sympy.S(5)/2 - I/2,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_283():
    f = exp(-atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2**(sympy.S(-1)/2 + I/2)*(1 - I)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(sympy.S.Half - I/2)*hyper((sympy.S.Half - I/2, sympy.S.Half - I/2), (sympy.S(3)/2 - I/2,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_284():
    f = exp(-atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -(-a*x + 1)*exp(-atan(a*x))/(2*a*c*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_285():
    f = exp(-atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -(-3*a*x + 1)*exp(-atan(a*x))/(10*a*c*(a**2*c*x**2 + c)**(sympy.S(3)/2)) - (-3*a*x + 3)*exp(-atan(a*x))/(10*a*c**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_286():
    f = exp(-atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -(-5*a*x + 1)*exp(-atan(a*x))/(26*a*c*(a**2*c*x**2 + c)**(sympy.S(5)/2)) - (-3*a*x + 1)*exp(-atan(a*x))/(13*a*c**2*(a**2*c*x**2 + c)**(sympy.S(3)/2)) - (-3*a*x + 3)*exp(-atan(a*x))/(13*a*c**3*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_287():
    f = (a**2*c*x**2 + c)**p*exp(-2*atan(a*x))
    F = 2**(p + I)*I*(-I*a*x + 1)**(p + 1 - I)*(a**2*c*x**2 + c)**p*hyper((p + 1 - I, -p - I), (p + 2 - I,), -I*a*x/2 + sympy.S.Half)/(a*(a**2*x**2 + 1)**p*(p + 1 - I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_288():
    f = (a**2*c*x**2 + c)**2*exp(-2*atan(a*x))
    F = -2**(1 + I)*c**2*(sympy.S(1)/5 - 3*I/5)*(-I*a*x + 1)**(3 - I)*hyper((-2 - I, 3 - I), (4 - I,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_289():
    f = (a**2*c*x**2 + c)*exp(-2*atan(a*x))
    F = -2**(1 + I)*c*(sympy.S(1)/5 - 2*I/5)*(-I*a*x + 1)**(2 - I)*hyper((-1 - I, 2 - I), (3 - I,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_290():
    f = exp(-2*atan(a*x))
    F = -2**(-1 + I)*(1 - I)*(-I*a*x + 1)**(1 - I)*hyper((-I, 1 - I), (2 - I,), -I*a*x/2 + sympy.S.Half)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_291():
    f = exp(-2*atan(a*x))/(a**2*c*x**2 + c)
    F = -exp(-2*atan(a*x))/(2*a*c)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_292():
    f = exp(-2*atan(a*x))/(a**2*c*x**2 + c)**2
    F = -(-a*x + 1)*exp(-2*atan(a*x))/(4*a*c**2*(a**2*x**2 + 1)) - exp(-2*atan(a*x))/(8*a*c**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_293():
    f = exp(-2*atan(a*x))/(a**2*c*x**2 + c)**3
    F = -(-3*a*x + 3)*exp(-2*atan(a*x))/(20*a*c**3*(a**2*x**2 + 1)) - (-2*a*x + 1)*exp(-2*atan(a*x))/(10*a*c**3*(a**2*x**2 + 1)**2) - 3*exp(-2*atan(a*x))/(40*a*c**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_294():
    f = exp(-2*atan(a*x))/(a**2*c*x**2 + c)**4
    F = -(-9*a*x + 9)*exp(-2*atan(a*x))/(80*a*c**4*(a**2*x**2 + 1)) - (-6*a*x + 3)*exp(-2*atan(a*x))/(40*a*c**4*(a**2*x**2 + 1)**2) - (-3*a*x + 1)*exp(-2*atan(a*x))/(20*a*c**4*(a**2*x**2 + 1)**3) - 9*exp(-2*atan(a*x))/(160*a*c**4)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_295():
    f = (a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(-2*atan(a*x))
    F = -2**(sympy.S(5)/2 + I)*c*(sympy.S(2)/29 - 5*I/29)*(-I*a*x + 1)**(sympy.S(5)/2 - I)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-3)/2 - I, sympy.S(5)/2 - I), (sympy.S(7)/2 - I,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_296():
    f = sqrt(a**2*c*x**2 + c)*exp(-2*atan(a*x))
    F = -2**(sympy.S(3)/2 + I)*(sympy.S(2)/13 - 3*I/13)*(-I*a*x + 1)**(sympy.S(3)/2 - I)*sqrt(a**2*c*x**2 + c)*hyper((sympy.S(-1)/2 - I, sympy.S(3)/2 - I), (sympy.S(5)/2 - I,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_297():
    f = exp(-2*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2**(sympy.S.Half + I)*(sympy.S(2)/5 - I/5)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(sympy.S.Half - I)*hyper((sympy.S.Half - I, sympy.S.Half - I), (sympy.S(3)/2 - I,), -I*a*x/2 + sympy.S.Half)/(a*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_298():
    f = exp(-2*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -(-a*x + 2)*exp(-2*atan(a*x))/(5*a*c*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_299():
    f = exp(-2*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = -(-3*a*x + 2)*exp(-2*atan(a*x))/(13*a*c*(a**2*c*x**2 + c)**(sympy.S(3)/2)) - (-6*a*x + 12)*exp(-2*atan(a*x))/(65*a*c**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_300():
    f = exp(-2*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(7)/2)
    F = -(-5*a*x + 2)*exp(-2*atan(a*x))/(29*a*c*(a**2*c*x**2 + c)**(sympy.S(5)/2)) - (-60*a*x + 40)*exp(-2*atan(a*x))/(377*a*c**2*(a**2*c*x**2 + c)**(sympy.S(3)/2)) - (-24*a*x + 48)*exp(-2*atan(a*x))/(377*a*c**3*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_301():
    f = exp(5*I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = I*log(a*x + I)/a + 4*I/(a*(-I*a*x + 1)) - 2*I/(a*(-I*a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_302():
    f = exp(4*I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = asinh(a*x)/a + 2*I*sqrt(I*a*x + 1)/(a*sqrt(-I*a*x + 1)) - 2*I*(I*a*x + 1)**(sympy.S(3)/2)/(3*a*(-I*a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_303():
    f = exp(3*I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = -I*log(a*x + I)/a + 2/(a*(a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_304():
    f = exp(2*I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = -asinh(a*x)/a - 2*I*sqrt(I*a*x + 1)/(a*sqrt(-I*a*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_305():
    f = exp(I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = I*log(a*x + I)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_306():
    f = exp(-I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = -I*log(-a*x + I)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_307():
    f = exp(-2*I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = 2*I*sqrt(-I*a*x + 1)/(a*sqrt(I*a*x + 1)) - asinh(a*x)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_308():
    f = exp(-3*I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = I*log(-a*x + I)/a - 2/(a*(-a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_309():
    f = exp(-4*I*atan(a*x))/sqrt(a**2*x**2 + 1)
    F = 2*I*(-I*a*x + 1)**(sympy.S(3)/2)/(3*a*(I*a*x + 1)**(sympy.S(3)/2)) - 2*I*sqrt(-I*a*x + 1)/(a*sqrt(I*a*x + 1)) + asinh(a*x)/a
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_310():
    f = exp(5*I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = I*sqrt(a**2*x**2 + 1)*log(a*x + I)/(a*sqrt(a**2*c*x**2 + c)) + 4*I*sqrt(a**2*x**2 + 1)/(a*(-I*a*x + 1)*sqrt(a**2*c*x**2 + c)) - 2*I*sqrt(a**2*x**2 + 1)/(a*(-I*a*x + 1)**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_311():
    f = exp(4*I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2*I*c*(I*a*x + 1)**3/(3*a*(a**2*c*x**2 + c)**(sympy.S(3)/2)) + 2*I*(I*a*x + 1)/(a*sqrt(a**2*c*x**2 + c)) + atanh(a*sqrt(c)*x/sqrt(a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_312():
    f = exp(3*I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -I*sqrt(a**2*x**2 + 1)*log(a*x + I)/(a*sqrt(a**2*c*x**2 + c)) + 2*sqrt(a**2*x**2 + 1)/(a*(a*x + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_313():
    f = exp(2*I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2*I*(I*a*x + 1)/(a*sqrt(a**2*c*x**2 + c)) - atanh(a*sqrt(c)*x/sqrt(a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_314():
    f = exp(I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = I*sqrt(a**2*x**2 + 1)*log(a*x + I)/(a*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_315():
    f = exp(-I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -I*sqrt(a**2*x**2 + 1)*log(-a*x + I)/(a*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_316():
    f = exp(-2*I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = 2*I*(-I*a*x + 1)/(a*sqrt(a**2*c*x**2 + c)) - atanh(a*sqrt(c)*x/sqrt(a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_317():
    f = exp(-3*I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = I*sqrt(a**2*x**2 + 1)*log(-a*x + I)/(a*sqrt(a**2*c*x**2 + c)) - 2*sqrt(a**2*x**2 + 1)/(a*(-a*x + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_318():
    f = exp(-4*I*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = 2*I*c*(-I*a*x + 1)**3/(3*a*(a**2*c*x**2 + c)**(sympy.S(3)/2)) - 2*I*(-I*a*x + 1)/(a*sqrt(a**2*c*x**2 + c)) + atanh(a*sqrt(c)*x/sqrt(a**2*c*x**2 + c))/(a*sqrt(c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_319():
    f = exp(5*I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = -I/(2*a*(a*x + I)**2) - 2/(3*a*(a*x + I)**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_320():
    f = exp(4*I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = -I*(I*a*x + 1)**(sympy.S(3)/2)/(15*a*(-I*a*x + 1)**(sympy.S(3)/2)) - I*(I*a*x + 1)**(sympy.S(3)/2)/(5*a*(-I*a*x + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_321():
    f = exp(3*I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = -I/(2*a*(-I*a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_322():
    f = exp(2*I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = -I*sqrt(I*a*x + 1)/(3*a*sqrt(-I*a*x + 1)) - I*sqrt(I*a*x + 1)/(3*a*(-I*a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_323():
    f = exp(I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = atan(a*x)/(2*a) + 1/(2*a*(a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_324():
    f = exp(-I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = atan(a*x)/(2*a) - 1/(2*a*(-a*x + I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_325():
    f = exp(-2*I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = I*sqrt(-I*a*x + 1)/(3*a*sqrt(I*a*x + 1)) + I*sqrt(-I*a*x + 1)/(3*a*(I*a*x + 1)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_326():
    f = exp(-3*I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = I/(2*a*(I*a*x + 1)**2)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_327():
    f = exp(-4*I*atan(a*x))/(a**2*x**2 + 1)**(sympy.S(3)/2)
    F = I*(-I*a*x + 1)**(sympy.S(3)/2)/(15*a*(I*a*x + 1)**(sympy.S(3)/2)) + I*(-I*a*x + 1)**(sympy.S(3)/2)/(5*a*(I*a*x + 1)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_328():
    f = exp(5*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -I*sqrt(a**2*x**2 + 1)/(2*a*c*(a*x + I)**2*sqrt(a**2*c*x**2 + c)) - 2*sqrt(a**2*x**2 + 1)/(3*a*c*(a*x + I)**3*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_329():
    f = exp(4*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = I*c*(I*a*x + 1)**5/(15*a*(a**2*c*x**2 + c)**(sympy.S(5)/2)) - I*c*(I*a*x + 1)**4/(3*a*(a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_330():
    f = exp(3*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -I*sqrt(a**2*x**2 + 1)/(2*a*c*(-I*a*x + 1)**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_331():
    f = exp(2*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x/(3*c*sqrt(a**2*c*x**2 + c)) - 2*I*(I*a*x + 1)/(3*a*(a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_332():
    f = exp(I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = sqrt(a**2*x**2 + 1)*atan(a*x)/(2*a*c*sqrt(a**2*c*x**2 + c)) + sqrt(a**2*x**2 + 1)/(2*a*c*(a*x + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_333():
    f = exp(-I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = sqrt(a**2*x**2 + 1)*atan(a*x)/(2*a*c*sqrt(a**2*c*x**2 + c)) - sqrt(a**2*x**2 + 1)/(2*a*c*(-a*x + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_334():
    f = exp(-2*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x/(3*c*sqrt(a**2*c*x**2 + c)) + 2*I*(-I*a*x + 1)/(3*a*(a**2*c*x**2 + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_335():
    f = exp(-3*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = I*sqrt(a**2*x**2 + 1)/(2*a*c*(I*a*x + 1)**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_336():
    f = exp(-4*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -I*c*(-I*a*x + 1)**5/(15*a*(a**2*c*x**2 + c)**(sympy.S(5)/2)) + I*c*(-I*a*x + 1)**4/(3*a*(a**2*c*x**2 + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_337():
    f = (a**2*c*x**2 + c)**2*exp(n*atan(a*x))
    F = -2**(-I*n/2 + 3)*c**2*(-I*a*x + 1)**(I*n/2 + 3)*hyper((I*n/2 - 2, I*n/2 + 3), (I*n/2 + 4,), -I*a*x/2 + sympy.S.Half)/(a*(-n + 6*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_338():
    f = (a**2*c*x**2 + c)*exp(n*atan(a*x))
    F = -2**(-I*n/2 + 2)*c*(-I*a*x + 1)**(I*n/2 + 2)*hyper((I*n/2 - 1, I*n/2 + 2), (I*n/2 + 3,), -I*a*x/2 + sympy.S.Half)/(a*(-n + 4*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_339():
    f = exp(n*atan(a*x))
    F = -2**(-I*n/2 + 1)*(-I*a*x + 1)**(I*n/2 + 1)*hyper((I*n/2, I*n/2 + 1), (I*n/2 + 2,), -I*a*x/2 + sympy.S.Half)/(a*(-n + 2*I))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_340():
    f = x**2*exp(n*atan(a*x))/(a**2*c*x**2 + c)
    F = 2**(-I*n/2 + 1)*I*(-I*a*x + 1)**(I*n/2)*hyper((I*n/2, I*n/2), (I*n/2 + 1,), -I*a*x/2 + sympy.S.Half)/(a**3*c) + x*(-I*a*x + 1)**(I*n/2)/(a**2*c*(I*a*x + 1)**(I*n/2)) - (I*n + 1)*(-I*a*x + 1)**(I*n/2)/(a**3*c*n*(I*a*x + 1)**(I*n/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_341():
    f = x*exp(n*atan(a*x))/(a**2*c*x**2 + c)
    F = -2**(-I*n/2 + 1)*I*(-I*a*x + 1)**(I*n/2)*hyper((I*n/2, I*n/2), (I*n/2 + 1,), -I*a*x/2 + sympy.S.Half)/(a**2*c*n) + I*(-I*a*x + 1)**(I*n/2)/(a**2*c*n*(I*a*x + 1)**(I*n/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_342():
    f = exp(n*atan(a*x))/(a**2*c*x**2 + c)
    F = exp(n*atan(a*x))/(a*c*n)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_343():
    f = exp(n*atan(a*x))/(a**2*c*x**2 + c)**4
    F = (6*a*x + n)*exp(n*atan(a*x))/(a*c**4*(n**2 + 36)*(a**2*x**2 + 1)**3) + 30*(4*a*x + n)*exp(n*atan(a*x))/(a*c**4*(n**2 + 16)*(n**2 + 36)*(a**2*x**2 + 1)**2) + 360*(2*a*x + n)*exp(n*atan(a*x))/(a*c**4*(n**2 + 4)*(n**2 + 16)*(n**2 + 36)*(a**2*x**2 + 1)) + 720*exp(n*atan(a*x))/(a*c**4*n*(n**2 + 4)*(n**2 + 16)*(n**2 + 36))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_344():
    f = (a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(n*atan(a*x))
    F = -2**(-I*n/2 + sympy.S(5)/2)*c*(-I*a*x + 1)**(I*n/2 + sympy.S(5)/2)*sqrt(a**2*c*x**2 + c)*hyper((I*n/2 + sympy.S(-3)/2, I*n/2 + sympy.S(5)/2), (I*n/2 + sympy.S(7)/2,), -I*a*x/2 + sympy.S.Half)/(a*(-n + 5*I)*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_345():
    f = sqrt(a**2*c*x**2 + c)*exp(n*atan(a*x))
    F = -2**(-I*n/2 + sympy.S(3)/2)*(-I*a*x + 1)**(I*n/2 + sympy.S(3)/2)*sqrt(a**2*c*x**2 + c)*hyper((I*n/2 + sympy.S(-1)/2, I*n/2 + sympy.S(3)/2), (I*n/2 + sympy.S(5)/2,), -I*a*x/2 + sympy.S.Half)/(a*(-n + 3*I)*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_346():
    f = exp(n*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2**(-I*n/2 + sympy.S.Half)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*hyper((I*n/2 + sympy.S.Half, I*n/2 + sympy.S.Half), (I*n/2 + sympy.S(3)/2,), -I*a*x/2 + sympy.S.Half)/(a*(-n + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_347():
    f = x**2*(a**2*c*x**2 + c)**(sympy.S(3)/2)*exp(n*atan(a*x))
    F = 2**(-I*n/2 + sympy.S(3)/2)*c*(5 - n**2)*(-I*a*x + 1)**(I*n/2 + sympy.S(5)/2)*sqrt(a**2*c*x**2 + c)*hyper((I*n/2 + sympy.S(-3)/2, I*n/2 + sympy.S(5)/2), (I*n/2 + sympy.S(7)/2,), -I*a*x/2 + sympy.S.Half)/(15*a**3*(-n + 5*I)*sqrt(a**2*x**2 + 1)) + c*x*(-I*a*x + 1)**(I*n/2 + sympy.S(5)/2)*(I*a*x + 1)**(-I*n/2 + sympy.S(5)/2)*sqrt(a**2*c*x**2 + c)/(6*a**2*sqrt(a**2*x**2 + 1)) - c*n*(-I*a*x + 1)**(I*n/2 + sympy.S(5)/2)*(I*a*x + 1)**(-I*n/2 + sympy.S(5)/2)*sqrt(a**2*c*x**2 + c)/(30*a**3*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_348():
    f = x**2*sqrt(a**2*c*x**2 + c)*exp(n*atan(a*x))
    F = 2**(-I*n/2 + sympy.S(-1)/2)*(3 - n**2)*(-I*a*x + 1)**(I*n/2 + sympy.S(3)/2)*sqrt(a**2*c*x**2 + c)*hyper((I*n/2 + sympy.S(-1)/2, I*n/2 + sympy.S(3)/2), (I*n/2 + sympy.S(5)/2,), -I*a*x/2 + sympy.S.Half)/(3*a**3*(-n + 3*I)*sqrt(a**2*x**2 + 1)) + x*(-I*a*x + 1)**(I*n/2 + sympy.S(3)/2)*(I*a*x + 1)**(-I*n/2 + sympy.S(3)/2)*sqrt(a**2*c*x**2 + c)/(4*a**2*sqrt(a**2*x**2 + 1)) - n*(-I*a*x + 1)**(I*n/2 + sympy.S(3)/2)*(I*a*x + 1)**(-I*n/2 + sympy.S(3)/2)*sqrt(a**2*c*x**2 + c)/(12*a**3*sqrt(a**2*x**2 + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_349():
    f = x**3*exp(n*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = 2**(-I*n/2 + sympy.S(-1)/2)*n*(5 - n**2)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S(3)/2)*hyper((I*n/2 + sympy.S.Half, I*n/2 + sympy.S(3)/2), (I*n/2 + sympy.S(5)/2,), -I*a*x/2 + sympy.S.Half)/(3*a**4*(4*n - I*(3 - n**2))*sqrt(a**2*c*x**2 + c)) + x**2*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)/(3*a**2*sqrt(a**2*c*x**2 + c)) - sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)*(a*n*x*(I*n + 1) - n**2 - I*n + 4)/(6*a**4*(I*n + 1)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_350():
    f = x**2*exp(n*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2**(-I*n/2 + sympy.S.Half)*I*(1 - n**2)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*hyper((I*n/2 + sympy.S(-1)/2, I*n/2 + sympy.S.Half), (I*n/2 + sympy.S(3)/2,), -I*a*x/2 + sympy.S.Half)/(a**3*(n**2 + 1)*sqrt(a**2*c*x**2 + c)) + x*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)/(2*a**2*sqrt(a**2*c*x**2 + c)) - (I*n + 1)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)/(2*a**3*(n + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_351():
    f = x*exp(n*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2**(-I*n/2 + sympy.S(3)/2)*I*n*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*hyper((I*n/2 + sympy.S(-1)/2, I*n/2 + sympy.S.Half), (I*n/2 + sympy.S(3)/2,), -I*a*x/2 + sympy.S.Half)/(a**2*(n**2 + 1)*sqrt(a**2*c*x**2 + c)) + sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)/(a**2*(-I*n + 1)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_352():
    f = exp(n*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = -2**(-I*n/2 + sympy.S.Half)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*hyper((I*n/2 + sympy.S.Half, I*n/2 + sympy.S.Half), (I*n/2 + sympy.S(3)/2,), -I*a*x/2 + sympy.S.Half)/(a*(-n + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_353():
    f = exp(n*atan(a*x))/(x*sqrt(a**2*c*x**2 + c))
    F = -2*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S(-1)/2)*hyper((1, I*n/2 + sympy.S.Half), (I*n/2 + sympy.S(3)/2,), (-I*a*x + 1)/(I*a*x + 1))/((I*n + 1)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_354():
    f = exp(n*atan(a*x))/(x**2*sqrt(a**2*c*x**2 + c))
    F = -2*a*n*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S(-1)/2)*hyper((1, I*n/2 + sympy.S.Half), (I*n/2 + sympy.S(3)/2,), (-I*a*x + 1)/(I*a*x + 1))/((I*n + 1)*sqrt(a**2*c*x**2 + c)) - sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)/(x*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_355():
    f = exp(n*atan(a*x))/(x**3*sqrt(a**2*c*x**2 + c))
    F = a**2*(1 - n**2)*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S(-1)/2)*hyper((1, I*n/2 + sympy.S.Half), (I*n/2 + sympy.S(3)/2,), (-I*a*x + 1)/(I*a*x + 1))/((I*n + 1)*sqrt(a**2*c*x**2 + c)) - a*n*sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)/(2*x*sqrt(a**2*c*x**2 + c)) - sqrt(a**2*x**2 + 1)*(-I*a*x + 1)**(I*n/2 + sympy.S.Half)*(I*a*x + 1)**(-I*n/2 + sympy.S.Half)/(2*x**2*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_356():
    f = (a**2*c*x**2 + c)**(sympy.S(1)/3)*exp(n*atan(a*x))
    F = -3*2**(-I*n/2 + sympy.S(4)/3)*(-I*a*x + 1)**(I*n/2 + sympy.S(4)/3)*(a**2*c*x**2 + c)**(sympy.S(1)/3)*hyper((I*n/2 + sympy.S(-1)/3, I*n/2 + sympy.S(4)/3), (I*n/2 + sympy.S(7)/3,), -I*a*x/2 + sympy.S.Half)/(a*(-3*n + 8*I)*(a**2*x**2 + 1)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_357():
    f = exp(n*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(1)/3)
    F = -3*2**(-I*n/2 + sympy.S(2)/3)*(a**2*x**2 + 1)**(sympy.S(1)/3)*(-I*a*x + 1)**(I*n/2 + sympy.S(2)/3)*hyper((I*n/2 + sympy.S(1)/3, I*n/2 + sympy.S(2)/3), (I*n/2 + sympy.S(5)/3,), -I*a*x/2 + sympy.S.Half)/(a*(-3*n + 4*I)*(a**2*c*x**2 + c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_358():
    f = exp(n*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(2)/3)
    F = -3*2**(-I*n/2 + sympy.S(1)/3)*(a**2*x**2 + 1)**(sympy.S(2)/3)*(-I*a*x + 1)**(I*n/2 + sympy.S(1)/3)*hyper((I*n/2 + sympy.S(1)/3, I*n/2 + sympy.S(2)/3), (I*n/2 + sympy.S(4)/3,), -I*a*x/2 + sympy.S.Half)/(a*(-3*n + 2*I)*(a**2*c*x**2 + c)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_359():
    f = exp(n*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(4)/3)
    F = 3*2**(-I*n/2 + sympy.S(-1)/3)*(a**2*x**2 + 1)**(sympy.S(1)/3)*(-I*a*x + 1)**(I*n/2 + sympy.S(-1)/3)*hyper((I*n/2 + sympy.S(-1)/3, I*n/2 + sympy.S(4)/3), (I*n/2 + sympy.S(2)/3,), -I*a*x/2 + sympy.S.Half)/(a*c*(3*n + 2*I)*(a**2*c*x**2 + c)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_360():
    f = x**m*(a**2*c*x**2 + c)*exp(n*atan(a*x))
    F = c*x**(m + 1)*appellf1(m + 1, -I*n/2 - 1, I*n/2 - 1, m + 2, I*a*x, -I*a*x)/(m + 1)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_361():
    f = x**m*exp(n*atan(a*x))/(a**2*c*x**2 + c)
    F = x**(m + 1)*appellf1(m + 1, -I*n/2 + 1, I*n/2 + 1, m + 2, I*a*x, -I*a*x)/(c*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_362():
    f = x**m*exp(n*atan(a*x))/(a**2*c*x**2 + c)**2
    F = x**(m + 1)*appellf1(m + 1, -I*n/2 + 2, I*n/2 + 2, m + 2, I*a*x, -I*a*x)/(c**2*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_363():
    f = x**m*exp(n*atan(a*x))/(a**2*c*x**2 + c)**3
    F = x**(m + 1)*appellf1(m + 1, -I*n/2 + 3, I*n/2 + 3, m + 2, I*a*x, -I*a*x)/(c**3*(m + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_364():
    f = x**m*exp(n*atan(a*x))/sqrt(a**2*c*x**2 + c)
    F = x**(m + 1)*sqrt(a**2*x**2 + 1)*appellf1(m + 1, -I*n/2 + sympy.S.Half, I*n/2 + sympy.S.Half, m + 2, I*a*x, -I*a*x)/((m + 1)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_365():
    f = x**m*exp(n*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = x**(m + 1)*sqrt(a**2*x**2 + 1)*appellf1(m + 1, -I*n/2 + sympy.S(3)/2, I*n/2 + sympy.S(3)/2, m + 2, I*a*x, -I*a*x)/(c*(m + 1)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_366():
    f = x**m*exp(n*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(5)/2)
    F = x**(m + 1)*sqrt(a**2*x**2 + 1)*appellf1(m + 1, -I*n/2 + sympy.S(5)/2, I*n/2 + sympy.S(5)/2, m + 2, I*a*x, -I*a*x)/(c**2*(m + 1)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_367():
    f = (a**2*c*x**2 + c)**p*exp(n*atan(a*x))
    F = 2**(-I*n/2 + p + 1)*(-I*a*x + 1)**(I*n/2 + p + 1)*(a**2*c*x**2 + c)**p*hyper((I*n/2 + p + 1, I*n/2 - p), (I*n/2 + p + 2,), -I*a*x/2 + sympy.S.Half)/(a*(n - 2*I*(p + 1))*(a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_368():
    f = (a**2*c*x**2 + c)**p*exp(-2*I*p*atan(a*x))
    F = I*(-I*a*x + 1)**(2*p + 1)*(a**2*c*x**2 + c)**p/(a*(2*p + 1)*(a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_369():
    f = (a**2*c*x**2 + c)**p*exp(2*I*p*atan(a*x))
    F = -I*(I*a*x + 1)**(2*p + 1)*(a**2*c*x**2 + c)**p/(a*(2*p + 1)*(a**2*x**2 + 1)**p)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_370():
    f = x**2*(a**2*c*x**2 + c)**(-n**2/2 - 1)*exp(I*n*atan(a*x))
    F = I*(-I*a*n*x + 1)*exp(I*n*atan(a*x))/(a**3*c*n*(1 - n**2)*(a**2*c*x**2 + c)**(n**2/2))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_371():
    f = x**2*exp(6*I*atan(a*x))/(a**2*c*x**2 + c)**19
    F = -(6*a*x + I)/(210*a**3*c**19*(-I*a*x + 1)**21*(I*a*x + 1)**15)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_372():
    f = x**2*exp(4*I*atan(a*x))/(a**2*c*x**2 + c)**9
    F = -(4*a*x + I)/(60*a**3*c**9*(-I*a*x + 1)**10*(I*a*x + 1)**6)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_373():
    f = x**2*exp(2*I*atan(a*x))/(a**2*c*x**2 + c)**3
    F = -(2*a*x + I)/(6*a**3*c**3*(-I*a*x + 1)**3*(I*a*x + 1))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_374():
    f = x**2*exp(-2*I*atan(a*x))/(a**2*c*x**2 + c)**3
    F = (-2*a*x + I)/(6*a**3*c**3*(-I*a*x + 1)*(I*a*x + 1)**3)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_375():
    f = x**2*exp(-4*I*atan(a*x))/(a**2*c*x**2 + c)**9
    F = (-4*a*x + I)/(60*a**3*c**9*(-I*a*x + 1)**6*(I*a*x + 1)**10)
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_376():
    f = x**2*exp(5*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(27)/2)
    F = -(5*a*x + I)*sqrt(a**2*x**2 + 1)/(120*a**3*c**13*(-I*a*x + 1)**15*(I*a*x + 1)**10*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_377():
    f = x**2*exp(3*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(11)/2)
    F = -(3*a*x + I)*sqrt(a**2*x**2 + 1)/(24*a**3*c**5*(-I*a*x + 1)**6*(I*a*x + 1)**3*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_378():
    f = x**2*exp(I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = I*sqrt(a**2*x**2 + 1)*log(-a*x + I)/(4*a**3*c*sqrt(a**2*c*x**2 + c)) + 3*I*sqrt(a**2*x**2 + 1)*log(a*x + I)/(4*a**3*c*sqrt(a**2*c*x**2 + c)) - sqrt(a**2*x**2 + 1)/(2*a**3*c*(a*x + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_379():
    f = x**2*exp(-I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(3)/2)
    F = -3*I*sqrt(a**2*x**2 + 1)*log(-a*x + I)/(4*a**3*c*sqrt(a**2*c*x**2 + c)) - I*sqrt(a**2*x**2 + 1)*log(a*x + I)/(4*a**3*c*sqrt(a**2*c*x**2 + c)) + sqrt(a**2*x**2 + 1)/(2*a**3*c*(-a*x + I)*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_380():
    f = x**2*exp(-3*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(11)/2)
    F = (-3*a*x + I)*sqrt(a**2*x**2 + 1)/(24*a**3*c**5*(-I*a*x + 1)**3*(I*a*x + 1)**6*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F


def test_integrate_5_Inverse_trig_functions_5_3_Inverse_tangent_5_3_6_Exponentials_of_inverse_tangent_381():
    f = x**2*exp(-5*I*atan(a*x))/(a**2*c*x**2 + c)**(sympy.S(27)/2)
    F = (-5*a*x + I)*sqrt(a**2*x**2 + 1)/(120*a**3*c**13*(-I*a*x + 1)**10*(I*a*x + 1)**15*sqrt(a**2*c*x**2 + c))
    assert integrate(f, x) == F

