"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.5 Hyperbolic secant/6.5.3 Hyperbolic secant functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, n, p = symbols('a b c d n p')

def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_1():
    f = sech(a + b*x)
    F = atan(sinh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_2():
    f = sech(a + b*x)**2
    F = tanh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_3():
    f = sech(a + b*x)**3
    F = tanh(a + b*x)*sech(a + b*x)/(2*b) + atan(sinh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_4():
    f = sech(a + b*x)**4
    F = -tanh(a + b*x)**3/(3*b) + tanh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_5():
    f = sech(a + b*x)**5
    F = tanh(a + b*x)*sech(a + b*x)**3/(4*b) + 3*tanh(a + b*x)*sech(a + b*x)/(8*b) + 3*atan(sinh(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_6():
    f = sech(a + b*x)**6
    F = tanh(a + b*x)**5/(5*b) - 2*tanh(a + b*x)**3/(3*b) + tanh(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_7():
    f = sech(7*x)**4
    F = -tanh(7*x)**3/21 + tanh(7*x)/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_8():
    f = sech(pi*x)**6
    F = tanh(pi*x)**5/(5*pi) - 2*tanh(pi*x)**3/(3*pi) + tanh(pi*x)/pi
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_9():
    f = sech(a + b*x)**(sympy.S(5)/2)
    F = 2*sinh(a + b*x)*sech(a + b*x)**(sympy.S(3)/2)/(3*b) - 2*I*sqrt(cosh(a + b*x))*elliptic_f(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_10():
    f = sech(a + b*x)**(sympy.S(3)/2)
    F = 2*sinh(a + b*x)*sqrt(sech(a + b*x))/b + 2*I*sqrt(cosh(a + b*x))*elliptic_e(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_11():
    f = sqrt(sech(a + b*x))
    F = -2*I*sqrt(cosh(a + b*x))*elliptic_f(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_12():
    f = 1/sqrt(sech(a + b*x))
    F = -2*I*sqrt(cosh(a + b*x))*elliptic_e(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_13():
    f = sech(a + b*x)**(sympy.S(-3)/2)
    F = 2*sinh(a + b*x)/(3*b*sqrt(sech(a + b*x))) - 2*I*sqrt(cosh(a + b*x))*elliptic_f(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_14():
    f = sech(a + b*x)**(sympy.S(-5)/2)
    F = 2*sinh(a + b*x)/(5*b*sech(a + b*x)**(sympy.S(3)/2)) - 6*I*sqrt(cosh(a + b*x))*elliptic_e(I*(a + b*x)/2, 2)*sqrt(sech(a + b*x))/(5*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_15():
    f = (b*sech(c + d*x))**(sympy.S(7)/2)
    F = 6*I*b**4*elliptic_e(I*(c + d*x)/2, 2)/(5*d*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x))) + 6*b**3*sqrt(b*sech(c + d*x))*sinh(c + d*x)/(5*d) + 2*b*(b*sech(c + d*x))**(sympy.S(5)/2)*sinh(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_16():
    f = (b*sech(c + d*x))**(sympy.S(5)/2)
    F = -2*I*b**2*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x))*elliptic_f(I*(c + d*x)/2, 2)/(3*d) + 2*b*(b*sech(c + d*x))**(sympy.S(3)/2)*sinh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_17():
    f = (b*sech(c + d*x))**(sympy.S(3)/2)
    F = 2*I*b**2*elliptic_e(I*(c + d*x)/2, 2)/(d*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x))) + 2*b*sqrt(b*sech(c + d*x))*sinh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_18():
    f = sqrt(b*sech(c + d*x))
    F = -2*I*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x))*elliptic_f(I*(c + d*x)/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_19():
    f = 1/sqrt(b*sech(c + d*x))
    F = -2*I*elliptic_e(I*(c + d*x)/2, 2)/(d*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_20():
    f = (b*sech(c + d*x))**(sympy.S(-3)/2)
    F = 2*sinh(c + d*x)/(3*b*d*sqrt(b*sech(c + d*x))) - 2*I*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x))*elliptic_f(I*(c + d*x)/2, 2)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_21():
    f = (b*sech(c + d*x))**(sympy.S(-5)/2)
    F = 2*sinh(c + d*x)/(5*b*d*(b*sech(c + d*x))**(sympy.S(3)/2)) - 6*I*elliptic_e(I*(c + d*x)/2, 2)/(5*b**2*d*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_22():
    f = (b*sech(c + d*x))**(sympy.S(-7)/2)
    F = 2*sinh(c + d*x)/(7*b*d*(b*sech(c + d*x))**(sympy.S(5)/2)) + 10*sinh(c + d*x)/(21*b**3*d*sqrt(b*sech(c + d*x))) - 10*I*sqrt(b*sech(c + d*x))*sqrt(cosh(c + d*x))*elliptic_f(I*(c + d*x)/2, 2)/(21*b**4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_23():
    f = (b*sech(c + d*x))**n
    F = -b*(b*sech(c + d*x))**(n - 1)*sinh(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cosh(c + d*x)**2)/(d*sqrt(-sinh(c + d*x)**2)*(1 - n))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_24():
    f = (sech(a + b*x)**2)**(sympy.S(7)/2)
    F = (sech(a + b*x)**2)**(sympy.S(5)/2)*tanh(a + b*x)/(6*b) + 5*(sech(a + b*x)**2)**(sympy.S(3)/2)*tanh(a + b*x)/(24*b) + 5*sqrt(sech(a + b*x)**2)*tanh(a + b*x)/(16*b) + 5*asin(tanh(a + b*x))/(16*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_25():
    f = (sech(a + b*x)**2)**(sympy.S(5)/2)
    F = (sech(a + b*x)**2)**(sympy.S(3)/2)*tanh(a + b*x)/(4*b) + 3*sqrt(sech(a + b*x)**2)*tanh(a + b*x)/(8*b) + 3*asin(tanh(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_26():
    f = (sech(a + b*x)**2)**(sympy.S(3)/2)
    F = sqrt(sech(a + b*x)**2)*tanh(a + b*x)/(2*b) + asin(tanh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_27():
    f = sqrt(sech(a + b*x)**2)
    F = asin(tanh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_28():
    f = 1/sqrt(sech(a + b*x)**2)
    F = tanh(a + b*x)/(b*sqrt(sech(a + b*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_29():
    f = (sech(a + b*x)**2)**(sympy.S(-3)/2)
    F = 2*tanh(a + b*x)/(3*b*sqrt(sech(a + b*x)**2)) + tanh(a + b*x)/(3*b*(sech(a + b*x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_30():
    f = (sech(a + b*x)**2)**(sympy.S(-5)/2)
    F = 8*tanh(a + b*x)/(15*b*sqrt(sech(a + b*x)**2)) + 4*tanh(a + b*x)/(15*b*(sech(a + b*x)**2)**(sympy.S(3)/2)) + tanh(a + b*x)/(5*b*(sech(a + b*x)**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_31():
    f = (sech(a + b*x)**2)**(sympy.S(-7)/2)
    F = 16*tanh(a + b*x)/(35*b*sqrt(sech(a + b*x)**2)) + 8*tanh(a + b*x)/(35*b*(sech(a + b*x)**2)**(sympy.S(3)/2)) + 6*tanh(a + b*x)/(35*b*(sech(a + b*x)**2)**(sympy.S(5)/2)) + tanh(a + b*x)/(7*b*(sech(a + b*x)**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_32():
    f = (a*sech(x)**2)**(sympy.S(5)/2)
    F = 3*a**(sympy.S(5)/2)*atan(sqrt(a)*tanh(x)/sqrt(a*sech(x)**2))/8 + 3*a**2*sqrt(a*sech(x)**2)*tanh(x)/8 + a*(a*sech(x)**2)**(sympy.S(3)/2)*tanh(x)/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_33():
    f = (a*sech(x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atan(sqrt(a)*tanh(x)/sqrt(a*sech(x)**2))/2 + a*sqrt(a*sech(x)**2)*tanh(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_34():
    f = sqrt(a*sech(x)**2)
    F = sqrt(a)*atan(sqrt(a)*tanh(x)/sqrt(a*sech(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_35():
    f = 1/sqrt(a*sech(x)**2)
    F = tanh(x)/sqrt(a*sech(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_36():
    f = (a*sech(x)**2)**(sympy.S(-3)/2)
    F = tanh(x)/(3*(a*sech(x)**2)**(sympy.S(3)/2)) + 2*tanh(x)/(3*a*sqrt(a*sech(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_37():
    f = (a*sech(x)**2)**(sympy.S(-5)/2)
    F = tanh(x)/(5*(a*sech(x)**2)**(sympy.S(5)/2)) + 4*tanh(x)/(15*a*(a*sech(x)**2)**(sympy.S(3)/2)) + 8*tanh(x)/(15*a**2*sqrt(a*sech(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_38():
    f = (a*sech(x)**2)**(sympy.S(-7)/2)
    F = tanh(x)/(7*(a*sech(x)**2)**(sympy.S(7)/2)) + 6*tanh(x)/(35*a*(a*sech(x)**2)**(sympy.S(5)/2)) + 8*tanh(x)/(35*a**2*(a*sech(x)**2)**(sympy.S(3)/2)) + 16*tanh(x)/(35*a**3*sqrt(a*sech(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_39():
    f = (a*sech(x)**3)**(sympy.S(5)/2)
    F = 154*a**2*sqrt(a*sech(x)**3)*sinh(x)*cosh(x)/195 + 154*I*a**2*sqrt(a*sech(x)**3)*cosh(x)**(sympy.S(3)/2)*elliptic_e(I*x/2, 2)/195 + 2*a**2*sqrt(a*sech(x)**3)*tanh(x)*sech(x)**4/13 + 22*a**2*sqrt(a*sech(x)**3)*tanh(x)*sech(x)**2/117 + 154*a**2*sqrt(a*sech(x)**3)*tanh(x)/585
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_40():
    f = (a*sech(x)**3)**(sympy.S(3)/2)
    F = 10*a*sqrt(a*sech(x)**3)*sinh(x)/21 - 10*I*a*sqrt(a*sech(x)**3)*cosh(x)**(sympy.S(3)/2)*elliptic_f(I*x/2, 2)/21 + 2*a*sqrt(a*sech(x)**3)*tanh(x)*sech(x)/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_41():
    f = sqrt(a*sech(x)**3)
    F = 2*sqrt(a*sech(x)**3)*sinh(x)*cosh(x) + 2*I*sqrt(a*sech(x)**3)*cosh(x)**(sympy.S(3)/2)*elliptic_e(I*x/2, 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_42():
    f = 1/sqrt(a*sech(x)**3)
    F = 2*tanh(x)/(3*sqrt(a*sech(x)**3)) - 2*I*elliptic_f(I*x/2, 2)/(3*sqrt(a*sech(x)**3)*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_43():
    f = (a*sech(x)**3)**(sympy.S(-3)/2)
    F = 2*sinh(x)*cosh(x)**2/(9*a*sqrt(a*sech(x)**3)) + 14*sinh(x)/(45*a*sqrt(a*sech(x)**3)) - 14*I*elliptic_e(I*x/2, 2)/(15*a*sqrt(a*sech(x)**3)*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_44():
    f = (a*sech(x)**3)**(sympy.S(-5)/2)
    F = 2*sinh(x)*cosh(x)**5/(15*a**2*sqrt(a*sech(x)**3)) + 26*sinh(x)*cosh(x)**3/(165*a**2*sqrt(a*sech(x)**3)) + 78*sinh(x)*cosh(x)/(385*a**2*sqrt(a*sech(x)**3)) + 26*tanh(x)/(77*a**2*sqrt(a*sech(x)**3)) - 26*I*elliptic_f(I*x/2, 2)/(77*a**2*sqrt(a*sech(x)**3)*cosh(x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_45():
    f = (a*sech(x)**4)**(sympy.S(7)/2)
    F = a**3*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**11/13 - 6*a**3*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**9/11 + 5*a**3*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**7/3 - 20*a**3*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**5/7 + 3*a**3*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**3 - 2*a**3*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x) + a**3*sqrt(a*sech(x)**4)*sinh(x)*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_46():
    f = (a*sech(x)**4)**(sympy.S(5)/2)
    F = a**2*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**7/9 - 4*a**2*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**5/7 + 6*a**2*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**3/5 - 4*a**2*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)/3 + a**2*sqrt(a*sech(x)**4)*sinh(x)*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_47():
    f = (a*sech(x)**4)**(sympy.S(3)/2)
    F = a*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)**3/5 - 2*a*sqrt(a*sech(x)**4)*sinh(x)**2*tanh(x)/3 + a*sqrt(a*sech(x)**4)*sinh(x)*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_48():
    f = sqrt(a*sech(x)**4)
    F = sqrt(a*sech(x)**4)*sinh(x)*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_49():
    f = 1/sqrt(a*sech(x)**4)
    F = x*sech(x)**2/(2*sqrt(a*sech(x)**4)) + tanh(x)/(2*sqrt(a*sech(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_50():
    f = (a*sech(x)**4)**(sympy.S(-3)/2)
    F = 5*x*sech(x)**2/(16*a*sqrt(a*sech(x)**4)) + sinh(x)*cosh(x)**3/(6*a*sqrt(a*sech(x)**4)) + 5*sinh(x)*cosh(x)/(24*a*sqrt(a*sech(x)**4)) + 5*tanh(x)/(16*a*sqrt(a*sech(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_51():
    f = (a*sech(x)**4)**(sympy.S(-5)/2)
    F = 63*x*sech(x)**2/(256*a**2*sqrt(a*sech(x)**4)) + sinh(x)*cosh(x)**7/(10*a**2*sqrt(a*sech(x)**4)) + 9*sinh(x)*cosh(x)**5/(80*a**2*sqrt(a*sech(x)**4)) + 21*sinh(x)*cosh(x)**3/(160*a**2*sqrt(a*sech(x)**4)) + 21*sinh(x)*cosh(x)/(128*a**2*sqrt(a*sech(x)**4)) + 63*tanh(x)/(256*a**2*sqrt(a*sech(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_52():
    f = sinh(x)**4/(a*sech(x) + a)
    F = -x/(8*a) - sinh(x)**3/(3*a) + sinh(x)*cosh(x)**3/(4*a) - sinh(x)*cosh(x)/(8*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_53():
    f = sinh(x)**3/(a*sech(x) + a)
    F = -sinh(x)**2/(2*a) + cosh(x)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_54():
    f = sinh(x)**2/(a*sech(x) + a)
    F = x/(2*a) + sinh(x)*cosh(x)/(2*a) - sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_55():
    f = sinh(x)/(a*sech(x) + a)
    F = -log(cosh(x) + 1)/a + cosh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_56():
    f = csch(x)/(a*sech(x) + a)
    F = -coth(x)*csch(x)/(2*a) - atanh(cosh(x))/(2*a) + csch(x)**2/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_57():
    f = csch(x)**2/(a*sech(x) + a)
    F = -coth(x)**3/(3*a) + csch(x)**3/(3*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_58():
    f = csch(x)**3/(a*sech(x) + a)
    F = -coth(x)*csch(x)**3/(4*a) - coth(x)*csch(x)/(8*a) + atanh(cosh(x))/(8*a) + csch(x)**4/(4*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_59():
    f = csch(x)**4/(a*sech(x) + a)
    F = -coth(x)**5/(5*a) + coth(x)**3/(3*a) + csch(x)**5/(5*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_60():
    f = sinh(x)**4/(a + b*sech(x))
    F = -(-3*a*cosh(x) + 4*b)*sinh(x)**3/(12*a**2) + (-a*(3*a**2 - 4*b**2)*cosh(x) + 8*b*(a**2 - b**2))*sinh(x)/(8*a**4) - 2*b*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/a**5 + x*(3*a**4 - 12*a**2*b**2 + 8*b**4)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_61():
    f = sinh(x)**3/(a + b*sech(x))
    F = cosh(x)**3/(3*a) - b*cosh(x)**2/(2*a**2) - (a**2 - b**2)*cosh(x)/a**3 + b*(a**2 - b**2)*log(a*cosh(x) + b)/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_62():
    f = sinh(x)**2/(a + b*sech(x))
    F = -(-a*cosh(x) + 2*b)*sinh(x)/(2*a**2) + 2*b*sqrt(a - b)*sqrt(a + b)*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/a**3 - x*(a**2 - 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_63():
    f = sinh(x)/(a + b*sech(x))
    F = cosh(x)/a - b*log(a*cosh(x) + b)/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_64():
    f = csch(x)/(a + b*sech(x))
    F = b*log(a*cosh(x) + b)/(a**2 - b**2) + log(1 - cosh(x))/(2*a + 2*b) - log(cosh(x) + 1)/(2*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_65():
    f = csch(x)**2/(a + b*sech(x))
    F = 2*a*b*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + (-a*cosh(x) + b)*csch(x)/(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_66():
    f = csch(x)**3/(a + b*sech(x))
    F = -a**2*b*log(a*cosh(x) + b)/(a**2 - b**2)**2 - a*log(1 - cosh(x))/(4*(a + b)**2) + a*log(cosh(x) + 1)/(4*(a - b)**2) + (-a*cosh(x) + b)*csch(x)**2/(2*a**2 - 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_67():
    f = csch(x)**4/(a + b*sech(x))
    F = -2*a**3*b*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/((a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + (-a*cosh(x) + b)*csch(x)**3/(3*a**2 - 3*b**2) - (3*a**2*b - a*(2*a**2 + b**2)*cosh(x))*csch(x)/(3*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_68():
    f = cosh(x)**4/(a*sech(x) + a)
    F = -sinh(x)*cosh(x)**3/(a*sech(x) + a) + 15*x/(8*a) - 4*sinh(x)**3/(3*a) + 5*sinh(x)*cosh(x)**3/(4*a) + 15*sinh(x)*cosh(x)/(8*a) - 4*sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_69():
    f = cosh(x)**3/(a*sech(x) + a)
    F = -sinh(x)*cosh(x)**2/(a*sech(x) + a) - 3*x/(2*a) + 4*sinh(x)**3/(3*a) - 3*sinh(x)*cosh(x)/(2*a) + 4*sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_70():
    f = cosh(x)**2/(a*sech(x) + a)
    F = -sinh(x)*cosh(x)/(a*sech(x) + a) + 3*x/(2*a) + 3*sinh(x)*cosh(x)/(2*a) - 2*sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_71():
    f = cosh(x)/(a*sech(x) + a)
    F = -sinh(x)/(a*sech(x) + a) - x/a + 2*sinh(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_72():
    f = sech(x)/(a*sech(x) + a)
    F = tanh(x)/(a*sech(x) + a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_73():
    f = sech(x)**2/(a*sech(x) + a)
    F = -tanh(x)/(a*sech(x) + a) + atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_74():
    f = sech(x)**3/(a*sech(x) + a)
    F = tanh(x)/(a*sech(x) + a) + tanh(x)/a - atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_75():
    f = sech(x)**4/(a*sech(x) + a)
    F = -tanh(x)*sech(x)**2/(a*sech(x) + a) + 3*tanh(x)*sech(x)/(2*a) - 2*tanh(x)/a + 3*atan(sinh(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_76():
    f = 1/(a*sech(c + d*x) + a)
    F = -tanh(c + d*x)/(d*(a*sech(c + d*x) + a)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_77():
    f = 1/(-a*sech(c + d*x) + a)
    F = -tanh(c + d*x)/(d*(-a*sech(c + d*x) + a)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_78():
    f = (a*sech(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(5)/2)*atanh(sqrt(a)*tanh(c + d*x)/sqrt(a*sech(c + d*x) + a))/d + 14*a**3*tanh(c + d*x)/(3*d*sqrt(a*sech(c + d*x) + a)) + 2*a**2*sqrt(a*sech(c + d*x) + a)*tanh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_79():
    f = (a*sech(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(3)/2)*atanh(sqrt(a)*tanh(c + d*x)/sqrt(a*sech(c + d*x) + a))/d + 2*a**2*tanh(c + d*x)/(d*sqrt(a*sech(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_80():
    f = sqrt(a*sech(c + d*x) + a)
    F = 2*sqrt(a)*atanh(sqrt(a)*tanh(c + d*x)/sqrt(a*sech(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_81():
    f = 1/sqrt(a*sech(c + d*x) + a)
    F = 2*atanh(sqrt(a)*tanh(c + d*x)/sqrt(a*sech(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*tanh(c + d*x)/(2*sqrt(a*sech(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_82():
    f = (a*sech(c + d*x) + a)**(sympy.S(-3)/2)
    F = -tanh(c + d*x)/(2*d*(a*sech(c + d*x) + a)**(sympy.S(3)/2)) + 2*atanh(sqrt(a)*tanh(c + d*x)/sqrt(a*sech(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*tanh(c + d*x)/(2*sqrt(a*sech(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_83():
    f = sqrt(-a*sech(c + d*x) + a)
    F = 2*sqrt(a)*atanh(sqrt(a)*tanh(c + d*x)/sqrt(-a*sech(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_84():
    f = 1/sqrt(-a*sech(c + d*x) + a)
    F = 2*atanh(sqrt(a)*tanh(c + d*x)/sqrt(-a*sech(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*tanh(c + d*x)/(2*sqrt(-a*sech(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_85():
    f = sqrt(3*sech(x) + 3)
    F = 2*sqrt(3)*atanh(tanh(x)/sqrt(sech(x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_86():
    f = sqrt(3 - 3*sech(x))
    F = 2*sqrt(3)*atanh(tanh(x)/sqrt(1 - sech(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_87():
    f = (a + b*sech(c + d*x))**4
    F = a**4*x + 4*a*b**3*tanh(c + d*x)*sech(c + d*x)/(3*d) + 2*a*b*(2*a**2 + b**2)*atan(sinh(c + d*x))/d + b**2*(a + b*sech(c + d*x))**2*tanh(c + d*x)/(3*d) + b**2*(17*a**2 + 2*b**2)*tanh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_88():
    f = (a + b*sech(c + d*x))**3
    F = a**3*x + 5*a*b**2*tanh(c + d*x)/(2*d) + b**2*(a + b*sech(c + d*x))*tanh(c + d*x)/(2*d) + b*(6*a**2 + b**2)*atan(sinh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_89():
    f = (a + b*sech(c + d*x))**2
    F = a**2*x + 2*a*b*atan(sinh(c + d*x))/d + b**2*tanh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_90():
    f = a + b*sech(c + d*x)
    F = a*x + b*atan(sinh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_91():
    f = 1/(a + b*sech(c + d*x))
    F = -2*b*atan(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(a*d*sqrt(a - b)*sqrt(a + b)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_92():
    f = (a + b*sech(c + d*x))**(-2)
    F = b**2*tanh(c + d*x)/(a*d*(a + b*sech(c + d*x))*(a**2 - b**2)) - 2*b*(2*a**2 - b**2)*atan(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_93():
    f = (a + b*sech(c + d*x))**(-3)
    F = b**2*tanh(c + d*x)/(2*a*d*(a + b*sech(c + d*x))**2*(a**2 - b**2)) + b**2*(5*a**2 - 2*b**2)*tanh(c + d*x)/(2*a**2*d*(a + b*sech(c + d*x))*(a**2 - b**2)**2) - b*(6*a**4 - 5*a**2*b**2 + 2*b**4)*atan(sqrt(a - b)*tanh(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_94():
    f = 1/sqrt(a + b*sech(c + d*x))
    F = (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_95():
    f = cosh(x)**4/(a + b*sech(x))
    F = sinh(x)*cosh(x)**3/(4*a) - b*sinh(x)*cosh(x)**2/(3*a**2) + (3*a**2 + 4*b**2)*sinh(x)*cosh(x)/(8*a**3) - b*(2*a**2 + 3*b**2)*sinh(x)/(3*a**4) - 2*b**5*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**5*sqrt(a - b)*sqrt(a + b)) + x*(3*a**4 + 4*a**2*b**2 + 8*b**4)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_96():
    f = cosh(x)**3/(a + b*sech(x))
    F = sinh(x)*cosh(x)**2/(3*a) - b*sinh(x)*cosh(x)/(2*a**2) + (2*a**2 + 3*b**2)*sinh(x)/(3*a**3) + 2*b**4*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**4*sqrt(a - b)*sqrt(a + b)) - b*x*(a**2 + 2*b**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_97():
    f = cosh(x)**2/(a + b*sech(x))
    F = sinh(x)*cosh(x)/(2*a) - b*sinh(x)/a**2 - 2*b**3*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**3*sqrt(a - b)*sqrt(a + b)) + x*(a**2 + 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_98():
    f = cosh(x)/(a + b*sech(x))
    F = sinh(x)/a + 2*b**2*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a**2*sqrt(a - b)*sqrt(a + b)) - b*x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_99():
    f = sech(x)/(a + b*sech(x))
    F = 2*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_100():
    f = sech(x)**2/(a + b*sech(x))
    F = -2*a*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b*sqrt(a - b)*sqrt(a + b)) + atan(sinh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_101():
    f = sech(x)**3/(a + b*sech(x))
    F = 2*a**2*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b**2*sqrt(a - b)*sqrt(a + b)) - a*atan(sinh(x))/b**2 + tanh(x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_102():
    f = sech(x)**4/(a + b*sech(x))
    F = -2*a**3*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(b**3*sqrt(a - b)*sqrt(a + b)) - a*tanh(x)/b**2 + tanh(x)*sech(x)/(2*b) + (2*a**2 + b**2)*atan(sinh(x))/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_103():
    f = tanh(x)**6/(a*sech(x) + a)
    F = x/a - (4 - 3*sech(x))*tanh(x)**3/(12*a) - (8 - 3*sech(x))*tanh(x)/(8*a) - 3*atan(sinh(x))/(8*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_104():
    f = tanh(x)**5/(a*sech(x) + a)
    F = log(cosh(x))/a - sech(x)**3/(3*a) + sech(x)**2/(2*a) + sech(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_105():
    f = tanh(x)**4/(a*sech(x) + a)
    F = x/a - (2 - sech(x))*tanh(x)/(2*a) - atan(sinh(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_106():
    f = tanh(x)**3/(a*sech(x) + a)
    F = log(cosh(x))/a + sech(x)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_107():
    f = tanh(x)**2/(a*sech(x) + a)
    F = x/a - atan(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_108():
    f = tanh(x)/(a*sech(x) + a)
    F = log(cosh(x) + 1)/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_109():
    f = coth(x)/(a*sech(x) + a)
    F = log(1 - cosh(x))/(4*a) + 3*log(cosh(x) + 1)/(4*a) + 1/(2*a*(cosh(x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_110():
    f = coth(x)**2/(a*sech(x) + a)
    F = x/a - (1 - sech(x))*coth(x)**3/(3*a) - (3 - 2*sech(x))*coth(x)/(3*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_111():
    f = coth(x)**3/(a*sech(x) + a)
    F = 5*log(1 - cosh(x))/(16*a) + 11*log(cosh(x) + 1)/(16*a) + 3/(4*a*(cosh(x) + 1)) - 1/(8*a*(cosh(x) + 1)**2) + 1/(8*a*(1 - cosh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_112():
    f = coth(x)**4/(a*sech(x) + a)
    F = x/a - (1 - sech(x))*coth(x)**5/(5*a) - (5 - 4*sech(x))*coth(x)**3/(15*a) - (15 - 8*sech(x))*coth(x)/(15*a)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_113():
    f = tanh(x)**7/(a + b*sech(x))
    F = -a*sech(x)**4/(4*b**2) - a*(a**2 - 3*b**2)*sech(x)**2/(2*b**4) + sech(x)**5/(5*b) + (a**2 - 3*b**2)*sech(x)**3/(3*b**3) + (a**4 - 3*a**2*b**2 + 3*b**4)*sech(x)/b**5 + log(cosh(x))/a - (a**2 - b**2)**3*log(a + b*sech(x))/(a*b**6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_114():
    f = tanh(x)**6/(a + b*sech(x))
    F = -a*tanh(x)**3/(3*b**2) + a*tanh(x)/b**2 + a*(a**2 - 3*b**2)*tanh(x)/b**4 - tanh(x)*sech(x)**3/(4*b) - 3*tanh(x)*sech(x)/(8*b) - 3*atan(sinh(x))/(8*b) - (a**2 - 3*b**2)*tanh(x)*sech(x)/(2*b**3) - (a**2 - 3*b**2)*atan(sinh(x))/(2*b**3) - (a**4 - 3*a**2*b**2 + 3*b**4)*atan(sinh(x))/b**5 + x/a + 2*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*b**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_115():
    f = tanh(x)**5/(a + b*sech(x))
    F = a*sech(x)**2/(2*b**2) - sech(x)**3/(3*b) - (a**2 - 2*b**2)*sech(x)/b**3 + log(cosh(x))/a + (a**2 - b**2)**2*log(a + b*sech(x))/(a*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_116():
    f = tanh(x)**4/(a + b*sech(x))
    F = -a*tanh(x)/b**2 + tanh(x)*sech(x)/(2*b) + (2*a**2 - 3*b**2)*atan(sinh(x))/(2*b**3) + x/a - 2*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_117():
    f = tanh(x)**3/(a + b*sech(x))
    F = sech(x)/b + (-a**2/b**2 + 1)*log(a + b*sech(x))/a + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_118():
    f = tanh(x)**2/(a + b*sech(x))
    F = -atan(sinh(x))/b + x/a + 2*sqrt(a - b)*sqrt(a + b)*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_119():
    f = tanh(x)/(a + b*sech(x))
    F = log(a + b*sech(x))/a + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_120():
    f = coth(x)/(a + b*sech(x))
    F = log(1 - sech(x))/(2*a + 2*b) + log(sech(x) + 1)/(2*a - 2*b) - b**2*log(a + b*sech(x))/(a*(a**2 - b**2)) + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_121():
    f = coth(x)**2/(a + b*sech(x))
    F = a*x/(a**2 - b**2) - a*coth(x)/(a**2 - b**2) + b*csch(x)/(a**2 - b**2) + 2*b**3*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - b**2*x/(a*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_122():
    f = coth(x)**3/(a + b*sech(x))
    F = -1/((4*a - 4*b)*(sech(x) + 1)) + (2*a + 3*b)*log(1 - sech(x))/(4*(a + b)**2) + (2*a - 3*b)*log(sech(x) + 1)/(4*(a - b)**2) - 1/((1 - sech(x))*(4*a + 4*b)) + b**4*log(a + b*sech(x))/(a*(a**2 - b**2)**2) + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_123():
    f = coth(x)**4/(a + b*sech(x))
    F = -a*b**2*x/(a**2 - b**2)**2 + a*b**2*coth(x)/(a**2 - b**2)**2 + a*x/(a**2 - b**2) - a*coth(x)**3/(3*a**2 - 3*b**2) - a*coth(x)/(a**2 - b**2) - b**3*csch(x)/(a**2 - b**2)**2 + b*csch(x)**3/(3*a**2 - 3*b**2) + b*csch(x)/(a**2 - b**2) - 2*b**5*atan(sqrt(a - b)*tanh(x/2)/sqrt(a + b))/(a*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + b**4*x/(a*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_124():
    f = coth(x)**5/(a + b*sech(x))
    F = -1/((16*a - 16*b)*(sech(x) + 1)**2) + (8*a**2 + 21*a*b + 15*b**2)*log(1 - sech(x))/(16*(a + b)**3) - (5*a - 7*b)/(16*(a - b)**2*(sech(x) + 1)) + (8*a**2 - 21*a*b + 15*b**2)*log(sech(x) + 1)/(16*(a - b)**3) - (5*a + 7*b)/(16*(1 - sech(x))*(a + b)**2) - 1/((1 - sech(x))**2*(16*a + 16*b)) - b**6*log(a + b*sech(x))/(a*(a**2 - b**2)**3) + log(cosh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_125():
    f = sqrt(a + b*sech(c + d*x))*tanh(c + d*x)**5
    F = 2*sqrt(a)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/d + 6*a*(a + b*sech(c + d*x))**(sympy.S(7)/2)/(7*b**4*d) + 2*a*(a + b*sech(c + d*x))**(sympy.S(3)/2)*(a**2 - 2*b**2)/(3*b**4*d) - 2*sqrt(a + b*sech(c + d*x))/d - 2*(a + b*sech(c + d*x))**(sympy.S(9)/2)/(9*b**4*d) - (a + b*sech(c + d*x))**(sympy.S(5)/2)*(6*a**2 - 4*b**2)/(5*b**4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_126():
    f = sqrt(a + b*sech(c + d*x))*tanh(c + d*x)**3
    F = 2*sqrt(a)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/d - 2*a*(a + b*sech(c + d*x))**(sympy.S(3)/2)/(3*b**2*d) - 2*sqrt(a + b*sech(c + d*x))/d + 2*(a + b*sech(c + d*x))**(sympy.S(5)/2)/(5*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_127():
    f = sqrt(a + b*sech(c + d*x))*tanh(c + d*x)
    F = 2*sqrt(a)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/d - 2*sqrt(a + b*sech(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_128():
    f = sqrt(a + b*sech(c + d*x))*coth(c + d*x)
    F = 2*sqrt(a)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/d - sqrt(a - b)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/d - sqrt(a + b)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_129():
    f = sqrt(a + b*sech(c + d*x))*coth(c + d*x)**3
    F = 2*sqrt(a)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/d - a*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) - a*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(d*sqrt(a - b)) - 3*b*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(4*d*sqrt(a + b)) + 3*b*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(4*d*sqrt(a - b)) - sqrt(a + b*sech(c + d*x))*coth(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_130():
    f = sqrt(a + b*sech(c + d*x))*tanh(c + d*x)**2
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_131():
    f = sqrt(a + b*sech(c + d*x))
    F = (Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Symbol('a') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_132():
    f = sqrt(a + b*sech(c + d*x))*coth(c + d*x)**2
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Symbol('a') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_133():
    f = tanh(c + d*x)**5/sqrt(a + b*sech(c + d*x))
    F = 6*a*(a + b*sech(c + d*x))**(sympy.S(5)/2)/(5*b**4*d) + 2*a*sqrt(a + b*sech(c + d*x))*(a**2 - 2*b**2)/(b**4*d) - 2*(a + b*sech(c + d*x))**(sympy.S(7)/2)/(7*b**4*d) - (a + b*sech(c + d*x))**(sympy.S(3)/2)*(6*a**2 - 4*b**2)/(3*b**4*d) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_134():
    f = tanh(c + d*x)**3/sqrt(a + b*sech(c + d*x))
    F = -2*a*sqrt(a + b*sech(c + d*x))/(b**2*d) + 2*(a + b*sech(c + d*x))**(sympy.S(3)/2)/(3*b**2*d) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_135():
    f = tanh(c + d*x)/sqrt(a + b*sech(c + d*x))
    F = 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_136():
    f = coth(c + d*x)/sqrt(a + b*sech(c + d*x))
    F = -atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) - atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(d*sqrt(a - b)) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_137():
    f = coth(c + d*x)**3/sqrt(a + b*sech(c + d*x))
    F = -b*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(4*d*(a + b)**(sympy.S(3)/2)) + b*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(4*d*(a - b)**(sympy.S(3)/2)) - sqrt(a + b*sech(c + d*x))/(d*(4*a - 4*b)*(sech(c + d*x) + 1)) - atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) - atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(d*sqrt(a - b)) - sqrt(a + b*sech(c + d*x))/(d*(1 - sech(c + d*x))*(4*a + 4*b)) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_138():
    f = tanh(c + d*x)**4/sqrt(a + b*sech(c + d*x))
    F = (Integer(-1) * ((Integer(4) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(9) * (Symbol('b'))**(Integer(2)))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b'))) + (Integer(9) * (Symbol('b'))**(Integer(2)))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(8) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sech((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_139():
    f = tanh(c + d*x)**2/sqrt(a + b*sech(c + d*x))
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_140():
    f = 1/sqrt(a + b*sech(c + d*x))
    F = (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_141():
    f = coth(c + d*x)**2/sqrt(a + b*sech(c + d*x))
    F = ((sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.coth((Symbol('c') + (Symbol('d') * x))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_142():
    f = tanh(c + d*x)**5/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = 2*a*(a + b*sech(c + d*x))**(sympy.S(3)/2)/(b**4*d) - 2*(a + b*sech(c + d*x))**(sympy.S(5)/2)/(5*b**4*d) - sqrt(a + b*sech(c + d*x))*(6*a**2 - 4*b**2)/(b**4*d) - 2*(a**2 - b**2)**2/(a*b**4*d*sqrt(a + b*sech(c + d*x))) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_143():
    f = tanh(c + d*x)**3/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(a + b*sech(c + d*x))/(b**2*d) + (2*a**2 - 2*b**2)/(a*b**2*d*sqrt(a + b*sech(c + d*x))) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_144():
    f = tanh(c + d*x)/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = -2/(a*d*sqrt(a + b*sech(c + d*x))) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_145():
    f = coth(c + d*x)/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(d*(a + b)**(sympy.S(3)/2)) - atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(d*(a - b)**(sympy.S(3)/2)) + 2*b**2/(a*d*sqrt(a + b*sech(c + d*x))*(a**2 - b**2)) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_146():
    f = coth(c + d*x)**3/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = -b*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(4*d*(a + b)**(sympy.S(5)/2)) + b*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(4*d*(a - b)**(sympy.S(5)/2)) - (2*a + 3*b)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a + b))/(2*d*(a + b)**(sympy.S(5)/2)) - sqrt(a + b*sech(c + d*x))/(4*d*(a - b)**2*(sech(c + d*x) + 1)) - (2*a - 3*b)*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a - b))/(2*d*(a - b)**(sympy.S(5)/2)) - sqrt(a + b*sech(c + d*x))/(4*d*(1 - sech(c + d*x))*(a + b)**2) - 2*b**4/(a*d*sqrt(a + b*sech(c + d*x))*(a**2 - b**2)**2) + 2*atanh(sqrt(a + b*sech(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_147():
    f = tanh(c + d*x)**4/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((Integer(4) * Symbol('a') * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((Integer(4) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(2) * Symbol('a')) + Symbol('b')) * ((Integer(4) * Symbol('a')) + Symbol('b')) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * Symbol('a') * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sech((Symbol('c') + (Symbol('d') * x))) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(2) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_148():
    f = tanh(c + d*x)**2/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_149():
    f = (a + b*sech(c + d*x))**(sympy.S(-3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_150():
    f = coth(c + d*x)**2/(a + b*sech(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(4) * Symbol('a') * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * Symbol('a')) + (Integer(-1) * Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.coth((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sech((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sech((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.coth((Symbol('c') + (Symbol('d') * x))) * ((Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tanh((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sech((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_151():
    f = (sech(a*c + b*c*x)**2)**(sympy.S(7)/2)*exp(c*(a + b*x))
    F = -64*sqrt(sech(a*c + b*c*x)**2)*cosh(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**3) + 48*sqrt(sech(a*c + b*c*x)**2)*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**4) - 192*sqrt(sech(a*c + b*c*x)**2)*cosh(a*c + b*c*x)/(5*b*c*(exp(2*c*(a + b*x)) + 1)**5) + 32*sqrt(sech(a*c + b*c*x)**2)*cosh(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_152():
    f = (sech(a*c + b*c*x)**2)**(sympy.S(5)/2)*exp(c*(a + b*x))
    F = -8*sqrt(sech(a*c + b*c*x)**2)*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**2) + 32*sqrt(sech(a*c + b*c*x)**2)*cosh(a*c + b*c*x)/(3*b*c*(exp(2*c*(a + b*x)) + 1)**3) - 4*sqrt(sech(a*c + b*c*x)**2)*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_153():
    f = (sech(a*c + b*c*x)**2)**(sympy.S(3)/2)*exp(c*(a + b*x))
    F = 2*sqrt(sech(a*c + b*c*x)**2)*exp(4*c*(a + b*x))*cosh(a*c + b*c*x)/(b*c*(exp(2*c*(a + b*x)) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_154():
    f = sqrt(sech(a*c + b*c*x)**2)*exp(c*(a + b*x))
    F = sqrt(sech(a*c + b*c*x)**2)*log(exp(2*c*(a + b*x)) + 1)*cosh(a*c + b*c*x)/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_155():
    f = exp(c*(a + b*x))/sqrt(sech(a*c + b*c*x)**2)
    F = x*sech(a*c + b*c*x)/(2*sqrt(sech(a*c + b*c*x)**2)) + exp(2*c*(a + b*x))*sech(a*c + b*c*x)/(4*b*c*sqrt(sech(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_156():
    f = exp(c*(a + b*x))/(sech(a*c + b*c*x)**2)**(sympy.S(3)/2)
    F = 3*x*sech(a*c + b*c*x)/(8*sqrt(sech(a*c + b*c*x)**2)) + exp(4*c*(a + b*x))*sech(a*c + b*c*x)/(32*b*c*sqrt(sech(a*c + b*c*x)**2)) + 3*exp(2*c*(a + b*x))*sech(a*c + b*c*x)/(16*b*c*sqrt(sech(a*c + b*c*x)**2)) - exp(-2*c*(a + b*x))*sech(a*c + b*c*x)/(16*b*c*sqrt(sech(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_157():
    f = exp(c*(a + b*x))/(sech(a*c + b*c*x)**2)**(sympy.S(5)/2)
    F = 5*x*sech(a*c + b*c*x)/(16*sqrt(sech(a*c + b*c*x)**2)) + exp(6*c*(a + b*x))*sech(a*c + b*c*x)/(192*b*c*sqrt(sech(a*c + b*c*x)**2)) + 5*exp(4*c*(a + b*x))*sech(a*c + b*c*x)/(128*b*c*sqrt(sech(a*c + b*c*x)**2)) + 5*exp(2*c*(a + b*x))*sech(a*c + b*c*x)/(32*b*c*sqrt(sech(a*c + b*c*x)**2)) - 5*exp(-2*c*(a + b*x))*sech(a*c + b*c*x)/(64*b*c*sqrt(sech(a*c + b*c*x)**2)) - exp(-4*c*(a + b*x))*sech(a*c + b*c*x)/(128*b*c*sqrt(sech(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_158():
    f = x**5/sqrt(sech(2*log(c*x)))
    F = x**6/(7*sqrt(sech(2*log(c*x)))) + 2*x**2/(21*c**4*sqrt(sech(2*log(c*x)))) + sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)/(21*c**5*x*(c**4 + x**(-4))*sqrt(sech(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_159():
    f = x**4/sqrt(sech(2*log(c*x)))
    F = x**5*(c**4 + x**(-4))/(6*c**4*sqrt(sech(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_160():
    f = x**3/sqrt(sech(2*log(c*x)))
    F = x**4/(5*sqrt(sech(2*log(c*x)))) + 2*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_e(2*acot(c*x), sympy.S.Half)/(5*c**3*x*(c**4 + x**(-4))*sqrt(sech(2*log(c*x)))) - sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)/(5*c**3*x*(c**4 + x**(-4))*sqrt(sech(2*log(c*x)))) + 2/(5*c**4*sqrt(sech(2*log(c*x)))) - 2/(5*c**4*x**2*(c**2 + x**(-2))*sqrt(sech(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_161():
    f = x**2/sqrt(sech(2*log(c*x)))
    F = x**3/(4*sqrt(sech(2*log(c*x)))) + atanh(sqrt(1 + 1/(c**4*x**4)))/(4*c**4*x*sqrt(1 + 1/(c**4*x**4))*sqrt(sech(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_162():
    f = x/sqrt(sech(2*log(c*x)))
    F = x**2/(3*sqrt(sech(2*log(c*x)))) - sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)/(3*c*x*(c**4 + x**(-4))*sqrt(sech(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_163():
    f = 1/sqrt(sech(2*log(c*x)))
    F = x/(2*sqrt(sech(2*log(c*x)))) - acsch(c**2*x**2)/(2*c**2*x*sqrt(1 + 1/(c**4*x**4))*sqrt(sech(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_164():
    f = sqrt(sech(2*log(c*x)))/x
    F = -I*sqrt(cosh(2*log(c*x)))*elliptic_f(I*log(c*x), 2)*sqrt(sech(2*log(c*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_165():
    f = sqrt(sech(2*log(c*x)))/x**2
    F = -c**2*x*sqrt(1 + 1/(c**4*x**4))*acsch(c**2*x**2)*sqrt(sech(2*log(c*x)))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_166():
    f = sqrt(sech(2*log(c*x)))/x**3
    F = c*x*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_e(2*acot(c*x), sympy.S.Half)*sqrt(sech(2*log(c*x))) - c*x*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)*sqrt(sech(2*log(c*x)))/2 - (c**4 + x**(-4))*sqrt(sech(2*log(c*x)))/(c**2 + x**(-2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_167():
    f = sqrt(sech(2*log(c*x)))/x**4
    F = x*(-c**4/2 - 1/(2*x**4))*sqrt(sech(2*log(c*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_168():
    f = sqrt(sech(2*log(c*x)))/x**5
    F = c**3*x*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)*sqrt(sech(2*log(c*x)))/6 + (-c**4/3 - 1/(3*x**4))*sqrt(sech(2*log(c*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_169():
    f = x**8/sech(2*log(c*x))**(sympy.S(3)/2)
    F = x**9/(12*sech(2*log(c*x))**(sympy.S(3)/2)) + x**5/((16*c**4 + 16/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) + x/(32*c**4*(c**4 + x**(-4))*sech(2*log(c*x))**(sympy.S(3)/2)) - atanh(sqrt(1 + 1/(c**4*x**4)))/(32*c**12*x**3*(1 + 1/(c**4*x**4))**(sympy.S(3)/2)*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_170():
    f = x**7/sech(2*log(c*x))**(sympy.S(3)/2)
    F = x**8/(11*sech(2*log(c*x))**(sympy.S(3)/2)) + 6*x**4/((77*c**4 + 77/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) + 4/(77*c**4*(c**4 + x**(-4))*sech(2*log(c*x))**(sympy.S(3)/2)) + 2*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)/(77*c**5*x**3*(c**4 + x**(-4))**2*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_171():
    f = x**6/sech(2*log(c*x))**(sympy.S(3)/2)
    F = x**7*(c**4 + x**(-4))/(10*c**4*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_172():
    f = x**5/sech(2*log(c*x))**(sympy.S(3)/2)
    F = x**6/(9*sech(2*log(c*x))**(sympy.S(3)/2)) + 2*x**2/((15*c**4 + 15/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) + 4*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_e(2*acot(c*x), sympy.S.Half)/(15*c**3*x**3*(c**4 + x**(-4))**2*sech(2*log(c*x))**(sympy.S(3)/2)) - 2*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)/(15*c**3*x**3*(c**4 + x**(-4))**2*sech(2*log(c*x))**(sympy.S(3)/2)) + 4/(15*c**4*x**2*(c**4 + x**(-4))*sech(2*log(c*x))**(sympy.S(3)/2)) - 4/(15*c**4*x**4*(c**2 + x**(-2))*(c**4 + x**(-4))*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_173():
    f = x**4/sech(2*log(c*x))**(sympy.S(3)/2)
    F = x**5/(8*sech(2*log(c*x))**(sympy.S(3)/2)) + 3*x/((16*c**4 + 16/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) + 3*atanh(sqrt(1 + 1/(c**4*x**4)))/(16*c**8*x**3*(1 + 1/(c**4*x**4))**(sympy.S(3)/2)*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_174():
    f = x**3/sech(2*log(c*x))**(sympy.S(3)/2)
    F = x**4/(7*sech(2*log(c*x))**(sympy.S(3)/2)) + 2/((7*c**4 + 7/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) - 2*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)/(7*c*x**3*(c**4 + x**(-4))**2*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_175():
    f = x**2/sech(2*log(c*x))**(sympy.S(3)/2)
    F = x**3/(6*sech(2*log(c*x))**(sympy.S(3)/2)) + 1/(x*(2*c**4 + 2/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) - acsch(c**2*x**2)/(2*c**6*x**3*(1 + 1/(c**4*x**4))**(sympy.S(3)/2)*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_176():
    f = x/sech(2*log(c*x))**(sympy.S(3)/2)
    F = 12*c*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_e(2*acot(c*x), sympy.S.Half)/(5*x**3*(c**4 + x**(-4))**2*sech(2*log(c*x))**(sympy.S(3)/2)) - 6*c*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*elliptic_f(2*acot(c*x), sympy.S.Half)/(5*x**3*(c**4 + x**(-4))**2*sech(2*log(c*x))**(sympy.S(3)/2)) + x**2/(5*sech(2*log(c*x))**(sympy.S(3)/2)) + 6/(x**2*(5*c**4 + 5/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) - 12/(x**4*(c**2 + x**(-2))*(5*c**4 + 5/x**4)*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_177():
    f = sech(2*log(c*x))**(sympy.S(-3)/2)
    F = x/(4*sech(2*log(c*x))**(sympy.S(3)/2)) - 3/(x**3*(4*c**4 + 4/x**4)*sech(2*log(c*x))**(sympy.S(3)/2)) + 3*atanh(sqrt(1 + 1/(c**4*x**4)))/(4*c**4*x**3*(1 + 1/(c**4*x**4))**(sympy.S(3)/2)*sech(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_178():
    f = sech(2*log(c*x))**(sympy.S(3)/2)/x
    F = sinh(2*log(c*x))*sqrt(sech(2*log(c*x))) + I*sqrt(cosh(2*log(c*x)))*elliptic_e(I*log(c*x), 2)*sqrt(sech(2*log(c*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_179():
    f = sech(2*log(c*x))**(sympy.S(3)/2)/x**2
    F = x**3*(c**4/2 + 1/(2*x**4))*sech(2*log(c*x))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_180():
    f = sech(2*log(c*x))**(sympy.S(3)/2)/x**3
    F = x**2*(c**4/2 + 1/(2*x**4))*sech(2*log(c*x))**(sympy.S(3)/2) - x**3*sqrt((c**4 + x**(-4))/(c**2 + x**(-2))**2)*(c**2 + x**(-2))*(c**4 + x**(-4))*elliptic_f(2*acot(c*x), sympy.S.Half)*sech(2*log(c*x))**(sympy.S(3)/2)/(4*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_181():
    f = sech(2*log(c*x))**(sympy.S(3)/2)/x**4
    F = -c**6*x**3*(1 + 1/(c**4*x**4))**(sympy.S(3)/2)*acsch(c**2*x**2)*sech(2*log(c*x))**(sympy.S(3)/2)/2 + x*(c**4/2 + 1/(2*x**4))*sech(2*log(c*x))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_182():
    f = sech(a + b*log(c*x**n))
    F = 2*x*(c*x**n)**b*exp(a)*hyper((1, (b + 1/n)/(2*b)), (sympy.S(3)/2 + 1/(2*b*n),), -(c*x**n)**(2*b)*exp(2*a))/(b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_183():
    f = sech(a + b*log(c*x**n))**2
    F = 4*x*(c*x**n)**(2*b)*exp(2*a)*hyper((2, 1 + 1/(2*b*n)), (2 + 1/(2*b*n),), -(c*x**n)**(2*b)*exp(2*a))/(2*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_184():
    f = sech(a + b*log(c*x**n))**3
    F = 8*x*(c*x**n)**(3*b)*exp(3*a)*hyper((3, (3*b + 1/n)/(2*b)), (sympy.S(5)/2 + 1/(2*b*n),), -(c*x**n)**(2*b)*exp(2*a))/(3*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_185():
    f = sech(a + b*log(c*x**n))**4
    F = 16*x*(c*x**n)**(4*b)*exp(4*a)*hyper((4, 2 + 1/(2*b*n)), (3 + 1/(2*b*n),), -(c*x**n)**(2*b)*exp(2*a))/(4*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_186():
    f = 2*b**2*n**2*sech(a + b*log(c*x**n))**3 + (-b**2*n**2 + 1)*sech(a + b*log(c*x**n))
    F = b*n*x*tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n)) + x*sech(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_187():
    f = sech(a + 2*log(c*sqrt(x)))**3
    F = 2*c**6*exp(-a)/(c**4 + exp(-2*a)/x**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_188():
    f = sech(a + 2*log(c/sqrt(x)))**3
    F = 2*c**2*exp(-3*a)/(c**4/x**2 + exp(-2*a))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_189():
    f = sech(a + log(c*x**n)/(n*(p - 2)))**p
    F = x*(2 - p)*((c*x**n)**(2/(n*(2 - p)))*exp(-2*a) + 1)*exp(2*a)*sech(a - log(c*x**n)/(n*(2 - p)))**p/((c*x**n)**(2/(n*(2 - p)))*(2 - 2*p))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_190():
    f = sech(a - log(c*x**n)/(n*(p - 2)))**p
    F = x*(1 + exp(-2*a)/(c*x**n)**(2/(n*(2 - p))))*(2 - p)*sech(a + log(c*x**n)/(n*(2 - p)))**p/(2 - 2*p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_191():
    f = sech(a + b*log(c*x**n))/x
    F = atan(sinh(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_192():
    f = sech(a + b*log(c*x**n))**2/x
    F = tanh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_193():
    f = sech(a + b*log(c*x**n))**3/x
    F = tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))/(2*b*n) + atan(sinh(a + b*log(c*x**n)))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_194():
    f = sech(a + b*log(c*x**n))**4/x
    F = -tanh(a + b*log(c*x**n))**3/(3*b*n) + tanh(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_195():
    f = sech(a + b*log(c*x**n))**5/x
    F = tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**3/(4*b*n) + 3*tanh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))/(8*b*n) + 3*atan(sinh(a + b*log(c*x**n)))/(8*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_196():
    f = sech(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = 2*sinh(a + b*log(c*x**n))*sech(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n) - 2*I*sqrt(cosh(a + b*log(c*x**n)))*elliptic_f(I*(a + b*log(c*x**n))/2, 2)*sqrt(sech(a + b*log(c*x**n)))/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_197():
    f = sech(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = 2*sinh(a + b*log(c*x**n))*sqrt(sech(a + b*log(c*x**n)))/(b*n) + 2*I*sqrt(cosh(a + b*log(c*x**n)))*elliptic_e(I*(a + b*log(c*x**n))/2, 2)*sqrt(sech(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_198():
    f = sqrt(sech(a + b*log(c*x**n)))/x
    F = -2*I*sqrt(cosh(a + b*log(c*x**n)))*elliptic_f(I*(a + b*log(c*x**n))/2, 2)*sqrt(sech(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_199():
    f = 1/(x*sqrt(sech(a + b*log(c*x**n))))
    F = -2*I*sqrt(cosh(a + b*log(c*x**n)))*elliptic_e(I*(a + b*log(c*x**n))/2, 2)*sqrt(sech(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_200():
    f = 1/(x*sech(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = 2*sinh(a + b*log(c*x**n))/(3*b*n*sqrt(sech(a + b*log(c*x**n)))) - 2*I*sqrt(cosh(a + b*log(c*x**n)))*elliptic_f(I*(a + b*log(c*x**n))/2, 2)*sqrt(sech(a + b*log(c*x**n)))/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_5_Hyperbolic_secant_6_5_3_Hyperbolic_secant_functions_201():
    f = 1/(x*sech(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 2*sinh(a + b*log(c*x**n))/(5*b*n*sech(a + b*log(c*x**n))**(sympy.S(3)/2)) - 6*I*sqrt(cosh(a + b*log(c*x**n)))*elliptic_e(I*(a + b*log(c*x**n))/2, 2)*sqrt(sech(a + b*log(c*x**n)))/(5*b*n)
    assert integrate(f, x) == F

