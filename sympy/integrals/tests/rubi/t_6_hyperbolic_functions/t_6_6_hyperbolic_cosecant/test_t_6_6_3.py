"""Generated from MathematicaSyntaxTestSuite.

Source: 6 Hyperbolic functions/6.6 Hyperbolic cosecant/6.6.3 Hyperbolic cosecant functions.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, n, p = symbols('a b c d n p')

def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_1():
    f = csch(a + b*x)
    F = -atanh(cosh(a + b*x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_2():
    f = csch(a + b*x)**2
    F = -coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_3():
    f = csch(a + b*x)**3
    F = -coth(a + b*x)*csch(a + b*x)/(2*b) + atanh(cosh(a + b*x))/(2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_4():
    f = csch(a + b*x)**4
    F = -coth(a + b*x)**3/(3*b) + coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_5():
    f = csch(a + b*x)**5
    F = -coth(a + b*x)*csch(a + b*x)**3/(4*b) + 3*coth(a + b*x)*csch(a + b*x)/(8*b) - 3*atanh(cosh(a + b*x))/(8*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_6():
    f = csch(a + b*x)**6
    F = -coth(a + b*x)**5/(5*b) + 2*coth(a + b*x)**3/(3*b) - coth(a + b*x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_7():
    f = csch(a + b*x)**(sympy.S(5)/2)
    F = 2*I*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(3*b) - 2*cosh(a + b*x)*csch(a + b*x)**(sympy.S(3)/2)/(3*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_8():
    f = csch(a + b*x)**(sympy.S(3)/2)
    F = -2*cosh(a + b*x)*sqrt(csch(a + b*x))/b - 2*I*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(b*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_9():
    f = sqrt(csch(a + b*x))
    F = -2*I*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_10():
    f = 1/sqrt(csch(a + b*x))
    F = -2*I*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(b*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_11():
    f = csch(a + b*x)**(sympy.S(-3)/2)
    F = 2*I*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x))*elliptic_f(I*a/2 + I*b*x/2 - pi/4, 2)/(3*b) + 2*cosh(a + b*x)/(3*b*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_12():
    f = csch(a + b*x)**(sympy.S(-5)/2)
    F = 2*cosh(a + b*x)/(5*b*csch(a + b*x)**(sympy.S(3)/2)) + 6*I*elliptic_e(I*a/2 + I*b*x/2 - pi/4, 2)/(5*b*sqrt(I*sinh(a + b*x))*sqrt(csch(a + b*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_13():
    f = (b*csch(c + d*x))**(sympy.S(7)/2)
    F = 6*I*b**4*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(5*d*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x))) + 6*b**3*sqrt(b*csch(c + d*x))*cosh(c + d*x)/(5*d) - 2*b*(b*csch(c + d*x))**(sympy.S(5)/2)*cosh(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_14():
    f = (b*csch(c + d*x))**(sympy.S(5)/2)
    F = 2*I*b**2*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(3*d) - 2*b*(b*csch(c + d*x))**(sympy.S(3)/2)*cosh(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_15():
    f = (b*csch(c + d*x))**(sympy.S(3)/2)
    F = -2*I*b**2*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(d*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x))) - 2*b*sqrt(b*csch(c + d*x))*cosh(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_16():
    f = sqrt(b*csch(c + d*x))
    F = -2*I*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_17():
    f = 1/sqrt(b*csch(c + d*x))
    F = -2*I*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(d*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_18():
    f = (b*csch(c + d*x))**(sympy.S(-3)/2)
    F = 2*cosh(c + d*x)/(3*b*d*sqrt(b*csch(c + d*x))) + 2*I*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_19():
    f = (b*csch(c + d*x))**(sympy.S(-5)/2)
    F = 2*cosh(c + d*x)/(5*b*d*(b*csch(c + d*x))**(sympy.S(3)/2)) + 6*I*elliptic_e(I*c/2 + I*d*x/2 - pi/4, 2)/(5*b**2*d*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_20():
    f = (b*csch(c + d*x))**(sympy.S(-7)/2)
    F = 2*cosh(c + d*x)/(7*b*d*(b*csch(c + d*x))**(sympy.S(5)/2)) - 10*cosh(c + d*x)/(21*b**3*d*sqrt(b*csch(c + d*x))) - 10*I*sqrt(I*sinh(c + d*x))*sqrt(b*csch(c + d*x))*elliptic_f(I*c/2 + I*d*x/2 - pi/4, 2)/(21*b**4*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_21():
    f = (b*csch(c + d*x))**n
    F = b*(b*csch(c + d*x))**(n - 1)*cosh(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), -sinh(c + d*x)**2)/(d*(1 - n)*sqrt(cosh(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_22():
    f = (-csch(x)**2)**(sympy.S(5)/2)
    F = (-csch(x)**2)**(sympy.S(3)/2)*coth(x)/4 + 3*sqrt(-csch(x)**2)*coth(x)/8 + 3*asin(coth(x))/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_23():
    f = (-csch(x)**2)**(sympy.S(3)/2)
    F = sqrt(-csch(x)**2)*coth(x)/2 + asin(coth(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_24():
    f = sqrt(-csch(x)**2)
    F = asin(coth(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_25():
    f = 1/sqrt(-csch(x)**2)
    F = coth(x)/sqrt(-csch(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_26():
    f = (-csch(x)**2)**(sympy.S(-3)/2)
    F = 2*coth(x)/(3*sqrt(-csch(x)**2)) + coth(x)/(3*(-csch(x)**2)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_27():
    f = (-csch(x)**2)**(sympy.S(-5)/2)
    F = 8*coth(x)/(15*sqrt(-csch(x)**2)) + 4*coth(x)/(15*(-csch(x)**2)**(sympy.S(3)/2)) + coth(x)/(5*(-csch(x)**2)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_28():
    f = (-csch(x)**2)**(sympy.S(-7)/2)
    F = 16*coth(x)/(35*sqrt(-csch(x)**2)) + 8*coth(x)/(35*(-csch(x)**2)**(sympy.S(3)/2)) + 6*coth(x)/(35*(-csch(x)**2)**(sympy.S(5)/2)) + coth(x)/(7*(-csch(x)**2)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_29():
    f = (a*csch(x)**2)**(sympy.S(5)/2)
    F = -3*a**(sympy.S(5)/2)*atanh(sqrt(a)*coth(x)/sqrt(a*csch(x)**2))/8 + 3*a**2*sqrt(a*csch(x)**2)*coth(x)/8 - a*(a*csch(x)**2)**(sympy.S(3)/2)*coth(x)/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_30():
    f = (a*csch(x)**2)**(sympy.S(3)/2)
    F = a**(sympy.S(3)/2)*atanh(sqrt(a)*coth(x)/sqrt(a*csch(x)**2))/2 - a*sqrt(a*csch(x)**2)*coth(x)/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_31():
    f = sqrt(a*csch(x)**2)
    F = -sqrt(a)*atanh(sqrt(a)*coth(x)/sqrt(a*csch(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_32():
    f = 1/sqrt(a*csch(x)**2)
    F = coth(x)/sqrt(a*csch(x)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_33():
    f = (a*csch(x)**2)**(sympy.S(-3)/2)
    F = coth(x)/(3*(a*csch(x)**2)**(sympy.S(3)/2)) - 2*coth(x)/(3*a*sqrt(a*csch(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_34():
    f = (a*csch(x)**2)**(sympy.S(-5)/2)
    F = coth(x)/(5*(a*csch(x)**2)**(sympy.S(5)/2)) - 4*coth(x)/(15*a*(a*csch(x)**2)**(sympy.S(3)/2)) + 8*coth(x)/(15*a**2*sqrt(a*csch(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_35():
    f = (a*csch(x)**2)**(sympy.S(-7)/2)
    F = coth(x)/(7*(a*csch(x)**2)**(sympy.S(7)/2)) - 6*coth(x)/(35*a*(a*csch(x)**2)**(sympy.S(5)/2)) + 8*coth(x)/(35*a**2*(a*csch(x)**2)**(sympy.S(3)/2)) - 16*coth(x)/(35*a**3*sqrt(a*csch(x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_36():
    f = (a*csch(x)**3)**(sympy.S(5)/2)
    F = 154*a**2*sqrt(a*csch(x)**3)*sinh(x)*cosh(x)/195 - 2*a**2*sqrt(a*csch(x)**3)*coth(x)*csch(x)**4/13 + 22*a**2*sqrt(a*csch(x)**3)*coth(x)*csch(x)**2/117 - 154*a**2*sqrt(a*csch(x)**3)*coth(x)/585 + 154*I*a**2*sqrt(a*csch(x)**3)*sinh(x)**2*elliptic_e(I*x/2 - pi/4, 2)/(195*sqrt(I*sinh(x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_37():
    f = (a*csch(x)**3)**(sympy.S(3)/2)
    F = -10*I*a*sqrt(I*sinh(x))*sqrt(a*csch(x)**3)*sinh(x)*elliptic_f(I*x/2 - pi/4, 2)/21 + 10*a*sqrt(a*csch(x)**3)*cosh(x)/21 - 2*a*sqrt(a*csch(x)**3)*coth(x)*csch(x)/7
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_38():
    f = 1/sqrt(a*csch(x)**3)
    F = 2*I*sqrt(I*sinh(x))*csch(x)**2*elliptic_f(I*x/2 - pi/4, 2)/(3*sqrt(a*csch(x)**3)) + 2*coth(x)/(3*sqrt(a*csch(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_39():
    f = (a*csch(x)**3)**(sympy.S(-3)/2)
    F = 2*sinh(x)**2*cosh(x)/(9*a*sqrt(a*csch(x)**3)) - 14*cosh(x)/(45*a*sqrt(a*csch(x)**3)) - 14*I*csch(x)*elliptic_e(I*x/2 - pi/4, 2)/(15*a*sqrt(I*sinh(x))*sqrt(a*csch(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_40():
    f = (a*csch(x)**3)**(sympy.S(-5)/2)
    F = -26*I*sqrt(I*sinh(x))*csch(x)**2*elliptic_f(I*x/2 - pi/4, 2)/(77*a**2*sqrt(a*csch(x)**3)) + 2*sinh(x)**5*cosh(x)/(15*a**2*sqrt(a*csch(x)**3)) - 26*sinh(x)**3*cosh(x)/(165*a**2*sqrt(a*csch(x)**3)) + 78*sinh(x)*cosh(x)/(385*a**2*sqrt(a*csch(x)**3)) - 26*coth(x)/(77*a**2*sqrt(a*csch(x)**3))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_41():
    f = (a*csch(x)**4)**(sympy.S(7)/2)
    F = -a**3*sqrt(a*csch(x)**4)*sinh(x)*cosh(x) - a**3*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**11/13 + 6*a**3*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**9/11 - 5*a**3*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**7/3 + 20*a**3*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**5/7 - 3*a**3*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**3 + 2*a**3*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_42():
    f = (a*csch(x)**4)**(sympy.S(5)/2)
    F = -a**2*sqrt(a*csch(x)**4)*sinh(x)*cosh(x) - a**2*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**7/9 + 4*a**2*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**5/7 - 6*a**2*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**3/5 + 4*a**2*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_43():
    f = (a*csch(x)**4)**(sympy.S(3)/2)
    F = -a*sqrt(a*csch(x)**4)*sinh(x)*cosh(x) - a*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)**3/5 + 2*a*sqrt(a*csch(x)**4)*cosh(x)**2*coth(x)/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_44():
    f = sqrt(a*csch(x)**4)
    F = -sqrt(a*csch(x)**4)*sinh(x)*cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_45():
    f = 1/sqrt(a*csch(x)**4)
    F = -x*csch(x)**2/(2*sqrt(a*csch(x)**4)) + coth(x)/(2*sqrt(a*csch(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_46():
    f = (a*csch(x)**4)**(sympy.S(-3)/2)
    F = -5*x*csch(x)**2/(16*a*sqrt(a*csch(x)**4)) + sinh(x)**3*cosh(x)/(6*a*sqrt(a*csch(x)**4)) - 5*sinh(x)*cosh(x)/(24*a*sqrt(a*csch(x)**4)) + 5*coth(x)/(16*a*sqrt(a*csch(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_47():
    f = (a*csch(x)**4)**(sympy.S(-5)/2)
    F = -63*x*csch(x)**2/(256*a**2*sqrt(a*csch(x)**4)) + sinh(x)**7*cosh(x)/(10*a**2*sqrt(a*csch(x)**4)) - 9*sinh(x)**5*cosh(x)/(80*a**2*sqrt(a*csch(x)**4)) + 21*sinh(x)**3*cosh(x)/(160*a**2*sqrt(a*csch(x)**4)) - 21*sinh(x)*cosh(x)/(128*a**2*sqrt(a*csch(x)**4)) + 63*coth(x)/(256*a**2*sqrt(a*csch(x)**4))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_48():
    f = 1/(I*a*csch(a + b*x) + a)
    F = -coth(a + b*x)/(b*(I*a*csch(a + b*x) + a)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_49():
    f = 1/(-I*a*csch(a + b*x) + a)
    F = -coth(a + b*x)/(b*(-I*a*csch(a + b*x) + a)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_50():
    f = (I*a*csch(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(5)/2)*atanh(sqrt(a)*coth(c + d*x)/sqrt(I*a*csch(c + d*x) + a))/d + 14*a**3*coth(c + d*x)/(3*d*sqrt(I*a*csch(c + d*x) + a)) + 2*a**2*sqrt(I*a*csch(c + d*x) + a)*coth(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_51():
    f = (I*a*csch(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(3)/2)*atanh(sqrt(a)*coth(c + d*x)/sqrt(I*a*csch(c + d*x) + a))/d + 2*a**2*coth(c + d*x)/(d*sqrt(I*a*csch(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_52():
    f = sqrt(I*a*csch(c + d*x) + a)
    F = 2*sqrt(a)*atanh(sqrt(a)*coth(c + d*x)/sqrt(I*a*csch(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_53():
    f = 1/sqrt(I*a*csch(c + d*x) + a)
    F = 2*atanh(sqrt(a)*coth(c + d*x)/sqrt(I*a*csch(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*coth(c + d*x)/(2*sqrt(I*a*csch(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_54():
    f = (I*a*csch(c + d*x) + a)**(sympy.S(-3)/2)
    F = -coth(c + d*x)/(2*d*(I*a*csch(c + d*x) + a)**(sympy.S(3)/2)) + 2*atanh(sqrt(a)*coth(c + d*x)/sqrt(I*a*csch(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*coth(c + d*x)/(2*sqrt(I*a*csch(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_55():
    f = sqrt(-I*a*csch(c + d*x) + a)
    F = 2*sqrt(a)*atanh(sqrt(a)*coth(c + d*x)/sqrt(-I*a*csch(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_56():
    f = 1/sqrt(-I*a*csch(c + d*x) + a)
    F = 2*atanh(sqrt(a)*coth(c + d*x)/sqrt(-I*a*csch(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*coth(c + d*x)/(2*sqrt(-I*a*csch(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_57():
    f = sqrt(3*I*csch(x) + 3)
    F = 2*sqrt(3)*atanh(coth(x)/sqrt(I*csch(x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_58():
    f = sqrt(-3*I*csch(x) + 3)
    F = 2*sqrt(3)*atanh(coth(x)/sqrt(-I*csch(x) + 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_59():
    f = sqrt(3*I*csch(x) - 3)
    F = -2*sqrt(3)*atan(coth(x)/sqrt(I*csch(x) - 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_60():
    f = sqrt(-3*I*csch(x) - 3)
    F = -2*sqrt(3)*atan(coth(x)/sqrt(-I*csch(x) - 1))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_61():
    f = sinh(x)**4/(csch(x) + I)
    F = -15*I*x/8 - 5*I*sinh(x)**3*cosh(x)/4 + 15*I*sinh(x)*cosh(x)/8 + 4*cosh(x)**3/3 - 4*cosh(x) - sinh(x)**3*cosh(x)/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_62():
    f = sinh(x)**3/(csch(x) + I)
    F = -3*x/2 + 3*sinh(x)*cosh(x)/2 - 4*I*cosh(x)**3/3 + 4*I*cosh(x) - sinh(x)**2*cosh(x)/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_63():
    f = sinh(x)**2/(csch(x) + I)
    F = 3*I*x/2 - 3*I*sinh(x)*cosh(x)/2 + 2*cosh(x) - sinh(x)*cosh(x)/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_64():
    f = sinh(x)/(csch(x) + I)
    F = x - 2*I*cosh(x) - cosh(x)/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_65():
    f = csch(x)/(csch(x) + I)
    F = I*coth(x)/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_66():
    f = csch(x)**2/(csch(x) + I)
    F = -atanh(cosh(x)) + coth(x)/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_67():
    f = csch(x)**3/(csch(x) + I)
    F = -coth(x) + I*atanh(cosh(x)) - I*coth(x)/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_68():
    f = csch(x)**4/(csch(x) + I)
    F = -3*coth(x)*csch(x)/2 + 2*I*coth(x) + 3*atanh(cosh(x))/2 + coth(x)*csch(x)**2/(csch(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_69():
    f = (a + b*csch(c + d*x))**4
    F = a**4*x - 4*a*b**3*coth(c + d*x)*csch(c + d*x)/(3*d) - 2*a*b*(2*a**2 - b**2)*atanh(cosh(c + d*x))/d - b**2*(a + b*csch(c + d*x))**2*coth(c + d*x)/(3*d) - b**2*(17*a**2 - 2*b**2)*coth(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_70():
    f = (a + b*csch(c + d*x))**3
    F = a**3*x - 5*a*b**2*coth(c + d*x)/(2*d) - b**2*(a + b*csch(c + d*x))*coth(c + d*x)/(2*d) - b*(6*a**2 - b**2)*atanh(cosh(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_71():
    f = (a + b*csch(c + d*x))**2
    F = a**2*x - 2*a*b*atanh(cosh(c + d*x))/d - b**2*coth(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_72():
    f = a + b*csch(c + d*x)
    F = a*x - b*atanh(cosh(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_73():
    f = 1/(a + b*csch(c + d*x))
    F = 2*b*atanh((a - b*tanh(c/2 + d*x/2))/sqrt(a**2 + b**2))/(a*d*sqrt(a**2 + b**2)) + x/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_74():
    f = (a + b*csch(c + d*x))**(-2)
    F = -b**2*coth(c + d*x)/(a*d*(a + b*csch(c + d*x))*(a**2 + b**2)) + 2*b*(2*a**2 + b**2)*atanh((a - b*tanh(c/2 + d*x/2))/sqrt(a**2 + b**2))/(a**2*d*(a**2 + b**2)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_75():
    f = (a + b*csch(c + d*x))**(-3)
    F = -b**2*coth(c + d*x)/(2*a*d*(a + b*csch(c + d*x))**2*(a**2 + b**2)) - b**2*(5*a**2 + 2*b**2)*coth(c + d*x)/(2*a**2*d*(a + b*csch(c + d*x))*(a**2 + b**2)**2) + b*(6*a**4 + 5*a**2*b**2 + 2*b**4)*atanh((a - b*tanh(c/2 + d*x/2))/sqrt(a**2 + b**2))/(a**3*d*(a**2 + b**2)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_76():
    f = sinh(x)**3/(a + b*csch(x))
    F = sinh(x)**2*cosh(x)/(3*a) - b*sinh(x)*cosh(x)/(2*a**2) - (2*a**2 - 3*b**2)*cosh(x)/(3*a**3) - 2*b**4*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a**4*sqrt(a**2 + b**2)) + b*x*(a**2 - 2*b**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_77():
    f = sinh(x)**2/(a + b*csch(x))
    F = sinh(x)*cosh(x)/(2*a) - b*cosh(x)/a**2 + 2*b**3*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a**3*sqrt(a**2 + b**2)) - x*(a**2 - 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_78():
    f = sinh(x)/(a + b*csch(x))
    F = cosh(x)/a - 2*b**2*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a**2*sqrt(a**2 + b**2)) - b*x/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_79():
    f = csch(x)/(a + b*csch(x))
    F = -2*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/sqrt(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_80():
    f = csch(x)**2/(a + b*csch(x))
    F = 2*a*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(b*sqrt(a**2 + b**2)) - atanh(cosh(x))/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_81():
    f = csch(x)**3/(a + b*csch(x))
    F = -2*a**2*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(b**2*sqrt(a**2 + b**2)) + a*atanh(cosh(x))/b**2 - coth(x)/b
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_82():
    f = csch(x)**4/(a + b*csch(x))
    F = 2*a**3*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(b**3*sqrt(a**2 + b**2)) + a*coth(x)/b**2 - coth(x)*csch(x)/(2*b) - (2*a**2 - b**2)*atanh(cosh(x))/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_83():
    f = cosh(x)**4/(csch(x) + I)
    F = I*x/8 - I*sinh(x)*cosh(x)**3/4 + I*sinh(x)*cosh(x)/8 + cosh(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_84():
    f = cosh(x)**3/(csch(x) + I)
    F = -I*sinh(x)**3/3 + sinh(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_85():
    f = cosh(x)**2/(csch(x) + I)
    F = I*x/2 - I*sinh(x)*cosh(x)/2 + cosh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_86():
    f = cosh(x)/(csch(x) + I)
    F = log(-sinh(x) + I) - I*sinh(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_87():
    f = sech(x)/(csch(x) + I)
    F = I*tanh(x)*sech(x)/2 - I*atan(sinh(x))/2 - sech(x)**2/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_88():
    f = sech(x)**2/(csch(x) + I)
    F = -I*tanh(x)**3/3 - sech(x)**3/3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_89():
    f = sech(x)**3/(csch(x) + I)
    F = I*tanh(x)*sech(x)**3/4 - I*tanh(x)*sech(x)/8 - I*atan(sinh(x))/8 - sech(x)**4/4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_90():
    f = sech(x)**4/(csch(x) + I)
    F = I*tanh(x)**5/5 - I*tanh(x)**3/3 - sech(x)**5/5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_91():
    f = cosh(x)**5/(a + b*csch(x))
    F = sinh(x)**5/(5*a) - b*sinh(x)**4/(4*a**2) + (2*a**2 + b**2)*sinh(x)**3/(3*a**3) - b*(2*a**2 + b**2)*sinh(x)**2/(2*a**4) + (a**2 + b**2)**2*sinh(x)/a**5 - b*(a**2 + b**2)**2*log(a*sinh(x) + b)/a**6
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_92():
    f = cosh(x)**4/(a + b*csch(x))
    F = -(-3*a*sinh(x) + 4*b)*cosh(x)**3/(12*a**2) - (-a*(3*a**2 + 4*b**2)*sinh(x) + 8*b*(a**2 + b**2))*cosh(x)/(8*a**4) + 2*b*(a**2 + b**2)**(sympy.S(3)/2)*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/a**5 + x*(3*a**4 + 12*a**2*b**2 + 8*b**4)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_93():
    f = cosh(x)**3/(a + b*csch(x))
    F = sinh(x)**3/(3*a) - b*sinh(x)**2/(2*a**2) + (a**2 + b**2)*sinh(x)/a**3 - b*(a**2 + b**2)*log(a*sinh(x) + b)/a**4
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_94():
    f = cosh(x)**2/(a + b*csch(x))
    F = -(-a*sinh(x) + 2*b)*cosh(x)/(2*a**2) + 2*b*sqrt(a**2 + b**2)*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/a**3 + x*(a**2 + 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_95():
    f = cosh(x)/(a + b*csch(x))
    F = sinh(x)/a - b*log(a*sinh(x) + b)/a**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_96():
    f = sech(x)/(a + b*csch(x))
    F = -b*log(a*sinh(x) + b)/(a**2 + b**2) + log(-sinh(x) + I)/(2*I*a + 2*b) - log(sinh(x) + I)/(2*I*a - 2*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_97():
    f = sech(x)**2/(a + b*csch(x))
    F = 2*a*b*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(3)/2) - (-a*sinh(x) + b)*sech(x)/(a**2 + b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_98():
    f = sech(x)**3/(a + b*csch(x))
    F = -a**2*b*log(a*sinh(x) + b)/(a**2 + b**2)**2 + I*a*log(sinh(x) + I)/(4*(a + I*b)**2) - I*a*log(-sinh(x) + I)/(4*(a - I*b)**2) - (-a*sinh(x) + b)*sech(x)**2/(2*a**2 + 2*b**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_99():
    f = sech(x)**4/(a + b*csch(x))
    F = 2*a**3*b*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a**2 + b**2)**(sympy.S(5)/2) - (-a*sinh(x) + b)*sech(x)**3/(3*a**2 + 3*b**2) - (3*a**2*b - a*(2*a**2 - b**2)*sinh(x))*sech(x)/(3*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_100():
    f = sech(x)**5/(a + b*csch(x))
    F = -a**4*b*log(a*sinh(x) + b)/(a**2 + b**2)**3 + a*(3*a + I*b)*log(sinh(x) + I)/(16*(I*a - b)**3) - a*(3*I*a + b)*log(-sinh(x) + I)/(16*(a - I*b)**3) - (-a*sinh(x) + b)*sech(x)**4/(4*a**2 + 4*b**2) - (4*a**2*b - a*(3*a**2 - b**2)*sinh(x))*sech(x)**2/(8*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_101():
    f = tanh(x)**5/(csch(x) + I)
    F = -21*I*log(-sinh(x) + I)/32 - 11*I*log(sinh(x) + I)/32 - 15*I/(16*I*sinh(x) + 16) + 9*I/(32*(I*sinh(x) + 1)**2) - I/(24*(I*sinh(x) + 1)**3) + I/(32*(-I*sinh(x) + 1)**2) - I/(-4*I*sinh(x) + 4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_102():
    f = tanh(x)**4/(csch(x) + I)
    F = -I*x + (-8*csch(x)/15 + I)*tanh(x) + (-4*csch(x)/15 + I/3)*tanh(x)**3 + (-csch(x)/5 + I/5)*tanh(x)**5
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_103():
    f = tanh(x)**3/(csch(x) + I)
    F = -11*I*log(-sinh(x) + I)/16 - 5*I*log(sinh(x) + I)/16 - 3*I/(4*I*sinh(x) + 4) + I/(8*(I*sinh(x) + 1)**2) - I/(-8*I*sinh(x) + 8)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_104():
    f = tanh(x)**2/(csch(x) + I)
    F = -I*x + (-2*csch(x)/3 + I)*tanh(x) + (-csch(x)/3 + I/3)*tanh(x)**3
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_105():
    f = tanh(x)/(csch(x) + I)
    F = -3*I*log(-sinh(x) + I)/4 - I*log(sinh(x) + I)/4 - I/(2*I*sinh(x) + 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_106():
    f = coth(x)/(csch(x) + I)
    F = -I*log(-sinh(x) + I)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_107():
    f = coth(x)**2/(csch(x) + I)
    F = -I*x - atanh(cosh(x))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_108():
    f = coth(x)**3/(csch(x) + I)
    F = -I*log(sinh(x)) - csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_109():
    f = coth(x)**4/(csch(x) + I)
    F = -I*x + (-csch(x) + 2*I)*coth(x)/2 - atanh(cosh(x))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_110():
    f = coth(x)**5/(csch(x) + I)
    F = -I*log(sinh(x)) - csch(x)**3/3 + I*csch(x)**2/2 - csch(x)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_111():
    f = coth(x)**6/(csch(x) + I)
    F = -I*x + (-3*csch(x) + 4*I)*coth(x)**3/12 + (-3*csch(x) + 8*I)*coth(x)/8 - 3*atanh(cosh(x))/8
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_112():
    f = tanh(x)**5/(a + b*csch(x))
    F = -a*(a**4 + 3*a**2*b**2 + 3*b**4)*log(tanh(x))/(a**2 + b**2)**3 - b**5*atan(sinh(x))/(a**2 + b**2)**3 - b**3*atan(sinh(x))/(2*(a**2 + b**2)**2) + 3*b*tanh(x)*sech(x)/(8*a**2 + 8*b**2) - 3*b*atan(sinh(x))/(8*a**2 + 8*b**2) - (a - b*csch(x))*tanh(x)**4/(4*a**2 + 4*b**2) - (a*(a**2 + 2*b**2) - b**3*csch(x))*tanh(x)**2/(2*(a**2 + b**2)**2) + b**6*log(a + b*csch(x))/(a*(a**2 + b**2)**3) + log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_113():
    f = tanh(x)**4/(a + b*csch(x))
    F = a*b**2*x/(a**2 + b**2)**2 - a*b**2*tanh(x)/(a**2 + b**2)**2 + a*x/(a**2 + b**2) - a*tanh(x)**3/(3*a**2 + 3*b**2) - a*tanh(x)/(a**2 + b**2) + b**3*sech(x)/(a**2 + b**2)**2 - b*sech(x)**3/(3*a**2 + 3*b**2) + b*sech(x)/(a**2 + b**2) + 2*b**5*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a*(a**2 + b**2)**(sympy.S(5)/2)) + b**4*x/(a*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_114():
    f = tanh(x)**3/(a + b*csch(x))
    F = -a*(a**2 + 2*b**2)*log(tanh(x))/(a**2 + b**2)**2 - b**3*atan(sinh(x))/(a**2 + b**2)**2 - b*atan(sinh(x))/(2*a**2 + 2*b**2) - (a - b*csch(x))*tanh(x)**2/(2*a**2 + 2*b**2) + b**4*log(a + b*csch(x))/(a*(a**2 + b**2)**2) + log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_115():
    f = tanh(x)**2/(a + b*csch(x))
    F = a*x/(a**2 + b**2) - a*tanh(x)/(a**2 + b**2) + b*sech(x)/(a**2 + b**2) + 2*b**3*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a*(a**2 + b**2)**(sympy.S(3)/2)) + b**2*x/(a*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_116():
    f = tanh(x)/(a + b*csch(x))
    F = -a*log(tanh(x))/(a**2 + b**2) - b*atan(sinh(x))/(a**2 + b**2) + b**2*log(a + b*csch(x))/(a*(a**2 + b**2)) + log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_117():
    f = coth(x)/(a + b*csch(x))
    F = log(a + b*csch(x))/a + log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_118():
    f = coth(x)**2/(a + b*csch(x))
    F = -atanh(cosh(x))/b + x/a + 2*sqrt(a**2 + b**2)*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a*b)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_119():
    f = coth(x)**3/(a + b*csch(x))
    F = (a/b**2 + 1/a)*log(a + b*csch(x)) - csch(x)/b + log(sinh(x))/a
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_120():
    f = coth(x)**4/(a + b*csch(x))
    F = a*coth(x)/b**2 - coth(x)*csch(x)/(2*b) - (2*a**2 + 3*b**2)*atanh(cosh(x))/(2*b**3) + x/a + 2*(a**2 + b**2)**(sympy.S(3)/2)*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a*b**3)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_121():
    f = coth(x)**5/(a + b*csch(x))
    F = a*csch(x)**2/(2*b**2) - csch(x)**3/(3*b) - (a**2 + 2*b**2)*csch(x)/b**3 + log(sinh(x))/a + (a**2 + b**2)**2*log(a + b*csch(x))/(a*b**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_122():
    f = coth(x)**6/(a + b*csch(x))
    F = a*coth(x)**3/(3*b**2) - a*coth(x)/b**2 + a*(a**2 + 3*b**2)*coth(x)/b**4 - coth(x)*csch(x)**3/(4*b) + 3*coth(x)*csch(x)/(8*b) - 3*atanh(cosh(x))/(8*b) - (a**2 + 3*b**2)*coth(x)*csch(x)/(2*b**3) + (a**2 + 3*b**2)*atanh(cosh(x))/(2*b**3) - (a**4 + 3*a**2*b**2 + 3*b**4)*atanh(cosh(x))/b**5 + x/a + 2*(a**2 + b**2)**(sympy.S(5)/2)*atanh((a - b*tanh(x/2))/sqrt(a**2 + b**2))/(a*b**5)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_123():
    f = coth(x)**7/(a + b*csch(x))
    F = a*csch(x)**4/(4*b**2) + a*(a**2 + 3*b**2)*csch(x)**2/(2*b**4) - csch(x)**5/(5*b) - (a**2 + 3*b**2)*csch(x)**3/(3*b**3) - (a**4 + 3*a**2*b**2 + 3*b**4)*csch(x)/b**5 + log(sinh(x))/a + (a**2 + b**2)**3*log(a + b*csch(x))/(a*b**6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_124():
    f = (csch(a*c + b*c*x)**2)**(sympy.S(7)/2)*exp(c*(a + b*x))
    F = 64*sqrt(csch(a*c + b*c*x)**2)*sinh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**3) - 48*sqrt(csch(a*c + b*c*x)**2)*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**4) + 192*sqrt(csch(a*c + b*c*x)**2)*sinh(a*c + b*c*x)/(5*b*c*(1 - exp(2*c*(a + b*x)))**5) - 32*sqrt(csch(a*c + b*c*x)**2)*sinh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**6)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_125():
    f = (csch(a*c + b*c*x)**2)**(sympy.S(5)/2)*exp(c*(a + b*x))
    F = -8*sqrt(csch(a*c + b*c*x)**2)*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**2) + 32*sqrt(csch(a*c + b*c*x)**2)*sinh(a*c + b*c*x)/(3*b*c*(1 - exp(2*c*(a + b*x)))**3) - 4*sqrt(csch(a*c + b*c*x)**2)*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**4)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_126():
    f = (csch(a*c + b*c*x)**2)**(sympy.S(3)/2)*exp(c*(a + b*x))
    F = -2*sqrt(csch(a*c + b*c*x)**2)*exp(4*c*(a + b*x))*sinh(a*c + b*c*x)/(b*c*(1 - exp(2*c*(a + b*x)))**2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_127():
    f = sqrt(csch(a*c + b*c*x)**2)*exp(c*(a + b*x))
    F = sqrt(csch(a*c + b*c*x)**2)*log(1 - exp(2*c*(a + b*x)))*sinh(a*c + b*c*x)/(b*c)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_128():
    f = exp(c*(a + b*x))/sqrt(csch(a*c + b*c*x)**2)
    F = -x*csch(a*c + b*c*x)/(2*sqrt(csch(a*c + b*c*x)**2)) + exp(2*c*(a + b*x))*csch(a*c + b*c*x)/(4*b*c*sqrt(csch(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_129():
    f = exp(c*(a + b*x))/(csch(a*c + b*c*x)**2)**(sympy.S(3)/2)
    F = 3*x*csch(a*c + b*c*x)/(8*sqrt(csch(a*c + b*c*x)**2)) + exp(4*c*(a + b*x))*csch(a*c + b*c*x)/(32*b*c*sqrt(csch(a*c + b*c*x)**2)) - 3*exp(2*c*(a + b*x))*csch(a*c + b*c*x)/(16*b*c*sqrt(csch(a*c + b*c*x)**2)) + exp(-2*c*(a + b*x))*csch(a*c + b*c*x)/(16*b*c*sqrt(csch(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_130():
    f = exp(c*(a + b*x))/(csch(a*c + b*c*x)**2)**(sympy.S(5)/2)
    F = -5*x*csch(a*c + b*c*x)/(16*sqrt(csch(a*c + b*c*x)**2)) + exp(6*c*(a + b*x))*csch(a*c + b*c*x)/(192*b*c*sqrt(csch(a*c + b*c*x)**2)) - 5*exp(4*c*(a + b*x))*csch(a*c + b*c*x)/(128*b*c*sqrt(csch(a*c + b*c*x)**2)) + 5*exp(2*c*(a + b*x))*csch(a*c + b*c*x)/(32*b*c*sqrt(csch(a*c + b*c*x)**2)) - 5*exp(-2*c*(a + b*x))*csch(a*c + b*c*x)/(64*b*c*sqrt(csch(a*c + b*c*x)**2)) + exp(-4*c*(a + b*x))*csch(a*c + b*c*x)/(128*b*c*sqrt(csch(a*c + b*c*x)**2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_131():
    f = x**5/sqrt(csch(2*log(c*x)))
    F = x**6/(7*sqrt(csch(2*log(c*x)))) - 2*x**2/(21*c**4*sqrt(csch(2*log(c*x)))) + 2*elliptic_f(acsc(c*x), -1)/(21*c**7*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_132():
    f = x**4/sqrt(csch(2*log(c*x)))
    F = x**5*(c**4 - 1/x**4)/(6*c**4*sqrt(csch(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_133():
    f = x**3/sqrt(csch(2*log(c*x)))
    F = x**4/(5*sqrt(csch(2*log(c*x)))) - 2/(5*c**4*sqrt(csch(2*log(c*x)))) - 2*elliptic_e(acsc(c*x), -1)/(5*c**5*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x)))) + 2*elliptic_f(acsc(c*x), -1)/(5*c**5*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_134():
    f = x**2/sqrt(csch(2*log(c*x)))
    F = x**3/(4*sqrt(csch(2*log(c*x)))) - atanh(sqrt(1 - 1/(c**4*x**4)))/(4*c**4*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_135():
    f = x/sqrt(csch(2*log(c*x)))
    F = x**2/(3*sqrt(csch(2*log(c*x)))) + 2*elliptic_f(acsc(c*x), -1)/(3*c**3*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_136():
    f = 1/sqrt(csch(2*log(c*x)))
    F = x/(2*sqrt(csch(2*log(c*x)))) + acsc(c**2*x**2)/(2*c**2*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_137():
    f = sqrt(csch(2*log(c*x)))/x
    F = -I*sqrt(I*sinh(2*log(c*x)))*sqrt(csch(2*log(c*x)))*elliptic_f(I*log(c*x) - pi/4, 2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_138():
    f = sqrt(csch(2*log(c*x)))/x**2
    F = -c**2*x*sqrt(1 - 1/(c**4*x**4))*acsc(c**2*x**2)*sqrt(csch(2*log(c*x)))/2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_139():
    f = sqrt(csch(2*log(c*x)))/x**3
    F = -c**3*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x)))*elliptic_e(acsc(c*x), -1) + c**3*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x)))*elliptic_f(acsc(c*x), -1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_140():
    f = sqrt(csch(2*log(c*x)))/x**4
    F = x*(c**4/2 - 1/(2*x**4))*sqrt(csch(2*log(c*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_141():
    f = sqrt(csch(2*log(c*x)))/x**5
    F = -c**5*x*sqrt(1 - 1/(c**4*x**4))*sqrt(csch(2*log(c*x)))*elliptic_f(acsc(c*x), -1)/3 + (c**4/3 - 1/(3*x**4))*sqrt(csch(2*log(c*x)))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_142():
    f = x**8/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**9/(12*csch(2*log(c*x))**(sympy.S(3)/2)) - x**5/((16*c**4 - 16/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) + x/(32*c**4*(c**4 - 1/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) + atanh(sqrt(1 - 1/(c**4*x**4)))/(32*c**12*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_143():
    f = x**7/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**8/(11*csch(2*log(c*x))**(sympy.S(3)/2)) - 6*x**4/((77*c**4 - 77/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) + 4/(77*c**4*(c**4 - 1/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) - 4*elliptic_f(acsc(c*x), -1)/(77*c**11*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_144():
    f = x**6/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**7*(c**4 - 1/x**4)/(10*c**4*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_145():
    f = x**5/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**6/(9*csch(2*log(c*x))**(sympy.S(3)/2)) - 2*x**2/((15*c**4 - 15/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) + 4/(15*c**4*x**2*(c**4 - 1/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) + 4*elliptic_e(acsc(c*x), -1)/(15*c**9*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2)) - 4*elliptic_f(acsc(c*x), -1)/(15*c**9*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_146():
    f = x**4/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**5/(8*csch(2*log(c*x))**(sympy.S(3)/2)) - 3*x/((16*c**4 - 16/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) + 3*atanh(sqrt(1 - 1/(c**4*x**4)))/(16*c**8*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_147():
    f = x**3/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**4/(7*csch(2*log(c*x))**(sympy.S(3)/2)) - 2/((7*c**4 - 7/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) - 4*elliptic_f(acsc(c*x), -1)/(7*c**7*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_148():
    f = x**2/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**3/(6*csch(2*log(c*x))**(sympy.S(3)/2)) - 1/(x*(2*c**4 - 2/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) - acsc(c**2*x**2)/(2*c**6*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_149():
    f = x/csch(2*log(c*x))**(sympy.S(3)/2)
    F = x**2/(5*csch(2*log(c*x))**(sympy.S(3)/2)) - 6/(x**2*(5*c**4 - 5/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) - 12*elliptic_e(acsc(c*x), -1)/(5*c**5*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2)) + 12*elliptic_f(acsc(c*x), -1)/(5*c**5*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_150():
    f = csch(2*log(c*x))**(sympy.S(-3)/2)
    F = x/(4*csch(2*log(c*x))**(sympy.S(3)/2)) + 3/(x**3*(4*c**4 - 4/x**4)*csch(2*log(c*x))**(sympy.S(3)/2)) - 3*atanh(sqrt(1 - 1/(c**4*x**4)))/(4*c**4*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_151():
    f = csch(2*log(c*x))**(sympy.S(3)/2)/x
    F = -cosh(2*log(c*x))*sqrt(csch(2*log(c*x))) - I*elliptic_e(I*log(c*x) - pi/4, 2)/(sqrt(I*sinh(2*log(c*x)))*sqrt(csch(2*log(c*x))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_152():
    f = csch(2*log(c*x))**(sympy.S(3)/2)/x**2
    F = x**3*(-c**4/2 + 1/(2*x**4))*csch(2*log(c*x))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_153():
    f = csch(2*log(c*x))**(sympy.S(3)/2)/x**3
    F = c**5*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*csch(2*log(c*x))**(sympy.S(3)/2)*elliptic_f(acsc(c*x), -1)/2 + x**2*(-c**4/2 + 1/(2*x**4))*csch(2*log(c*x))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_154():
    f = csch(2*log(c*x))**(sympy.S(3)/2)/x**4
    F = c**6*x**3*(1 - 1/(c**4*x**4))**(sympy.S(3)/2)*acsc(c**2*x**2)*csch(2*log(c*x))**(sympy.S(3)/2)/2 + x*(-c**4/2 + 1/(2*x**4))*csch(2*log(c*x))**(sympy.S(3)/2)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_155():
    f = csch(a + b*log(c*x**n))
    F = -2*x*(c*x**n)**b*exp(a)*hyper((1, (b + 1/n)/(2*b)), (sympy.S(3)/2 + 1/(2*b*n),), (c*x**n)**(2*b)*exp(2*a))/(b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_156():
    f = csch(a + b*log(c*x**n))**2
    F = 4*x*(c*x**n)**(2*b)*exp(2*a)*hyper((2, 1 + 1/(2*b*n)), (2 + 1/(2*b*n),), (c*x**n)**(2*b)*exp(2*a))/(2*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_157():
    f = csch(a + b*log(c*x**n))**3
    F = -8*x*(c*x**n)**(3*b)*exp(3*a)*hyper((3, (3*b + 1/n)/(2*b)), (sympy.S(5)/2 + 1/(2*b*n),), (c*x**n)**(2*b)*exp(2*a))/(3*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_158():
    f = csch(a + b*log(c*x**n))**4
    F = 16*x*(c*x**n)**(4*b)*exp(4*a)*hyper((4, 2 + 1/(2*b*n)), (3 + 1/(2*b*n),), (c*x**n)**(2*b)*exp(2*a))/(4*b*n + 1)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_159():
    f = 2*b**2*n**2*csch(a + b*log(c*x**n))**3 - (-b**2*n**2 + 1)*csch(a + b*log(c*x**n))
    F = -b*n*x*coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n)) - x*csch(a + b*log(c*x**n))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_160():
    f = csch(a + 2*log(c*sqrt(x)))**3
    F = -2*c**6*exp(-a)/(c**4 - exp(-2*a)/x**2)**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_161():
    f = csch(a + 2*log(c/sqrt(x)))**3
    F = 2*c**2*exp(-3*a)/(-c**4/x**2 + exp(-2*a))**2
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_162():
    f = csch(a + log(c*x**n)/(n*(p - 2)))**p
    F = -x*(2 - p)*(-(c*x**n)**(2/(n*(2 - p)))*exp(-2*a) + 1)*exp(2*a)*csch(a - log(c*x**n)/(n*(2 - p)))**p/((c*x**n)**(2/(n*(2 - p)))*(2 - 2*p))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_163():
    f = csch(a - log(c*x**n)/(n*(p - 2)))**p
    F = x*(1 - exp(-2*a)/(c*x**n)**(2/(n*(2 - p))))*(2 - p)*csch(a + log(c*x**n)/(n*(2 - p)))**p/(2 - 2*p)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_164():
    f = csch(a + b*log(c*x**n))/x
    F = -atanh(cosh(a + b*log(c*x**n)))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_165():
    f = csch(a + b*log(c*x**n))**2/x
    F = -coth(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_166():
    f = csch(a + b*log(c*x**n))**3/x
    F = -coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n))/(2*b*n) + atanh(cosh(a + b*log(c*x**n)))/(2*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_167():
    f = csch(a + b*log(c*x**n))**4/x
    F = -coth(a + b*log(c*x**n))**3/(3*b*n) + coth(a + b*log(c*x**n))/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_168():
    f = csch(a + b*log(c*x**n))**5/x
    F = -coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**3/(4*b*n) + 3*coth(a + b*log(c*x**n))*csch(a + b*log(c*x**n))/(8*b*n) - 3*atanh(cosh(a + b*log(c*x**n)))/(8*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_169():
    f = csch(a + b*log(c*x**n))**(sympy.S(5)/2)/x
    F = 2*I*sqrt(I*sinh(a + b*log(c*x**n)))*sqrt(csch(a + b*log(c*x**n)))*elliptic_f(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(3*b*n) - 2*cosh(a + b*log(c*x**n))*csch(a + b*log(c*x**n))**(sympy.S(3)/2)/(3*b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_170():
    f = csch(a + b*log(c*x**n))**(sympy.S(3)/2)/x
    F = -2*cosh(a + b*log(c*x**n))*sqrt(csch(a + b*log(c*x**n)))/(b*n) - 2*I*elliptic_e(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(b*n*sqrt(I*sinh(a + b*log(c*x**n)))*sqrt(csch(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_171():
    f = sqrt(csch(a + b*log(c*x**n)))/x
    F = -2*I*sqrt(I*sinh(a + b*log(c*x**n)))*sqrt(csch(a + b*log(c*x**n)))*elliptic_f(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(b*n)
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_172():
    f = 1/(x*sqrt(csch(a + b*log(c*x**n))))
    F = -2*I*elliptic_e(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(b*n*sqrt(I*sinh(a + b*log(c*x**n)))*sqrt(csch(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_173():
    f = 1/(x*csch(a + b*log(c*x**n))**(sympy.S(3)/2))
    F = 2*I*sqrt(I*sinh(a + b*log(c*x**n)))*sqrt(csch(a + b*log(c*x**n)))*elliptic_f(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(3*b*n) + 2*cosh(a + b*log(c*x**n))/(3*b*n*sqrt(csch(a + b*log(c*x**n))))
    assert integrate(f, x) == F


def test_integrate_6_Hyperbolic_functions_6_6_Hyperbolic_cosecant_6_6_3_Hyperbolic_cosecant_functions_174():
    f = 1/(x*csch(a + b*log(c*x**n))**(sympy.S(5)/2))
    F = 2*cosh(a + b*log(c*x**n))/(5*b*n*csch(a + b*log(c*x**n))**(sympy.S(3)/2)) + 6*I*elliptic_e(I*a/2 + I*b*log(c*x**n)/2 - pi/4, 2)/(5*b*n*sqrt(I*sinh(a + b*log(c*x**n)))*sqrt(csch(a + b*log(c*x**n))))
    assert integrate(f, x) == F

