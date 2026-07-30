"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.1 Sine/4.1.1.1 (a+b sin)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, n = symbols('a b c d n')

def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_1():
    f = (a*sin(c + d*x) + a)**(sympy.S(7)/2)
    F = -256*a**4*cos(c + d*x)/(35*d*sqrt(a*sin(c + d*x) + a)) - 64*a**3*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(35*d) - 24*a**2*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(35*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_2():
    f = (a*sin(c + d*x) + a)**(sympy.S(5)/2)
    F = -64*a**3*cos(c + d*x)/(15*d*sqrt(a*sin(c + d*x) + a)) - 16*a**2*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(15*d) - 2*a*(a*sin(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_3():
    f = (a*sin(c + d*x) + a)**(sympy.S(3)/2)
    F = -8*a**2*cos(c + d*x)/(3*d*sqrt(a*sin(c + d*x) + a)) - 2*a*sqrt(a*sin(c + d*x) + a)*cos(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_4():
    f = sqrt(a*sin(c + d*x) + a)
    F = -2*a*cos(c + d*x)/(d*sqrt(a*sin(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_5():
    f = 1/sqrt(a*sin(c + d*x) + a)
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_6():
    f = (a*sin(c + d*x) + a)**(sympy.S(-3)/2)
    F = -cos(c + d*x)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_7():
    f = (a*sin(c + d*x) + a)**(sympy.S(-5)/2)
    F = -cos(c + d*x)/(4*d*(a*sin(c + d*x) + a)**(sympy.S(5)/2)) - 3*cos(c + d*x)/(16*a*d*(a*sin(c + d*x) + a)**(sympy.S(3)/2)) - 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*cos(c + d*x)/(2*sqrt(a*sin(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_8():
    f = (a*sin(c + d*x) + a)**(sympy.S(4)/3)
    F = -2*2**(sympy.S(5)/6)*a*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(sin(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_9():
    f = (a*sin(c + d*x) + a)**(sympy.S(2)/3)
    F = -2*2**(sympy.S(1)/6)*(a*sin(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(sin(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_10():
    f = (a*sin(c + d*x) + a)**(sympy.S(1)/3)
    F = -2**(sympy.S(5)/6)*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(sin(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_11():
    f = (a*sin(c + d*x) + a)**(sympy.S(-1)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_12():
    f = (a*sin(c + d*x) + a)**(sympy.S(-2)/3)
    F = -2**(sympy.S(5)/6)*(sin(c + d*x) + 1)**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(2*d*(a*sin(c + d*x) + a)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_13():
    f = (a*sin(c + d*x) + a)**(sympy.S(-4)/3)
    F = -2**(sympy.S(1)/6)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S(11)/6), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(2*a*d*(a*sin(c + d*x) + a)**(sympy.S(1)/3)*(sin(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_14():
    f = (a*sin(c + d*x) + a)**n
    F = -2**(n + sympy.S.Half)*(a*sin(c + d*x) + a)**n*(sin(c + d*x) + 1)**(-n + sympy.S(-1)/2)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_15():
    f = (-a*sin(c + d*x) + a)**n
    F = 2**(n + sympy.S.Half)*(1 - sin(c + d*x))**(-n + sympy.S(-1)/2)*(-a*sin(c + d*x) + a)**n*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sin(c + d*x)/2 + sympy.S.Half)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_16():
    f = (2*sin(c + d*x) + 2)**n
    F = -2**(2*n + sympy.S.Half)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sympy.S.Half - sin(c + d*x)/2)/(d*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_17():
    f = (2 - 2*sin(c + d*x))**n
    F = 2**(2*n + sympy.S.Half)*cos(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), sin(c + d*x)/2 + sympy.S.Half)/(d*sqrt(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_18():
    f = 1/(3*sin(c + d*x) + 5)
    F = x/4 + atan(cos(c + d*x)/(sin(c + d*x) + 3))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_19():
    f = (3*sin(c + d*x) + 5)**(-2)
    F = 5*x/64 + 5*atan(cos(c + d*x)/(sin(c + d*x) + 3))/(32*d) + 3*cos(c + d*x)/(16*d*(3*sin(c + d*x) + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_20():
    f = (3*sin(c + d*x) + 5)**(-3)
    F = 59*x/2048 + 59*atan(cos(c + d*x)/(sin(c + d*x) + 3))/(1024*d) + 45*cos(c + d*x)/(512*d*(3*sin(c + d*x) + 5)) + 3*cos(c + d*x)/(32*d*(3*sin(c + d*x) + 5)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_21():
    f = (3*sin(c + d*x) + 5)**(-4)
    F = 385*x/32768 + 385*atan(cos(c + d*x)/(sin(c + d*x) + 3))/(16384*d) + 311*cos(c + d*x)/(8192*d*(3*sin(c + d*x) + 5)) + 25*cos(c + d*x)/(512*d*(3*sin(c + d*x) + 5)**2) + cos(c + d*x)/(16*d*(3*sin(c + d*x) + 5)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_22():
    f = 1/(5 - 3*sin(c + d*x))
    F = x/4 - atan(cos(c + d*x)/(3 - sin(c + d*x)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_23():
    f = (5 - 3*sin(c + d*x))**(-2)
    F = 5*x/64 - 5*atan(cos(c + d*x)/(3 - sin(c + d*x)))/(32*d) - 3*cos(c + d*x)/(16*d*(5 - 3*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_24():
    f = (5 - 3*sin(c + d*x))**(-3)
    F = 59*x/2048 - 59*atan(cos(c + d*x)/(3 - sin(c + d*x)))/(1024*d) - 45*cos(c + d*x)/(512*d*(5 - 3*sin(c + d*x))) - 3*cos(c + d*x)/(32*d*(5 - 3*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_25():
    f = (5 - 3*sin(c + d*x))**(-4)
    F = 385*x/32768 - 385*atan(cos(c + d*x)/(3 - sin(c + d*x)))/(16384*d) - 311*cos(c + d*x)/(8192*d*(5 - 3*sin(c + d*x))) - 25*cos(c + d*x)/(512*d*(5 - 3*sin(c + d*x))**2) - cos(c + d*x)/(16*d*(5 - 3*sin(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_26():
    f = 1/(3*sin(c + d*x) - 5)
    F = -x/4 + atan(cos(c + d*x)/(3 - sin(c + d*x)))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_27():
    f = (3*sin(c + d*x) - 5)**(-2)
    F = 5*x/64 - 5*atan(cos(c + d*x)/(3 - sin(c + d*x)))/(32*d) - 3*cos(c + d*x)/(16*d*(5 - 3*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_28():
    f = (3*sin(c + d*x) - 5)**(-3)
    F = -59*x/2048 + 59*atan(cos(c + d*x)/(3 - sin(c + d*x)))/(1024*d) + 45*cos(c + d*x)/(512*d*(5 - 3*sin(c + d*x))) + 3*cos(c + d*x)/(32*d*(5 - 3*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_29():
    f = (3*sin(c + d*x) - 5)**(-4)
    F = 385*x/32768 - 385*atan(cos(c + d*x)/(3 - sin(c + d*x)))/(16384*d) - 311*cos(c + d*x)/(8192*d*(5 - 3*sin(c + d*x))) - 25*cos(c + d*x)/(512*d*(5 - 3*sin(c + d*x))**2) - cos(c + d*x)/(16*d*(5 - 3*sin(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_30():
    f = 1/(-3*sin(c + d*x) - 5)
    F = -x/4 - atan(cos(c + d*x)/(sin(c + d*x) + 3))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_31():
    f = (-3*sin(c + d*x) - 5)**(-2)
    F = 5*x/64 + 5*atan(cos(c + d*x)/(sin(c + d*x) + 3))/(32*d) + 3*cos(c + d*x)/(16*d*(3*sin(c + d*x) + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_32():
    f = (-3*sin(c + d*x) - 5)**(-3)
    F = -59*x/2048 - 59*atan(cos(c + d*x)/(sin(c + d*x) + 3))/(1024*d) - 45*cos(c + d*x)/(512*d*(3*sin(c + d*x) + 5)) - 3*cos(c + d*x)/(32*d*(3*sin(c + d*x) + 5)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_33():
    f = (-3*sin(c + d*x) - 5)**(-4)
    F = 385*x/32768 + 385*atan(cos(c + d*x)/(sin(c + d*x) + 3))/(16384*d) + 311*cos(c + d*x)/(8192*d*(3*sin(c + d*x) + 5)) + 25*cos(c + d*x)/(512*d*(3*sin(c + d*x) + 5)**2) + cos(c + d*x)/(16*d*(3*sin(c + d*x) + 5)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_34():
    f = 1/(5*sin(c + d*x) + 3)
    F = -log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(4*d) + log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_35():
    f = (5*sin(c + d*x) + 3)**(-2)
    F = 3*log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(64*d) - 3*log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(64*d) - 5*cos(c + d*x)/(16*d*(5*sin(c + d*x) + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_36():
    f = (5*sin(c + d*x) + 3)**(-3)
    F = -43*log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(2048*d) + 43*log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(2048*d) + 45*cos(c + d*x)/(512*d*(5*sin(c + d*x) + 3)) - 5*cos(c + d*x)/(32*d*(5*sin(c + d*x) + 3)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_37():
    f = (5*sin(c + d*x) + 3)**(-4)
    F = 279*log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(32768*d) - 279*log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(32768*d) - 995*cos(c + d*x)/(24576*d*(5*sin(c + d*x) + 3)) + 25*cos(c + d*x)/(512*d*(5*sin(c + d*x) + 3)**2) - 5*cos(c + d*x)/(48*d*(5*sin(c + d*x) + 3)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_38():
    f = 1/(3 - 5*sin(c + d*x))
    F = -log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(4*d) + log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_39():
    f = (3 - 5*sin(c + d*x))**(-2)
    F = 3*log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(64*d) - 3*log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(64*d) + 5*cos(c + d*x)/(16*d*(3 - 5*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_40():
    f = (3 - 5*sin(c + d*x))**(-3)
    F = -43*log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(2048*d) + 43*log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(2048*d) - 45*cos(c + d*x)/(512*d*(3 - 5*sin(c + d*x))) + 5*cos(c + d*x)/(32*d*(3 - 5*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_41():
    f = (3 - 5*sin(c + d*x))**(-4)
    F = 279*log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(32768*d) - 279*log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(32768*d) + 995*cos(c + d*x)/(24576*d*(3 - 5*sin(c + d*x))) - 25*cos(c + d*x)/(512*d*(3 - 5*sin(c + d*x))**2) + 5*cos(c + d*x)/(48*d*(3 - 5*sin(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_42():
    f = 1/(5*sin(c + d*x) - 3)
    F = log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(4*d) - log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_43():
    f = (5*sin(c + d*x) - 3)**(-2)
    F = 3*log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(64*d) - 3*log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(64*d) + 5*cos(c + d*x)/(16*d*(3 - 5*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_44():
    f = (5*sin(c + d*x) - 3)**(-3)
    F = 43*log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(2048*d) - 43*log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(2048*d) + 45*cos(c + d*x)/(512*d*(3 - 5*sin(c + d*x))) - 5*cos(c + d*x)/(32*d*(3 - 5*sin(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_45():
    f = (5*sin(c + d*x) - 3)**(-4)
    F = 279*log(-3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(32768*d) - 279*log(-sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(32768*d) + 995*cos(c + d*x)/(24576*d*(3 - 5*sin(c + d*x))) - 25*cos(c + d*x)/(512*d*(3 - 5*sin(c + d*x))**2) + 5*cos(c + d*x)/(48*d*(3 - 5*sin(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_46():
    f = 1/(-5*sin(c + d*x) - 3)
    F = log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(4*d) - log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_47():
    f = (-5*sin(c + d*x) - 3)**(-2)
    F = 3*log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(64*d) - 3*log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(64*d) - 5*cos(c + d*x)/(16*d*(5*sin(c + d*x) + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_48():
    f = (-5*sin(c + d*x) - 3)**(-3)
    F = 43*log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(2048*d) - 43*log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(2048*d) - 45*cos(c + d*x)/(512*d*(5*sin(c + d*x) + 3)) + 5*cos(c + d*x)/(32*d*(5*sin(c + d*x) + 3)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_49():
    f = (-5*sin(c + d*x) - 3)**(-4)
    F = 279*log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(32768*d) - 279*log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(32768*d) - 995*cos(c + d*x)/(24576*d*(5*sin(c + d*x) + 3)) + 25*cos(c + d*x)/(512*d*(5*sin(c + d*x) + 3)**2) - 5*cos(c + d*x)/(48*d*(5*sin(c + d*x) + 3)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_50():
    f = (a + b*sin(c + d*x))**(sympy.S(7)/2)
    F = -24*a*b*(a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(35*d) + 32*a*sqrt(a + b*sin(c + d*x))*(11*a**2 + 13*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(105*d*sqrt((a + b*sin(c + d*x))/(a + b))) - 2*b*(a + b*sin(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)/(7*d) - 2*b*sqrt(a + b*sin(c + d*x))*(71*a**2 + 25*b**2)*cos(c + d*x)/(105*d) - sqrt((a + b*sin(c + d*x))/(a + b))*(142*a**4 - 92*a**2*b**2 - 50*b**4)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(105*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_51():
    f = (a + b*sin(c + d*x))**(sympy.S(5)/2)
    F = -16*a*b*sqrt(a + b*sin(c + d*x))*cos(c + d*x)/(15*d) - 16*a*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15*d*sqrt(a + b*sin(c + d*x))) - 2*b*(a + b*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*d) + sqrt(a + b*sin(c + d*x))*(46*a**2 + 18*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15*d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_52():
    f = (a + b*sin(c + d*x))**(sympy.S(3)/2)
    F = 8*a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt((a + b*sin(c + d*x))/(a + b))) - 2*b*sqrt(a + b*sin(c + d*x))*cos(c + d*x)/(3*d) - sqrt((a + b*sin(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_53():
    f = sqrt(a + b*sin(c + d*x))
    F = 2*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_54():
    f = 1/sqrt(a + b*sin(c + d*x))
    F = 2*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_55():
    f = (a + b*sin(c + d*x))**(sympy.S(-3)/2)
    F = 2*b*cos(c + d*x)/(d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)) + 2*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_56():
    f = (a + b*sin(c + d*x))**(sympy.S(-5)/2)
    F = 8*a*b*cos(c + d*x)/(3*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) + 8*a*sqrt(a + b*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(3*d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)**2) + 2*b*cos(c + d*x)/(d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) - 2*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(d*sqrt(a + b*sin(c + d*x))*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_57():
    f = (a + b*sin(c + d*x))**(sympy.S(-7)/2)
    F = 16*a*b*cos(c + d*x)/(15*d*(a + b*sin(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)**2) - 16*a*sqrt((a + b*sin(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**2) + 2*b*(23*a**2 + 9*b**2)*cos(c + d*x)/(15*d*sqrt(a + b*sin(c + d*x))*(a**2 - b**2)**3) + 2*b*cos(c + d*x)/(d*(a + b*sin(c + d*x))**(sympy.S(5)/2)*(5*a**2 - 5*b**2)) + sqrt(a + b*sin(c + d*x))*(46*a**2 + 18*b**2)*elliptic_e(c/2 + d*x/2 - pi/4, 2*b/(a + b))/(15*d*sqrt((a + b*sin(c + d*x))/(a + b))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_58():
    f = (a + b*sin(c + d*x))**(sympy.S(4)/3)
    F = -sqrt(2)*(a + b)*(a + b*sin(c + d*x))**(sympy.S(1)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sin(c + d*x))/(a + b), sympy.S.Half - sin(c + d*x)/2)/(d*((a + b*sin(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_59():
    f = (a + b*sin(c + d*x))**(sympy.S(2)/3)
    F = -sqrt(2)*(a + b*sin(c + d*x))**(sympy.S(2)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sin(c + d*x))/(a + b), sympy.S.Half - sin(c + d*x)/2)/(d*((a + b*sin(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_60():
    f = (a + b*sin(c + d*x))**(sympy.S(1)/3)
    F = -sqrt(2)*(a + b*sin(c + d*x))**(sympy.S(1)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sin(c + d*x))/(a + b), sympy.S.Half - sin(c + d*x)/2)/(d*((a + b*sin(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_61():
    f = (a + b*sin(c + d*x))**(sympy.S(-1)/3)
    F = -sqrt(2)*((a + b*sin(c + d*x))/(a + b))**(sympy.S(1)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sin(c + d*x))/(a + b), sympy.S.Half - sin(c + d*x)/2)/(d*(a + b*sin(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_62():
    f = (a + b*sin(c + d*x))**(sympy.S(-2)/3)
    F = -sqrt(2)*((a + b*sin(c + d*x))/(a + b))**(sympy.S(2)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(d*(a + b*sin(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_63():
    f = (a + b*sin(c + d*x))**(sympy.S(-4)/3)
    F = -sqrt(2)*((a + b*sin(c + d*x))/(a + b))**(sympy.S(1)/3)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(4)/3, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(d*(a + b)*(a + b*sin(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_64():
    f = (a + b*sin(c + d*x))**n
    F = -sqrt(2)*(a + b*sin(c + d*x))**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, b*(1 - sin(c + d*x))/(a + b))/(d*((a + b*sin(c + d*x))/(a + b))**n*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_65():
    f = (4*sin(c + d*x) + 3)**n
    F = -sqrt(2)*7**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, sympy.S(4)/7 - 4*sin(c + d*x)/7)/(d*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_66():
    f = (3 - 4*sin(c + d*x))**n
    F = sqrt(2)*7**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sin(c + d*x)/2 + sympy.S.Half, 4*sin(c + d*x)/7 + sympy.S(4)/7)/(d*sqrt(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_67():
    f = (3*sin(c + d*x) + 4)**n
    F = sqrt(2)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sin(c + d*x)/2 + sympy.S.Half, -3*sin(c + d*x) - 3)/(d*sqrt(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_68():
    f = (4 - 3*sin(c + d*x))**n
    F = sqrt(2)*7**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sin(c + d*x)/2 + sympy.S.Half, 3*sin(c + d*x)/7 + sympy.S(3)/7)/(d*sqrt(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_69():
    f = (4*sin(c + d*x) - 3)**n
    F = -sqrt(2)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sympy.S.Half - sin(c + d*x)/2, 4 - 4*sin(c + d*x))/(d*sqrt(sin(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_70():
    f = (-4*sin(c + d*x) - 3)**n
    F = sqrt(2)*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sin(c + d*x)/2 + sympy.S.Half, 4*sin(c + d*x) + 4)/(d*sqrt(1 - sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_71():
    f = (3*sin(c + d*x) - 4)**n
    F = sqrt(2)*7**n*(3*sin(c + d*x) - 4)**n*cos(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, -n, sympy.S(3)/2, sin(c + d*x)/2 + sympy.S.Half, 3*sin(c + d*x)/7 + sympy.S(3)/7)/(d*sqrt(1 - sin(c + d*x))*(4 - 3*sin(c + d*x))**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_1_Sine_4_1_1_1_a_plus_b_sin_pow_n_72():
    f = (-3*sin(c + d*x) - 4)**n
    F = -sqrt(7)*(-3*sin(c + d*x) - 4)**(n + 1)*sqrt(-sin(c + d*x) - 1)*cos(c + d*x)*appellf1(n + 1, sympy.S.Half, sympy.S.Half, n + 2, 3*sin(c + d*x)/7 + sympy.S(4)/7, 3*sin(c + d*x) + 4)/(7*d*sqrt(1 - sin(c + d*x))*(n + 1)*(sin(c + d*x) + 1))
    assert integrate(f, x) == F

