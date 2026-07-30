"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.6 Cosecant/4.6.1.2 (d csc)^n (a+b csc)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_1():
    f = csc(x)**5/(a*csc(x) + a)
    F = cot(x)*csc(x)**3/(a*csc(x) + a) - 4*cot(x)**3/(3*a) + 3*cot(x)*csc(x)/(2*a) - 4*cot(x)/a + 3*atanh(cos(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_2():
    f = csc(x)**4/(a*csc(x) + a)
    F = cot(x)*csc(x)**2/(a*csc(x) + a) - 3*cot(x)*csc(x)/(2*a) + 2*cot(x)/a - 3*atanh(cos(x))/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_3():
    f = csc(x)**3/(a*csc(x) + a)
    F = -cot(x)/(a*csc(x) + a) - cot(x)/a + atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_4():
    f = csc(x)**2/(a*csc(x) + a)
    F = cot(x)/(a*csc(x) + a) - atanh(cos(x))/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_5():
    f = csc(x)/(a*csc(x) + a)
    F = -cot(x)/(a*csc(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_6():
    f = 1/(a*csc(c + d*x) + a)
    F = cot(c + d*x)/(d*(a*csc(c + d*x) + a)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_7():
    f = sin(x)/(a*csc(x) + a)
    F = cos(x)/(a*csc(x) + a) - x/a - 2*cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_8():
    f = sin(x)**2/(a*csc(x) + a)
    F = sin(x)*cos(x)/(a*csc(x) + a) + 3*x/(2*a) - 3*sin(x)*cos(x)/(2*a) + 2*cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_9():
    f = sin(x)**3/(a*csc(x) + a)
    F = sin(x)**2*cos(x)/(a*csc(x) + a) - 3*x/(2*a) + 3*sin(x)*cos(x)/(2*a) + 4*cos(x)**3/(3*a) - 4*cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_10():
    f = sin(x)**4/(a*csc(x) + a)
    F = sin(x)**3*cos(x)/(a*csc(x) + a) + 15*x/(8*a) - 5*sin(x)**3*cos(x)/(4*a) - 15*sin(x)*cos(x)/(8*a) - 4*cos(x)**3/(3*a) + 4*cos(x)/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_11():
    f = (a*csc(c + d*x) + a)**(-2)
    F = cot(c + d*x)/(3*d*(a*csc(c + d*x) + a)**2) + x/a**2 + 4*cot(c + d*x)/(3*a**2*d*(csc(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_12():
    f = (a*csc(c + d*x) + a)**(-3)
    F = 22*cot(c + d*x)/(15*d*(a**3*csc(c + d*x) + a**3)) + cot(c + d*x)/(5*d*(a*csc(c + d*x) + a)**3) + 7*cot(c + d*x)/(15*a*d*(a*csc(c + d*x) + a)**2) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_13():
    f = (a*csc(x) + a)**(sympy.S(5)/2)
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(a)*cot(x)/sqrt(a*csc(x) + a)) - 14*a**3*cot(x)/(3*sqrt(a*csc(x) + a)) - 2*a**2*sqrt(a*csc(x) + a)*cot(x)/3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_14():
    f = (a*csc(x) + a)**(sympy.S(3)/2)
    F = -2*a**(sympy.S(3)/2)*atan(sqrt(a)*cot(x)/sqrt(a*csc(x) + a)) - 2*a**2*cot(x)/sqrt(a*csc(x) + a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_15():
    f = sqrt(a*csc(x) + a)
    F = -2*sqrt(a)*atan(sqrt(a)*cot(x)/sqrt(a*csc(x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_16():
    f = 1/sqrt(a*csc(x) + a)
    F = -2*atan(sqrt(a)*cot(x)/sqrt(a*csc(x) + a))/sqrt(a) + sqrt(2)*atan(sqrt(2)*sqrt(a)*cot(x)/(2*sqrt(a*csc(x) + a)))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_17():
    f = (a*csc(x) + a)**(sympy.S(-3)/2)
    F = cot(x)/(2*(a*csc(x) + a)**(sympy.S(3)/2)) - 2*atan(sqrt(a)*cot(x)/sqrt(a*csc(x) + a))/a**(sympy.S(3)/2) + 5*sqrt(2)*atan(sqrt(2)*sqrt(a)*cot(x)/(2*sqrt(a*csc(x) + a)))/(4*a**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_18():
    f = (a*csc(x) + a)**(sympy.S(-5)/2)
    F = cot(x)/(4*(a*csc(x) + a)**(sympy.S(5)/2)) + 11*cot(x)/(16*a*(a*csc(x) + a)**(sympy.S(3)/2)) - 2*atan(sqrt(a)*cot(x)/sqrt(a*csc(x) + a))/a**(sympy.S(5)/2) + 43*sqrt(2)*atan(sqrt(2)*sqrt(a)*cot(x)/(2*sqrt(a*csc(x) + a)))/(32*a**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_19():
    f = sqrt(a*csc(e + f*x) + a)*sqrt(csc(e + f*x))
    F = -2*sqrt(a)*asinh(sqrt(a)*cot(e + f*x)/sqrt(a*csc(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_20():
    f = sqrt(-csc(e + f*x))*sqrt(-a*csc(e + f*x) + a)
    F = -2*sqrt(a)*asinh(sqrt(a)*cot(e + f*x)/sqrt(-a*csc(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_21():
    f = sqrt(a*csc(c + d*x) + a)*csc(c + d*x)**(sympy.S(4)/3)
    F = -4*3**(sympy.S(3)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*cot(c + d*x)*elliptic_f(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(5*d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) - 6*a*cos(c + d*x)*csc(c + d*x)**(sympy.S(4)/3)/(5*d*sqrt(a*csc(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_22():
    f = sqrt(a*csc(c + d*x) + a)*csc(c + d*x)**(sympy.S(1)/3)
    F = -2*3**(sympy.S(3)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*cot(c + d*x)*elliptic_f(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_23():
    f = sqrt(a*csc(c + d*x) + a)/csc(c + d*x)**(sympy.S(2)/3)
    F = -3**(sympy.S(3)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*sqrt(sqrt(3) + 2)*cot(c + d*x)*elliptic_f(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(2*d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) - 3*a*cos(c + d*x)*csc(c + d*x)**(sympy.S(1)/3)/(2*d*sqrt(a*csc(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_24():
    f = sqrt(a*csc(c + d*x) + a)*csc(c + d*x)**(sympy.S(5)/3)
    F = -12*3**(sympy.S(1)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*cot(c + d*x)*elliptic_e(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(7*d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) + 8*sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*cot(c + d*x)*elliptic_f(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(7*d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) - 6*a*cos(c + d*x)*csc(c + d*x)**(sympy.S(5)/3)/(7*d*sqrt(a*csc(c + d*x) + a)) + 24*a*cot(c + d*x)/(7*d*sqrt(a*csc(c + d*x) + a)*(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_25():
    f = sqrt(a*csc(c + d*x) + a)*csc(c + d*x)**(sympy.S(2)/3)
    F = -3*3**(sympy.S(1)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*cot(c + d*x)*elliptic_e(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) + 2*sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*cot(c + d*x)*elliptic_f(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) + 6*a*cot(c + d*x)/(d*sqrt(a*csc(c + d*x) + a)*(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_26():
    f = sqrt(a*csc(c + d*x) + a)/csc(c + d*x)**(sympy.S(1)/3)
    F = 3*3**(sympy.S(1)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*cot(c + d*x)*elliptic_e(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(2*d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) - sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*cot(c + d*x)*elliptic_f(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) - 3*a*cos(c + d*x)*csc(c + d*x)**(sympy.S(2)/3)/(d*sqrt(a*csc(c + d*x) + a)) - 3*a*cot(c + d*x)/(d*sqrt(a*csc(c + d*x) + a)*(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_27():
    f = sqrt(a*csc(c + d*x) + a)/csc(c + d*x)**(sympy.S(4)/3)
    F = 15*3**(sympy.S(1)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*sqrt(2 - sqrt(3))*cot(c + d*x)*elliptic_e(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(16*d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) - 5*sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt((csc(c + d*x)**(sympy.S(2)/3) + csc(c + d*x)**(sympy.S(1)/3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(1 - csc(c + d*x)**(sympy.S(1)/3))*cot(c + d*x)*elliptic_f(asin((-csc(c + d*x)**(sympy.S(1)/3) - sqrt(3) + 1)/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))), -7 - 4*sqrt(3))/(8*d*sqrt((1 - csc(c + d*x)**(sympy.S(1)/3))/(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3))**2)*(-a*csc(c + d*x) + a)*sqrt(a*csc(c + d*x) + a)) - 15*a*cos(c + d*x)*csc(c + d*x)**(sympy.S(2)/3)/(8*d*sqrt(a*csc(c + d*x) + a)) - 3*a*cos(c + d*x)/(4*d*sqrt(a*csc(c + d*x) + a)*csc(c + d*x)**(sympy.S(1)/3)) - 15*a*cot(c + d*x)/(8*d*sqrt(a*csc(c + d*x) + a)*(-csc(c + d*x)**(sympy.S(1)/3) + 1 + sqrt(3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_28():
    f = sqrt(a*csc(c + d*x) + a)*csc(c + d*x)**n
    F = -2*a*cot(c + d*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - csc(c + d*x))/(d*sqrt(a*csc(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_29():
    f = sqrt(-a*csc(c + d*x) + a)*csc(c + d*x)**n
    F = -2*a*cos(c + d*x)*csc(c + d*x)**(n + 1)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), csc(c + d*x) + 1)/(d*(-csc(c + d*x))**n*sqrt(-a*csc(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_30():
    f = (a*csc(e + f*x) + a)**m*csc(e + f*x)**3
    F = -2**(m + sympy.S.Half)*(a*csc(e + f*x) + a)**m*(csc(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(m**2 + m + 1)*cot(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - csc(e + f*x)/2)/(f*(m + 1)*(m + 2)) + (a*csc(e + f*x) + a)**m*cot(e + f*x)/(f*(m**2 + 3*m + 2)) - (a*csc(e + f*x) + a)**(m + 1)*cot(e + f*x)/(a*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_31():
    f = (a*csc(e + f*x) + a)**m*csc(e + f*x)**2
    F = -2**(m + sympy.S.Half)*m*(a*csc(e + f*x) + a)**m*(csc(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cot(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - csc(e + f*x)/2)/(f*(m + 1)) - (a*csc(e + f*x) + a)**m*cot(e + f*x)/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_32():
    f = (a*csc(e + f*x) + a)**m*csc(e + f*x)
    F = -2**(m + sympy.S.Half)*(a*csc(e + f*x) + a)**m*(csc(e + f*x) + 1)**(-m + sympy.S(-1)/2)*cot(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - csc(e + f*x)/2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_33():
    f = (a*csc(e + f*x) + a)**m
    F = -sqrt(2)*(a*csc(e + f*x) + a)**m*cot(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1, m + sympy.S(3)/2, csc(e + f*x)/2 + sympy.S.Half, csc(e + f*x) + 1)/(f*sqrt(1 - csc(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_34():
    f = (a*csc(e + f*x) + a)**m*sin(e + f*x)
    F = sqrt(2)*(a*csc(e + f*x) + a)**m*cot(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 2, m + sympy.S(3)/2, csc(e + f*x)/2 + sympy.S.Half, csc(e + f*x) + 1)/(f*sqrt(1 - csc(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_35():
    f = (a*csc(e + f*x) + a)**m*sin(e + f*x)**2
    F = -sqrt(2)*(a*csc(e + f*x) + a)**m*cot(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 3, m + sympy.S(3)/2, csc(e + f*x)/2 + sympy.S.Half, csc(e + f*x) + 1)/(f*sqrt(1 - csc(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_36():
    f = (a + b*csc(c + d*x))**4
    F = a**4*x - 4*a*b**3*cot(c + d*x)*csc(c + d*x)/(3*d) - 2*a*b*(2*a**2 + b**2)*atanh(cos(c + d*x))/d - b**2*(a + b*csc(c + d*x))**2*cot(c + d*x)/(3*d) - b**2*(17*a**2 + 2*b**2)*cot(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_37():
    f = (a + b*csc(c + d*x))**3
    F = a**3*x - 5*a*b**2*cot(c + d*x)/(2*d) - b**2*(a + b*csc(c + d*x))*cot(c + d*x)/(2*d) - b*(6*a**2 + b**2)*atanh(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_38():
    f = (a + b*csc(c + d*x))**2
    F = a**2*x - 2*a*b*atanh(cos(c + d*x))/d - b**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_39():
    f = csc(x)**5/(a + b*csc(x))
    F = -2*a**4*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(b**4*sqrt(a**2 - b**2)) + a*cot(x)*csc(x)/(2*b**2) + a*(2*a**2 + b**2)*atanh(cos(x))/(2*b**4) - cot(x)*csc(x)**2/(3*b) - (3*a**2 + 2*b**2)*cot(x)/(3*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_40():
    f = csc(x)**4/(a + b*csc(x))
    F = 2*a**3*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(b**3*sqrt(a**2 - b**2)) + a*cot(x)/b**2 - cot(x)*csc(x)/(2*b) - (2*a**2 + b**2)*atanh(cos(x))/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_41():
    f = csc(x)**3/(a + b*csc(x))
    F = -2*a**2*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(b**2*sqrt(a**2 - b**2)) + a*atanh(cos(x))/b**2 - cot(x)/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_42():
    f = csc(x)**2/(a + b*csc(x))
    F = 2*a*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(b*sqrt(a**2 - b**2)) - atanh(cos(x))/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_43():
    f = csc(x)/(a + b*csc(x))
    F = -2*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/sqrt(a**2 - b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_44():
    f = 1/(a + b*csc(c + d*x))
    F = 2*b*atanh((a + b*tan(c/2 + d*x/2))/sqrt(a**2 - b**2))/(a*d*sqrt(a**2 - b**2)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_45():
    f = sin(x)/(a + b*csc(x))
    F = -cos(x)/a - 2*b**2*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a**2*sqrt(a**2 - b**2)) - b*x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_46():
    f = sin(x)**2/(a + b*csc(x))
    F = -sin(x)*cos(x)/(2*a) + b*cos(x)/a**2 + 2*b**3*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a**3*sqrt(a**2 - b**2)) + x*(a**2 + 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_47():
    f = sin(x)**3/(a + b*csc(x))
    F = -sin(x)**2*cos(x)/(3*a) + b*sin(x)*cos(x)/(2*a**2) - (2*a**2 + 3*b**2)*cos(x)/(3*a**3) - 2*b**4*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a**4*sqrt(a**2 - b**2)) - b*x*(a**2 + 2*b**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_48():
    f = sin(x)**4/(a + b*csc(x))
    F = -sin(x)**3*cos(x)/(4*a) + b*sin(x)**2*cos(x)/(3*a**2) - (3*a**2 + 4*b**2)*sin(x)*cos(x)/(8*a**3) + b*(2*a**2 + 3*b**2)*cos(x)/(3*a**4) + 2*b**5*atanh((a + b*tan(x/2))/sqrt(a**2 - b**2))/(a**5*sqrt(a**2 - b**2)) + x*(3*a**4 + 4*a**2*b**2 + 8*b**4)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_49():
    f = (a + b*csc(c + d*x))**(-2)
    F = -b**2*cot(c + d*x)/(a*d*(a + b*csc(c + d*x))*(a**2 - b**2)) + 2*b*(2*a**2 - b**2)*atanh((a + b*tan(c/2 + d*x/2))/sqrt(a**2 - b**2))/(a**2*d*(a**2 - b**2)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_50():
    f = (a + b*csc(c + d*x))**(-3)
    F = -b**2*cot(c + d*x)/(2*a*d*(a + b*csc(c + d*x))**2*(a**2 - b**2)) - b**2*(5*a**2 - 2*b**2)*cot(c + d*x)/(2*a**2*d*(a + b*csc(c + d*x))*(a**2 - b**2)**2) + b*(6*a**4 - 5*a**2*b**2 + 2*b**4)*atanh((a + b*tan(c/2 + d*x/2))/sqrt(a**2 - b**2))/(a**3*d*(a**2 - b**2)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_51():
    f = (a + b*csc(c + d*x))**(-4)
    F = -b**2*cot(c + d*x)/(3*a*d*(a + b*csc(c + d*x))**3*(a**2 - b**2)) - b**2*(8*a**2 - 3*b**2)*cot(c + d*x)/(6*a**2*d*(a + b*csc(c + d*x))**2*(a**2 - b**2)**2) - b**2*(26*a**4 - 17*a**2*b**2 + 6*b**4)*cot(c + d*x)/(6*a**3*d*(a + b*csc(c + d*x))*(a**2 - b**2)**3) + b*(8*a**6 - 8*a**4*b**2 + 7*a**2*b**4 - 2*b**6)*atanh((a + b*tan(c/2 + d*x/2))/sqrt(a**2 - b**2))/(a**4*d*(a**2 - b**2)**(sympy.S(7)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_52():
    f = 1/(5*csc(c + d*x) + 3)
    F = -x/12 - 5*atan(cos(c + d*x)/(sin(c + d*x) + 3))/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_53():
    f = 1/(3*csc(c + d*x) + 5)
    F = x/5 + 3*log(sin(c/2 + d*x/2) + 3*cos(c/2 + d*x/2))/(20*d) - 3*log(3*sin(c/2 + d*x/2) + cos(c/2 + d*x/2))/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_54():
    f = (a + b*csc(e + f*x))**m*csc(e + f*x)**3
    F = sqrt(2)*a*(a + b)*(a + b*csc(e + f*x))**m*cot(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - csc(e + f*x)/2, b*(1 - csc(e + f*x))/(a + b))/(b**2*f*((a + b*csc(e + f*x))/(a + b))**m*(m + 2)*sqrt(csc(e + f*x) + 1)) - (a + b*csc(e + f*x))**(m + 1)*cot(e + f*x)/(b*f*(m + 2)) - sqrt(2)*(a + b*csc(e + f*x))**m*(a**2 + b**2*(m + 1))*cot(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - csc(e + f*x)/2, b*(1 - csc(e + f*x))/(a + b))/(b**2*f*((a + b*csc(e + f*x))/(a + b))**m*(m + 2)*sqrt(csc(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_55():
    f = (a + b*csc(e + f*x))**m*csc(e + f*x)**2
    F = sqrt(2)*a*(a + b*csc(e + f*x))**m*cot(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - csc(e + f*x)/2, b*(1 - csc(e + f*x))/(a + b))/(b*f*((a + b*csc(e + f*x))/(a + b))**m*sqrt(csc(e + f*x) + 1)) - sqrt(2)*(a + b)*(a + b*csc(e + f*x))**m*cot(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - csc(e + f*x)/2, b*(1 - csc(e + f*x))/(a + b))/(b*f*((a + b*csc(e + f*x))/(a + b))**m*sqrt(csc(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_56():
    f = (a + b*csc(e + f*x))**m*csc(e + f*x)
    F = -sqrt(2)*(a + b*csc(e + f*x))**m*cot(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - csc(e + f*x)/2, b*(1 - csc(e + f*x))/(a + b))/(f*((a + b*csc(e + f*x))/(a + b))**m*sqrt(csc(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_57():
    f = (a + b*csc(e + f*x))**m
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.csc((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_58():
    f = (a + b*csc(e + f*x))**m*sin(e + f*x)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.csc((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * sympy.sin((Symbol('e') + (Symbol('f') * x)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_6_Cosecant_4_6_1_2_d_csc_pow_n_a_plus_b_csc_pow_m_59():
    f = (a + b*csc(e + f*x))**m*sin(e + f*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.csc((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')) * (sympy.sin((Symbol('e') + (Symbol('f') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F

