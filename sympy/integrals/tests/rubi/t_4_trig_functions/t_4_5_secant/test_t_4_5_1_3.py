"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.1.3 (d sin)^n (a+b sec)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_1():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**9
    F = -a*log(cos(c + d*x))/d - a*cos(c + d*x)**9/(9*d) - a*cos(c + d*x)**8/(8*d) + 4*a*cos(c + d*x)**7/(7*d) + 2*a*cos(c + d*x)**6/(3*d) - 6*a*cos(c + d*x)**5/(5*d) - 3*a*cos(c + d*x)**4/(2*d) + 4*a*cos(c + d*x)**3/(3*d) + 2*a*cos(c + d*x)**2/d - a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_2():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**7
    F = -a*log(cos(c + d*x))/d + a*cos(c + d*x)**7/(7*d) + a*cos(c + d*x)**6/(6*d) - 3*a*cos(c + d*x)**5/(5*d) - 3*a*cos(c + d*x)**4/(4*d) + a*cos(c + d*x)**3/d + 3*a*cos(c + d*x)**2/(2*d) - a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_3():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**5
    F = -a*log(cos(c + d*x))/d - a*cos(c + d*x)**5/(5*d) - a*cos(c + d*x)**4/(4*d) + 2*a*cos(c + d*x)**3/(3*d) + a*cos(c + d*x)**2/d - a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_4():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**3
    F = -a*log(cos(c + d*x))/d + a*cos(c + d*x)**3/(3*d) + a*cos(c + d*x)**2/(2*d) - a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_5():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)
    F = -a*log(cos(c + d*x))/d - a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_6():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)
    F = a*log(1 - cos(c + d*x))/d - a*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_7():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**3
    F = -a**2/(2*d*(-a*cos(c + d*x) + a)) + 3*a*log(1 - cos(c + d*x))/(4*d) + a*log(cos(c + d*x) + 1)/(4*d) - a*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_8():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**5
    F = -a**3/(8*d*(-a*cos(c + d*x) + a)**2) - a**2/(8*d*(a*cos(c + d*x) + a)) - a**2/(2*d*(-a*cos(c + d*x) + a)) + 11*a*log(1 - cos(c + d*x))/(16*d) + 5*a*log(cos(c + d*x) + 1)/(16*d) - a*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_9():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**7
    F = -a**4/(24*d*(-a*cos(c + d*x) + a)**3) - a**3/(32*d*(a*cos(c + d*x) + a)**2) - 5*a**3/(32*d*(-a*cos(c + d*x) + a)**2) - 3*a**2/(16*d*(a*cos(c + d*x) + a)) - a**2/(2*d*(-a*cos(c + d*x) + a)) + 21*a*log(1 - cos(c + d*x))/(32*d) + 11*a*log(cos(c + d*x) + 1)/(32*d) - a*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_10():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**8
    F = 35*a*x/128 - a*sin(c + d*x)**7*cos(c + d*x)/(8*d) - a*sin(c + d*x)**7/(7*d) - 7*a*sin(c + d*x)**5*cos(c + d*x)/(48*d) - a*sin(c + d*x)**5/(5*d) - 35*a*sin(c + d*x)**3*cos(c + d*x)/(192*d) - a*sin(c + d*x)**3/(3*d) - 35*a*sin(c + d*x)*cos(c + d*x)/(128*d) - a*sin(c + d*x)/d + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_11():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**6
    F = 5*a*x/16 - a*sin(c + d*x)**5*cos(c + d*x)/(6*d) - a*sin(c + d*x)**5/(5*d) - 5*a*sin(c + d*x)**3*cos(c + d*x)/(24*d) - a*sin(c + d*x)**3/(3*d) - 5*a*sin(c + d*x)*cos(c + d*x)/(16*d) - a*sin(c + d*x)/d + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_12():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**4
    F = 3*a*x/8 - a*sin(c + d*x)**3*cos(c + d*x)/(4*d) - a*sin(c + d*x)**3/(3*d) - 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - a*sin(c + d*x)/d + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_13():
    f = (a*sec(c + d*x) + a)*sin(c + d*x)**2
    F = a*x/2 - a*sin(c + d*x)*cos(c + d*x)/(2*d) - a*sin(c + d*x)/d + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_14():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**2
    F = -a*cot(c + d*x)/d + a*atanh(sin(c + d*x))/d - a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_15():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**4
    F = -a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + a*atanh(sin(c + d*x))/d - a*csc(c + d*x)**3/(3*d) - a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_16():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**6
    F = -a*cot(c + d*x)**5/(5*d) - 2*a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + a*atanh(sin(c + d*x))/d - a*csc(c + d*x)**5/(5*d) - a*csc(c + d*x)**3/(3*d) - a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_17():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**8
    F = -a*cot(c + d*x)**7/(7*d) - 3*a*cot(c + d*x)**5/(5*d) - a*cot(c + d*x)**3/d - a*cot(c + d*x)/d + a*atanh(sin(c + d*x))/d - a*csc(c + d*x)**7/(7*d) - a*csc(c + d*x)**5/(5*d) - a*csc(c + d*x)**3/(3*d) - a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_18():
    f = (a*sec(c + d*x) + a)*csc(c + d*x)**10
    F = -a*cot(c + d*x)**9/(9*d) - 4*a*cot(c + d*x)**7/(7*d) - 6*a*cot(c + d*x)**5/(5*d) - 4*a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + a*atanh(sin(c + d*x))/d - a*csc(c + d*x)**9/(9*d) - a*csc(c + d*x)**7/(7*d) - a*csc(c + d*x)**5/(5*d) - a*csc(c + d*x)**3/(3*d) - a*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_19():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**9
    F = -2*a**2*log(cos(c + d*x))/d - a**2*cos(c + d*x)**9/(9*d) - a**2*cos(c + d*x)**8/(4*d) + 3*a**2*cos(c + d*x)**7/(7*d) + 4*a**2*cos(c + d*x)**6/(3*d) - 2*a**2*cos(c + d*x)**5/(5*d) - 3*a**2*cos(c + d*x)**4/d - 2*a**2*cos(c + d*x)**3/(3*d) + 4*a**2*cos(c + d*x)**2/d + 3*a**2*cos(c + d*x)/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_20():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**7
    F = -2*a**2*log(cos(c + d*x))/d + a**2*cos(c + d*x)**7/(7*d) + a**2*cos(c + d*x)**6/(3*d) - 2*a**2*cos(c + d*x)**5/(5*d) - 3*a**2*cos(c + d*x)**4/(2*d) + 3*a**2*cos(c + d*x)**2/d + 2*a**2*cos(c + d*x)/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_21():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**5
    F = -2*a**2*log(cos(c + d*x))/d - a**2*cos(c + d*x)**5/(5*d) - a**2*cos(c + d*x)**4/(2*d) + a**2*cos(c + d*x)**3/(3*d) + 2*a**2*cos(c + d*x)**2/d + a**2*cos(c + d*x)/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_22():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**3
    F = -2*a**2*log(cos(c + d*x))/d + a**2*cos(c + d*x)**3/(3*d) + a**2*cos(c + d*x)**2/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_23():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)
    F = -2*a**2*log(cos(c + d*x))/d - a**2*cos(c + d*x)/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_24():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)
    F = 2*a**2*log(1 - cos(c + d*x))/d - 2*a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_25():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**3
    F = -a**3/(d*(-a*cos(c + d*x) + a)) + 2*a**2*log(1 - cos(c + d*x))/d - 2*a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_26():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**5
    F = -a**4/(4*d*(-a*cos(c + d*x) + a)**2) - 5*a**3/(4*d*(-a*cos(c + d*x) + a)) + 17*a**2*log(1 - cos(c + d*x))/(8*d) - a**2*log(cos(c + d*x) + 1)/(8*d) - 2*a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_27():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**7
    F = -a**5/(12*d*(-a*cos(c + d*x) + a)**3) - 3*a**4/(8*d*(-a*cos(c + d*x) + a)**2) + a**3/(16*d*(a*cos(c + d*x) + a)) - 23*a**3/(16*d*(-a*cos(c + d*x) + a)) + 9*a**2*log(1 - cos(c + d*x))/(4*d) - a**2*log(cos(c + d*x) + 1)/(4*d) - 2*a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_28():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**9
    F = -a**6/(32*d*(-a*cos(c + d*x) + a)**4) - 7*a**5/(48*d*(-a*cos(c + d*x) + a)**3) + a**4/(64*d*(a*cos(c + d*x) + a)**2) - 15*a**4/(32*d*(-a*cos(c + d*x) + a)**2) + 9*a**3/(64*d*(a*cos(c + d*x) + a)) - 51*a**3/(32*d*(-a*cos(c + d*x) + a)) + 303*a**2*log(1 - cos(c + d*x))/(128*d) - 47*a**2*log(cos(c + d*x) + 1)/(128*d) - 2*a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_29():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**8
    F = -245*a**2*x/128 - 2*a**2*sin(c + d*x)**7/(7*d) - 2*a**2*sin(c + d*x)**5/(5*d) - 2*a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)*cos(c + d*x)**7/(8*d) - 17*a**2*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 11*a**2*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 139*a**2*sin(c + d*x)*cos(c + d*x)/(128*d) - 2*a**2*sin(c + d*x)/d + a**2*tan(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_30():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**6
    F = -25*a**2*x/16 - 2*a**2*sin(c + d*x)**5/(5*d) - 2*a**2*sin(c + d*x)**3/(3*d) - a**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 7*a**2*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 7*a**2*sin(c + d*x)*cos(c + d*x)/(16*d) - 2*a**2*sin(c + d*x)/d + a**2*tan(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_31():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**4
    F = -9*a**2*x/8 - 2*a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) - a**2*sin(c + d*x)*cos(c + d*x)/(8*d) - 2*a**2*sin(c + d*x)/d + a**2*tan(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_32():
    f = (a*sec(c + d*x) + a)**2*sin(c + d*x)**2
    F = -a**2*x/2 - a**2*sin(c + d*x)*cos(c + d*x)/(2*d) - 2*a**2*sin(c + d*x)/d + a**2*tan(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_33():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**2
    F = a**2*tan(c + d*x)/d - 2*a**2*cot(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d - 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_34():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**4
    F = -a**4*tan(c + d*x)/(3*d*(-a*cos(c + d*x) + a)**2) + 10*a**2*tan(c + d*x)/(3*d) + 2*a**2*atanh(sin(c + d*x))/d - 2*a**2*tan(c + d*x)/(d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_35():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**6
    F = a**2*tan(c + d*x)/d - 2*a**2*cot(c + d*x)**5/(5*d) - 5*a**2*cot(c + d*x)**3/(3*d) - 4*a**2*cot(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d - 2*a**2*csc(c + d*x)**5/(5*d) - 2*a**2*csc(c + d*x)**3/(3*d) - 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_36():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**8
    F = a**2*tan(c + d*x)/d - 2*a**2*cot(c + d*x)**7/(7*d) - 7*a**2*cot(c + d*x)**5/(5*d) - 3*a**2*cot(c + d*x)**3/d - 5*a**2*cot(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d - 2*a**2*csc(c + d*x)**7/(7*d) - 2*a**2*csc(c + d*x)**5/(5*d) - 2*a**2*csc(c + d*x)**3/(3*d) - 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_37():
    f = (a*sec(c + d*x) + a)**2*csc(c + d*x)**10
    F = a**2*tan(c + d*x)/d - 2*a**2*cot(c + d*x)**9/(9*d) - 9*a**2*cot(c + d*x)**7/(7*d) - 16*a**2*cot(c + d*x)**5/(5*d) - 14*a**2*cot(c + d*x)**3/(3*d) - 6*a**2*cot(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d - 2*a**2*csc(c + d*x)**9/(9*d) - 2*a**2*csc(c + d*x)**7/(7*d) - 2*a**2*csc(c + d*x)**5/(5*d) - 2*a**2*csc(c + d*x)**3/(3*d) - 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_38():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**9
    F = a**3*log(cos(c + d*x))/d - a**3*cos(c + d*x)**9/(9*d) - 3*a**3*cos(c + d*x)**8/(8*d) + a**3*cos(c + d*x)**7/(7*d) + 11*a**3*cos(c + d*x)**6/(6*d) + 6*a**3*cos(c + d*x)**5/(5*d) - 7*a**3*cos(c + d*x)**4/(2*d) - 14*a**3*cos(c + d*x)**3/(3*d) + 3*a**3*cos(c + d*x)**2/d + 11*a**3*cos(c + d*x)/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_39():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**7
    F = a**3*cos(c + d*x)**7/(7*d) + a**3*cos(c + d*x)**6/(2*d) - 2*a**3*cos(c + d*x)**4/d - 2*a**3*cos(c + d*x)**3/d + 3*a**3*cos(c + d*x)**2/d + 8*a**3*cos(c + d*x)/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_40():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**5
    F = -a**3*log(cos(c + d*x))/d - a**3*cos(c + d*x)**5/(5*d) - 3*a**3*cos(c + d*x)**4/(4*d) - a**3*cos(c + d*x)**3/(3*d) + 5*a**3*cos(c + d*x)**2/(2*d) + 5*a**3*cos(c + d*x)/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_41():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**3
    F = -2*a**3*log(cos(c + d*x))/d + a**3*cos(c + d*x)**3/(3*d) + 3*a**3*cos(c + d*x)**2/(2*d) + 2*a**3*cos(c + d*x)/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_42():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)
    F = -3*a**3*log(cos(c + d*x))/d - a**3*cos(c + d*x)/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_43():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)
    F = 4*a**3*log(1 - cos(c + d*x))/d - 4*a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_44():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**3
    F = -2*a**4/(d*(-a*cos(c + d*x) + a)) + 5*a**3*log(1 - cos(c + d*x))/d - 5*a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_45():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**5
    F = -a**5/(2*d*(-a*cos(c + d*x) + a)**2) - 3*a**4/(d*(-a*cos(c + d*x) + a)) + 6*a**3*log(1 - cos(c + d*x))/d - 6*a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_46():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**7
    F = -a**6/(6*d*(-a*cos(c + d*x) + a)**3) - 7*a**5/(8*d*(-a*cos(c + d*x) + a)**2) - 31*a**4/(8*d*(-a*cos(c + d*x) + a)) + 111*a**3*log(1 - cos(c + d*x))/(16*d) + a**3*log(cos(c + d*x) + 1)/(16*d) - 7*a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_47():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**9
    F = -a**7/(16*d*(-a*cos(c + d*x) + a)**4) - a**6/(3*d*(-a*cos(c + d*x) + a)**3) - 39*a**5/(32*d*(-a*cos(c + d*x) + a)**2) - a**4/(32*d*(a*cos(c + d*x) + a)) - 75*a**4/(16*d*(-a*cos(c + d*x) + a)) + 501*a**3*log(1 - cos(c + d*x))/(64*d) + 11*a**3*log(cos(c + d*x) + 1)/(64*d) - 8*a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_48():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**8
    F = -805*a**3*x/128 - 3*a**3*sin(c + d*x)**7/(7*d) - 2*a**3*sin(c + d*x)**5/(5*d) - a**3*sin(c + d*x)**3/(3*d) + a**3*sin(c + d*x)*cos(c + d*x)**7/(8*d) - a**3*sin(c + d*x)*cos(c + d*x)**5/(48*d) - 293*a**3*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 603*a**3*sin(c + d*x)*cos(c + d*x)/(128*d) + a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 3*a**3*tan(c + d*x)/d - a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_49():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**6
    F = -85*a**3*x/16 - 3*a**3*sin(c + d*x)**5/(5*d) - 2*a**3*sin(c + d*x)**3/(3*d) - a**3*sin(c + d*x)*cos(c + d*x)**5/(6*d) - 5*a**3*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 43*a**3*sin(c + d*x)*cos(c + d*x)/(16*d) - a**3*sin(c + d*x)/d + a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 3*a**3*tan(c + d*x)/d + a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_50():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**4
    F = -33*a**3*x/8 - a**3*sin(c + d*x)**3/d + a**3*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 7*a**3*sin(c + d*x)*cos(c + d*x)/(8*d) - 2*a**3*sin(c + d*x)/d + a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 3*a**3*tan(c + d*x)/d + 3*a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_51():
    f = (a*sec(c + d*x) + a)**3*sin(c + d*x)**2
    F = -5*a**3*x/2 - a**3*sin(c + d*x)*cos(c + d*x)/(2*d) - 3*a**3*sin(c + d*x)/d + a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 3*a**3*tan(c + d*x)/d + 5*a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_52():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**2
    F = a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 3*a**3*tan(c + d*x)/d + 9*a**3*atanh(sin(c + d*x))/(2*d) - 4*a**3*sin(c + d*x)/(d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_53():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**4
    F = a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 3*a**3*tan(c + d*x)/d + 11*a**3*atanh(sin(c + d*x))/(2*d) - 17*a**3*sin(c + d*x)/(3*d*(1 - cos(c + d*x))) - 2*a**3*sin(c + d*x)/(3*d*(1 - cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_54():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**6
    F = -76*a**6*tan(c + d*x)*sec(c + d*x)/(15*d*(-a**3*cos(c + d*x) + a**3)) - a**6*tan(c + d*x)*sec(c + d*x)/(5*d*(-a*cos(c + d*x) + a)**3) - 11*a**5*tan(c + d*x)*sec(c + d*x)/(15*d*(-a*cos(c + d*x) + a)**2) + 13*a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 152*a**3*tan(c + d*x)/(15*d) + 13*a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_55():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**8
    F = 3*a**3*tan(c + d*x)/d - 4*a**3*cot(c + d*x)**7/(7*d) - 3*a**3*cot(c + d*x)**5/d - 7*a**3*cot(c + d*x)**3/d - 13*a**3*cot(c + d*x)/d + 15*a**3*atanh(sin(c + d*x))/(2*d) + a**3*csc(c + d*x)**7*sec(c + d*x)**2/(2*d) - 15*a**3*csc(c + d*x)**7/(14*d) - 3*a**3*csc(c + d*x)**5/(2*d) - 5*a**3*csc(c + d*x)**3/(2*d) - 15*a**3*csc(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_56():
    f = (a*sec(c + d*x) + a)**3*csc(c + d*x)**10
    F = 3*a**3*tan(c + d*x)/d - 4*a**3*cot(c + d*x)**9/(9*d) - 19*a**3*cot(c + d*x)**7/(7*d) - 36*a**3*cot(c + d*x)**5/(5*d) - 34*a**3*cot(c + d*x)**3/(3*d) - 16*a**3*cot(c + d*x)/d + 17*a**3*atanh(sin(c + d*x))/(2*d) + a**3*csc(c + d*x)**9*sec(c + d*x)**2/(2*d) - 17*a**3*csc(c + d*x)**9/(18*d) - 17*a**3*csc(c + d*x)**7/(14*d) - 17*a**3*csc(c + d*x)**5/(10*d) - 17*a**3*csc(c + d*x)**3/(6*d) - 17*a**3*csc(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_57():
    f = sin(c + d*x)**9/(a*sec(c + d*x) + a)
    F = sin(c + d*x)**8/(8*a*d) - cos(c + d*x)**9/(9*a*d) + 3*cos(c + d*x)**7/(7*a*d) - 3*cos(c + d*x)**5/(5*a*d) + cos(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_58():
    f = sin(c + d*x)**7/(a*sec(c + d*x) + a)
    F = sin(c + d*x)**6/(6*a*d) + cos(c + d*x)**7/(7*a*d) - 2*cos(c + d*x)**5/(5*a*d) + cos(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_59():
    f = sin(c + d*x)**5/(a*sec(c + d*x) + a)
    F = sin(c + d*x)**4/(4*a*d) - cos(c + d*x)**5/(5*a*d) + cos(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_60():
    f = sin(c + d*x)**3/(a*sec(c + d*x) + a)
    F = sin(c + d*x)**2/(2*a*d) + cos(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_61():
    f = sin(c + d*x)/(a*sec(c + d*x) + a)
    F = log(cos(c + d*x) + 1)/(a*d) - cos(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_62():
    f = csc(c + d*x)/(a*sec(c + d*x) + a)
    F = cot(c + d*x)*csc(c + d*x)/(2*a*d) - atanh(cos(c + d*x))/(2*a*d) - csc(c + d*x)**2/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_63():
    f = csc(c + d*x)**3/(a*sec(c + d*x) + a)
    F = cot(c + d*x)*csc(c + d*x)**3/(4*a*d) - cot(c + d*x)*csc(c + d*x)/(8*a*d) - atanh(cos(c + d*x))/(8*a*d) - csc(c + d*x)**4/(4*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_64():
    f = csc(c + d*x)**5/(a*sec(c + d*x) + a)
    F = cot(c + d*x)*csc(c + d*x)**5/(6*a*d) - cot(c + d*x)*csc(c + d*x)**3/(24*a*d) - cot(c + d*x)*csc(c + d*x)/(16*a*d) - atanh(cos(c + d*x))/(16*a*d) - csc(c + d*x)**6/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_65():
    f = sin(c + d*x)**8/(a*sec(c + d*x) + a)
    F = -5*x/(128*a) + sin(c + d*x)**7/(7*a*d) + sin(c + d*x)**5*cos(c + d*x)**3/(8*a*d) + 5*sin(c + d*x)**3*cos(c + d*x)**3/(48*a*d) + 5*sin(c + d*x)*cos(c + d*x)**3/(64*a*d) - 5*sin(c + d*x)*cos(c + d*x)/(128*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_66():
    f = sin(c + d*x)**6/(a*sec(c + d*x) + a)
    F = -x/(16*a) + sin(c + d*x)**5/(5*a*d) + sin(c + d*x)**3*cos(c + d*x)**3/(6*a*d) + sin(c + d*x)*cos(c + d*x)**3/(8*a*d) - sin(c + d*x)*cos(c + d*x)/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_67():
    f = sin(c + d*x)**4/(a*sec(c + d*x) + a)
    F = -x/(8*a) + sin(c + d*x)**3/(3*a*d) + sin(c + d*x)*cos(c + d*x)**3/(4*a*d) - sin(c + d*x)*cos(c + d*x)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_68():
    f = sin(c + d*x)**2/(a*sec(c + d*x) + a)
    F = -x/(2*a) - sin(c + d*x)*cos(c + d*x)/(2*a*d) + sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_69():
    f = csc(c + d*x)**2/(a*sec(c + d*x) + a)
    F = cot(c + d*x)**3/(3*a*d) - csc(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_70():
    f = csc(c + d*x)**4/(a*sec(c + d*x) + a)
    F = cot(c + d*x)**5/(5*a*d) + cot(c + d*x)**3/(3*a*d) - csc(c + d*x)**5/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_71():
    f = csc(c + d*x)**6/(a*sec(c + d*x) + a)
    F = cot(c + d*x)**7/(7*a*d) + 2*cot(c + d*x)**5/(5*a*d) + cot(c + d*x)**3/(3*a*d) - csc(c + d*x)**7/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_72():
    f = csc(c + d*x)**8/(a*sec(c + d*x) + a)
    F = cot(c + d*x)**9/(9*a*d) + 3*cot(c + d*x)**7/(7*a*d) + 3*cot(c + d*x)**5/(5*a*d) + cot(c + d*x)**3/(3*a*d) - csc(c + d*x)**9/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_73():
    f = csc(c + d*x)**10/(a*sec(c + d*x) + a)
    F = cot(c + d*x)**11/(11*a*d) + 4*cot(c + d*x)**9/(9*a*d) + 6*cot(c + d*x)**7/(7*a*d) + 4*cot(c + d*x)**5/(5*a*d) + cot(c + d*x)**3/(3*a*d) - csc(c + d*x)**11/(11*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_74():
    f = sin(c + d*x)**11/(a*sec(c + d*x) + a)**2
    F = 4*(-a*cos(c + d*x) + a)**6/(3*a**8*d) - 4*(-a*cos(c + d*x) + a)**7/(a**9*d) + 19*(-a*cos(c + d*x) + a)**8/(4*a**10*d) - 25*(-a*cos(c + d*x) + a)**9/(9*a**11*d) + 4*(-a*cos(c + d*x) + a)**10/(5*a**12*d) - (-a*cos(c + d*x) + a)**11/(11*a**13*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_75():
    f = sin(c + d*x)**9/(a*sec(c + d*x) + a)**2
    F = 4*(-a*cos(c + d*x) + a)**5/(5*a**7*d) - 2*(-a*cos(c + d*x) + a)**6/(a**8*d) + 13*(-a*cos(c + d*x) + a)**7/(7*a**9*d) - 3*(-a*cos(c + d*x) + a)**8/(4*a**10*d) + (-a*cos(c + d*x) + a)**9/(9*a**11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_76():
    f = sin(c + d*x)**7/(a*sec(c + d*x) + a)**2
    F = cos(c + d*x)**7/(7*a**2*d) - cos(c + d*x)**6/(3*a**2*d) + cos(c + d*x)**4/(2*a**2*d) - cos(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_77():
    f = sin(c + d*x)**5/(a*sec(c + d*x) + a)**2
    F = -cos(c + d*x)**5/(5*a**2*d) + cos(c + d*x)**4/(2*a**2*d) - cos(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_78():
    f = sin(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = -2*log(cos(c + d*x) + 1)/(a**2*d) + cos(c + d*x)**3/(3*a**2*d) - cos(c + d*x)**2/(a**2*d) + 2*cos(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_79():
    f = sin(c + d*x)/(a*sec(c + d*x) + a)**2
    F = 1/(d*(a**2*cos(c + d*x) + a**2)) + 2*log(cos(c + d*x) + 1)/(a**2*d) - cos(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_80():
    f = csc(c + d*x)/(a*sec(c + d*x) + a)**2
    F = -3/(4*d*(a**2*cos(c + d*x) + a**2)) + 1/(4*d*(a*cos(c + d*x) + a)**2) - atanh(cos(c + d*x))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_81():
    f = csc(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = -(2*a*cos(c + d*x) + a)/(6*d*(1 - cos(c + d*x))*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_82():
    f = csc(c + d*x)**5/(a*sec(c + d*x) + a)**2
    F = a**2/(32*d*(a*cos(c + d*x) + a)**4) - a/(48*d*(a*cos(c + d*x) + a)**3) - 1/(32*d*(a**2*cos(c + d*x) + a**2)) - 1/(64*d*(-a**2*cos(c + d*x) + a**2)) - 1/(32*d*(a*cos(c + d*x) + a)**2) - 1/(64*d*(-a*cos(c + d*x) + a)**2) + atanh(cos(c + d*x))/(64*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_83():
    f = sin(c + d*x)**8/(a*sec(c + d*x) + a)**2
    F = 11*x/(128*a**2) + 2*sin(c + d*x)**7/(7*a**2*d) - 2*sin(c + d*x)**5/(5*a**2*d) - sin(c + d*x)**3*cos(c + d*x)**5/(8*a**2*d) - sin(c + d*x)**3*cos(c + d*x)**3/(6*a**2*d) - sin(c + d*x)*cos(c + d*x)**5/(16*a**2*d) - 7*sin(c + d*x)*cos(c + d*x)**3/(64*a**2*d) + 11*sin(c + d*x)*cos(c + d*x)/(128*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_84():
    f = sin(c + d*x)**6/(a*sec(c + d*x) + a)**2
    F = 3*x/(16*a**2) - sin(c + d*x)**5/(10*a**2*d) - sin(c + d*x)**3*cos(c + d*x)/(8*a**2*d) - 3*sin(c + d*x)*cos(c + d*x)/(16*a**2*d) - (-a*cos(c + d*x) + a)**3*sin(c + d*x)**3/(6*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_85():
    f = sin(c + d*x)**4/(a*sec(c + d*x) + a)**2
    F = 7*x/(8*a**2) + 2*sin(c + d*x)**3/(3*a**2*d) + sin(c + d*x)*cos(c + d*x)**3/(4*a**2*d) + 7*sin(c + d*x)*cos(c + d*x)/(8*a**2*d) - 2*sin(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_86():
    f = sin(c + d*x)**2/(a*sec(c + d*x) + a)**2
    F = -5*x/(2*a**2) - sin(c + d*x)*cos(c + d*x)/(2*a**2*d) + 2*sin(c + d*x)/(a**2*d) + 2*sin(c + d*x)/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_87():
    f = csc(c + d*x)**2/(a*sec(c + d*x) + a)**2
    F = -2*cot(c + d*x)**5/(5*a**2*d) - cot(c + d*x)**3/(3*a**2*d) + 2*csc(c + d*x)**5/(5*a**2*d) - 2*csc(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_88():
    f = csc(c + d*x)**4/(a*sec(c + d*x) + a)**2
    F = -2*cot(c + d*x)**7/(7*a**2*d) - 3*cot(c + d*x)**5/(5*a**2*d) - cot(c + d*x)**3/(3*a**2*d) + 2*csc(c + d*x)**7/(7*a**2*d) - 2*csc(c + d*x)**5/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_89():
    f = csc(c + d*x)**6/(a*sec(c + d*x) + a)**2
    F = -2*cot(c + d*x)**9/(9*a**2*d) - 5*cot(c + d*x)**7/(7*a**2*d) - 4*cot(c + d*x)**5/(5*a**2*d) - cot(c + d*x)**3/(3*a**2*d) + 2*csc(c + d*x)**9/(9*a**2*d) - 2*csc(c + d*x)**7/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_90():
    f = csc(c + d*x)**8/(a*sec(c + d*x) + a)**2
    F = -2*cot(c + d*x)**11/(11*a**2*d) - 7*cot(c + d*x)**9/(9*a**2*d) - 9*cot(c + d*x)**7/(7*a**2*d) - cot(c + d*x)**5/(a**2*d) - cot(c + d*x)**3/(3*a**2*d) + 2*csc(c + d*x)**11/(11*a**2*d) - 2*csc(c + d*x)**9/(9*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_91():
    f = sin(c + d*x)**11/(a*sec(c + d*x) + a)**3
    F = 2*(-a*cos(c + d*x) + a)**6/(3*a**9*d) - 16*(-a*cos(c + d*x) + a)**7/(7*a**10*d) + 25*(-a*cos(c + d*x) + a)**8/(8*a**11*d) - 19*(-a*cos(c + d*x) + a)**9/(9*a**12*d) + 7*(-a*cos(c + d*x) + a)**10/(10*a**13*d) - (-a*cos(c + d*x) + a)**11/(11*a**14*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_92():
    f = sin(c + d*x)**9/(a*sec(c + d*x) + a)**3
    F = -cos(c + d*x)**9/(9*a**3*d) + 3*cos(c + d*x)**8/(8*a**3*d) - 2*cos(c + d*x)**7/(7*a**3*d) - cos(c + d*x)**6/(3*a**3*d) + 3*cos(c + d*x)**5/(5*a**3*d) - cos(c + d*x)**4/(4*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_93():
    f = sin(c + d*x)**7/(a*sec(c + d*x) + a)**3
    F = cos(c + d*x)**7/(7*a**3*d) - cos(c + d*x)**6/(2*a**3*d) + 3*cos(c + d*x)**5/(5*a**3*d) - cos(c + d*x)**4/(4*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_94():
    f = sin(c + d*x)**5/(a*sec(c + d*x) + a)**3
    F = 4*log(cos(c + d*x) + 1)/(a**3*d) - cos(c + d*x)**5/(5*a**3*d) + 3*cos(c + d*x)**4/(4*a**3*d) - 4*cos(c + d*x)**3/(3*a**3*d) + 2*cos(c + d*x)**2/(a**3*d) - 4*cos(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_95():
    f = sin(c + d*x)**3/(a*sec(c + d*x) + a)**3
    F = -2/(d*(a**3*cos(c + d*x) + a**3)) - 7*log(cos(c + d*x) + 1)/(a**3*d) + cos(c + d*x)**3/(3*a**3*d) - 3*cos(c + d*x)**2/(2*a**3*d) + 5*cos(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_96():
    f = sin(c + d*x)/(a*sec(c + d*x) + a)**3
    F = 3/(d*(a**3*cos(c + d*x) + a**3)) - 1/(2*a*d*(a*cos(c + d*x) + a)**2) + 3*log(cos(c + d*x) + 1)/(a**3*d) - cos(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_97():
    f = csc(c + d*x)/(a*sec(c + d*x) + a)**3
    F = -7/(8*d*(a**3*cos(c + d*x) + a**3)) - 1/(6*d*(a*cos(c + d*x) + a)**3) + 5/(8*a*d*(a*cos(c + d*x) + a)**2) - atanh(cos(c + d*x))/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_98():
    f = csc(c + d*x)**3/(a*sec(c + d*x) + a)**3
    F = -a/(16*d*(a*cos(c + d*x) + a)**4) - 1/(16*d*(a**3*cos(c + d*x) + a**3)) - 1/(32*d*(-a**3*cos(c + d*x) + a**3)) + 1/(6*d*(a*cos(c + d*x) + a)**3) - 3/(32*a*d*(a*cos(c + d*x) + a)**2) + atanh(cos(c + d*x))/(32*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_99():
    f = csc(c + d*x)**5/(a*sec(c + d*x) + a)**3
    F = -a**2/(40*d*(a*cos(c + d*x) + a)**5) + 3*a/(64*d*(a*cos(c + d*x) + a)**4) - 3/(128*d*(a**3*cos(c + d*x) + a**3)) - 1/(64*a*d*(a*cos(c + d*x) + a)**2) - 1/(128*a*d*(-a*cos(c + d*x) + a)**2) + 3*atanh(cos(c + d*x))/(128*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_100():
    f = sin(c + d*x)**8/(a*sec(c + d*x) + a)**3
    F = -29*x/(128*a**3) + 3*sin(c + d*x)**7/(7*a**3*d) - 7*sin(c + d*x)**5/(5*a**3*d) + 4*sin(c + d*x)**3/(3*a**3*d) + sin(c + d*x)*cos(c + d*x)**7/(8*a**3*d) + 23*sin(c + d*x)*cos(c + d*x)**5/(48*a**3*d) - 29*sin(c + d*x)*cos(c + d*x)**3/(192*a**3*d) - 29*sin(c + d*x)*cos(c + d*x)/(128*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_101():
    f = sin(c + d*x)**6/(a*sec(c + d*x) + a)**3
    F = -23*x/(16*a**3) + 3*sin(c + d*x)**5/(5*a**3*d) - 7*sin(c + d*x)**3/(3*a**3*d) - sin(c + d*x)*cos(c + d*x)**5/(6*a**3*d) - 23*sin(c + d*x)*cos(c + d*x)**3/(24*a**3*d) - 23*sin(c + d*x)*cos(c + d*x)/(16*a**3*d) + 4*sin(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_102():
    f = sin(c + d*x)**4/(a*sec(c + d*x) + a)**3
    F = 51*x/(8*a**3) + sin(c + d*x)**3/(a**3*d) + sin(c + d*x)*cos(c + d*x)**3/(4*a**3*d) + 19*sin(c + d*x)*cos(c + d*x)/(8*a**3*d) - 7*sin(c + d*x)/(a**3*d) - 4*sin(c + d*x)/(a**3*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_103():
    f = sin(c + d*x)**2/(a*sec(c + d*x) + a)**3
    F = -11*x/(2*a**3) - sin(c + d*x)*cos(c + d*x)/(2*a**3*d) + 3*sin(c + d*x)/(a**3*d) + 19*sin(c + d*x)/(3*a**3*d*(cos(c + d*x) + 1)) - 2*sin(c + d*x)/(3*a**3*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_104():
    f = csc(c + d*x)**2/(a*sec(c + d*x) + a)**3
    F = 4*cot(c + d*x)**7/(7*a**3*d) + 3*cot(c + d*x)**5/(5*a**3*d) - 4*csc(c + d*x)**7/(7*a**3*d) + 7*csc(c + d*x)**5/(5*a**3*d) - csc(c + d*x)**3/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_105():
    f = csc(c + d*x)**4/(a*sec(c + d*x) + a)**3
    F = 4*cot(c + d*x)**9/(9*a**3*d) + cot(c + d*x)**7/(a**3*d) + 3*cot(c + d*x)**5/(5*a**3*d) - 4*csc(c + d*x)**9/(9*a**3*d) + csc(c + d*x)**7/(a**3*d) - 3*csc(c + d*x)**5/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_106():
    f = csc(c + d*x)**6/(a*sec(c + d*x) + a)**3
    F = 4*cot(c + d*x)**11/(11*a**3*d) + 11*cot(c + d*x)**9/(9*a**3*d) + 10*cot(c + d*x)**7/(7*a**3*d) + 3*cot(c + d*x)**5/(5*a**3*d) - 4*csc(c + d*x)**11/(11*a**3*d) + 7*csc(c + d*x)**9/(9*a**3*d) - 3*csc(c + d*x)**7/(7*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_107():
    f = csc(c + d*x)**8/(a*sec(c + d*x) + a)**3
    F = 4*cot(c + d*x)**13/(13*a**3*d) + 15*cot(c + d*x)**11/(11*a**3*d) + 7*cot(c + d*x)**9/(3*a**3*d) + 13*cot(c + d*x)**7/(7*a**3*d) + 3*cot(c + d*x)**5/(5*a**3*d) - 4*csc(c + d*x)**13/(13*a**3*d) + 7*csc(c + d*x)**11/(11*a**3*d) - csc(c + d*x)**9/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_108():
    f = (e*sin(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)
    F = -a*e**(sympy.S(5)/2)*atan(sqrt(e*sin(c + d*x))/sqrt(e))/d + a*e**(sympy.S(5)/2)*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/d + 6*a*e**2*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*d*sqrt(sin(c + d*x))) - 2*a*e*(e*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*d) - 2*a*e*(e*sin(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_109():
    f = (e*sin(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)
    F = a*e**(sympy.S(3)/2)*atan(sqrt(e*sin(c + d*x))/sqrt(e))/d + a*e**(sympy.S(3)/2)*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/d + 2*a*e**2*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d*sqrt(e*sin(c + d*x))) - 2*a*e*sqrt(e*sin(c + d*x))*cos(c + d*x)/(3*d) - 2*a*e*sqrt(e*sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_110():
    f = sqrt(e*sin(c + d*x))*(a*sec(c + d*x) + a)
    F = -a*sqrt(e)*atan(sqrt(e*sin(c + d*x))/sqrt(e))/d + a*sqrt(e)*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/d + 2*a*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(d*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_111():
    f = (a*sec(c + d*x) + a)/sqrt(e*sin(c + d*x))
    F = 2*a*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(d*sqrt(e*sin(c + d*x))) + a*atan(sqrt(e*sin(c + d*x))/sqrt(e))/(d*sqrt(e)) + a*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/(d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_112():
    f = (a*sec(c + d*x) + a)/(e*sin(c + d*x))**(sympy.S(3)/2)
    F = -2*a*cos(c + d*x)/(d*e*sqrt(e*sin(c + d*x))) - 2*a/(d*e*sqrt(e*sin(c + d*x))) - 2*a*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(d*e**2*sqrt(sin(c + d*x))) - a*atan(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2)) + a*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_113():
    f = (a*sec(c + d*x) + a)/(e*sin(c + d*x))**(sympy.S(5)/2)
    F = -2*a*cos(c + d*x)/(3*d*e*(e*sin(c + d*x))**(sympy.S(3)/2)) - 2*a/(3*d*e*(e*sin(c + d*x))**(sympy.S(3)/2)) + 2*a*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d*e**2*sqrt(e*sin(c + d*x))) + a*atan(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(5)/2)) + a*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_114():
    f = (e*sin(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)**2
    F = -2*a**2*e**(sympy.S(5)/2)*atan(sqrt(e*sin(c + d*x))/sqrt(e))/d + 2*a**2*e**(sympy.S(5)/2)*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/d - 9*a**2*e**2*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*d*sqrt(sin(c + d*x))) - 2*a**2*e*(e*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*d) + a**2*e*(e*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)/d - 4*a**2*e*(e*sin(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_115():
    f = (e*sin(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)**2
    F = 2*a**2*e**(sympy.S(3)/2)*atan(sqrt(e*sin(c + d*x))/sqrt(e))/d + 2*a**2*e**(sympy.S(3)/2)*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/d - a**2*e**2*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d*sqrt(e*sin(c + d*x))) - 2*a**2*e*sqrt(e*sin(c + d*x))*cos(c + d*x)/(3*d) + a**2*e*sqrt(e*sin(c + d*x))*sec(c + d*x)/d - 4*a**2*e*sqrt(e*sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_116():
    f = sqrt(e*sin(c + d*x))*(a*sec(c + d*x) + a)**2
    F = -2*a**2*sqrt(e)*atan(sqrt(e*sin(c + d*x))/sqrt(e))/d + 2*a**2*sqrt(e)*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/d + a**2*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(d*sqrt(sin(c + d*x))) + a**2*(e*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)/(d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_117():
    f = (a*sec(c + d*x) + a)**2/sqrt(e*sin(c + d*x))
    F = 3*a**2*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(d*sqrt(e*sin(c + d*x))) + a**2*sqrt(e*sin(c + d*x))*sec(c + d*x)/(d*e) + 2*a**2*atan(sqrt(e*sin(c + d*x))/sqrt(e))/(d*sqrt(e)) + 2*a**2*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/(d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_118():
    f = (a*sec(c + d*x) + a)**2/(e*sin(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*cos(c + d*x)/(d*e*sqrt(e*sin(c + d*x))) - 2*a**2*sec(c + d*x)/(d*e*sqrt(e*sin(c + d*x))) - 4*a**2/(d*e*sqrt(e*sin(c + d*x))) - 5*a**2*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(d*e**2*sqrt(sin(c + d*x))) + 3*a**2*(e*sin(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)/(d*e**3) - 2*a**2*atan(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2)) + 2*a**2*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_119():
    f = (a*sec(c + d*x) + a)**2/(e*sin(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*cos(c + d*x)/(3*d*e*(e*sin(c + d*x))**(sympy.S(3)/2)) - 2*a**2*sec(c + d*x)/(3*d*e*(e*sin(c + d*x))**(sympy.S(3)/2)) - 4*a**2/(3*d*e*(e*sin(c + d*x))**(sympy.S(3)/2)) + 7*a**2*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d*e**2*sqrt(e*sin(c + d*x))) + 5*a**2*sqrt(e*sin(c + d*x))*sec(c + d*x)/(3*d*e**3) + 2*a**2*atan(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(5)/2)) + 2*a**2*atanh(sqrt(e*sin(c + d*x))/sqrt(e))/(d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_120():
    f = (e*sin(c + d*x))**(sympy.S(7)/2)/(a*sec(c + d*x) + a)
    F = -4*e**4*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a*d*sqrt(e*sin(c + d*x))) + 2*e**3*sqrt(e*sin(c + d*x))*cos(c + d*x)**3/(7*a*d) - 2*e**3*sqrt(e*sin(c + d*x))*cos(c + d*x)/(21*a*d) + 2*e*(e*sin(c + d*x))**(sympy.S(5)/2)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_121():
    f = (e*sin(c + d*x))**(sympy.S(5)/2)/(a*sec(c + d*x) + a)
    F = -4*e**2*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a*d*sqrt(sin(c + d*x))) - 2*e*(e*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*a*d) + 2*e*(e*sin(c + d*x))**(sympy.S(3)/2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_122():
    f = (e*sin(c + d*x))**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = -4*e**2*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*a*d*sqrt(e*sin(c + d*x))) - 2*e*sqrt(e*sin(c + d*x))*cos(c + d*x)/(3*a*d) + 2*e*sqrt(e*sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_123():
    f = sqrt(e*sin(c + d*x))/(a*sec(c + d*x) + a)
    F = 2*e*cos(c + d*x)/(a*d*sqrt(e*sin(c + d*x))) - 2*e/(a*d*sqrt(e*sin(c + d*x))) + 4*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(a*d*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_124():
    f = 1/(sqrt(e*sin(c + d*x))*(a*sec(c + d*x) + a))
    F = 2*e*cos(c + d*x)/(3*a*d*(e*sin(c + d*x))**(sympy.S(3)/2)) - 2*e/(3*a*d*(e*sin(c + d*x))**(sympy.S(3)/2)) + 4*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*a*d*sqrt(e*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_125():
    f = 1/((e*sin(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a))
    F = 2*e*cos(c + d*x)/(5*a*d*(e*sin(c + d*x))**(sympy.S(5)/2)) - 2*e/(5*a*d*(e*sin(c + d*x))**(sympy.S(5)/2)) - 4*cos(c + d*x)/(5*a*d*e*sqrt(e*sin(c + d*x))) - 4*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a*d*e**2*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_126():
    f = 1/((e*sin(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a))
    F = 2*e*cos(c + d*x)/(7*a*d*(e*sin(c + d*x))**(sympy.S(7)/2)) - 2*e/(7*a*d*(e*sin(c + d*x))**(sympy.S(7)/2)) - 4*cos(c + d*x)/(21*a*d*e*(e*sin(c + d*x))**(sympy.S(3)/2)) + 4*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a*d*e**2*sqrt(e*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_127():
    f = (e*sin(c + d*x))**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**2
    F = 52*e**4*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a**2*d*sqrt(e*sin(c + d*x))) + 2*e**3*sqrt(e*sin(c + d*x))*cos(c + d*x)**3/(7*a**2*d) + 26*e**3*sqrt(e*sin(c + d*x))*cos(c + d*x)/(21*a**2*d) - 4*e**3*sqrt(e*sin(c + d*x))/(a**2*d) + 4*e*(e*sin(c + d*x))**(sympy.S(5)/2)/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_128():
    f = (e*sin(c + d*x))**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**2
    F = -2*e**3*cos(c + d*x)**3/(a**2*d*sqrt(e*sin(c + d*x))) - 2*e**3*cos(c + d*x)/(a**2*d*sqrt(e*sin(c + d*x))) + 4*e**3/(a**2*d*sqrt(e*sin(c + d*x))) - 44*e**2*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a**2*d*sqrt(sin(c + d*x))) - 12*e*(e*sin(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*a**2*d) + 4*e*(e*sin(c + d*x))**(sympy.S(3)/2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_129():
    f = (e*sin(c + d*x))**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**2
    F = -2*e**3*cos(c + d*x)**3/(3*a**2*d*(e*sin(c + d*x))**(sympy.S(3)/2)) - 2*e**3*cos(c + d*x)/(3*a**2*d*(e*sin(c + d*x))**(sympy.S(3)/2)) + 4*e**3/(3*a**2*d*(e*sin(c + d*x))**(sympy.S(3)/2)) - 4*e**2*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(a**2*d*sqrt(e*sin(c + d*x))) - 4*e*sqrt(e*sin(c + d*x))*cos(c + d*x)/(3*a**2*d) + 4*e*sqrt(e*sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_130():
    f = sqrt(e*sin(c + d*x))/(a*sec(c + d*x) + a)**2
    F = -2*e**3*cos(c + d*x)**3/(5*a**2*d*(e*sin(c + d*x))**(sympy.S(5)/2)) - 2*e**3*cos(c + d*x)/(5*a**2*d*(e*sin(c + d*x))**(sympy.S(5)/2)) + 4*e**3/(5*a**2*d*(e*sin(c + d*x))**(sympy.S(5)/2)) + 16*e*cos(c + d*x)/(5*a**2*d*sqrt(e*sin(c + d*x))) - 4*e/(a**2*d*sqrt(e*sin(c + d*x))) + 28*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a**2*d*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_131():
    f = 1/(sqrt(e*sin(c + d*x))*(a*sec(c + d*x) + a)**2)
    F = -2*e**3*cos(c + d*x)**3/(7*a**2*d*(e*sin(c + d*x))**(sympy.S(7)/2)) - 2*e**3*cos(c + d*x)/(7*a**2*d*(e*sin(c + d*x))**(sympy.S(7)/2)) + 4*e**3/(7*a**2*d*(e*sin(c + d*x))**(sympy.S(7)/2)) + 16*e*cos(c + d*x)/(21*a**2*d*(e*sin(c + d*x))**(sympy.S(3)/2)) - 4*e/(3*a**2*d*(e*sin(c + d*x))**(sympy.S(3)/2)) + 20*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a**2*d*sqrt(e*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_132():
    f = 1/((e*sin(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)**2)
    F = -2*e**3*cos(c + d*x)**3/(9*a**2*d*(e*sin(c + d*x))**(sympy.S(9)/2)) - 2*e**3*cos(c + d*x)/(9*a**2*d*(e*sin(c + d*x))**(sympy.S(9)/2)) + 4*e**3/(9*a**2*d*(e*sin(c + d*x))**(sympy.S(9)/2)) + 16*e*cos(c + d*x)/(45*a**2*d*(e*sin(c + d*x))**(sympy.S(5)/2)) - 4*e/(5*a**2*d*(e*sin(c + d*x))**(sympy.S(5)/2)) - 4*cos(c + d*x)/(15*a**2*d*e*sqrt(e*sin(c + d*x))) - 4*sqrt(e*sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(15*a**2*d*e**2*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_133():
    f = 1/((e*sin(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)**2)
    F = -2*e**3*cos(c + d*x)**3/(11*a**2*d*(e*sin(c + d*x))**(sympy.S(11)/2)) - 2*e**3*cos(c + d*x)/(11*a**2*d*(e*sin(c + d*x))**(sympy.S(11)/2)) + 4*e**3/(11*a**2*d*(e*sin(c + d*x))**(sympy.S(11)/2)) + 16*e*cos(c + d*x)/(77*a**2*d*(e*sin(c + d*x))**(sympy.S(7)/2)) - 4*e/(7*a**2*d*(e*sin(c + d*x))**(sympy.S(7)/2)) - 4*cos(c + d*x)/(231*a**2*d*e*(e*sin(c + d*x))**(sympy.S(3)/2)) + 4*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(231*a**2*d*e**2*sqrt(e*sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_134():
    f = (e*sin(c + d*x))**m*(a*sec(c + d*x) + a)**3
    F = 3*a**3*(e*sin(c + d*x))**(m + 1)*sqrt(cos(c + d*x)**2)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*e*(m + 1)) + 3*a**3*(e*sin(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)) + a**3*(e*sin(c + d*x))**(m + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)) + a**3*(e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_135():
    f = (e*sin(c + d*x))**m*(a*sec(c + d*x) + a)**2
    F = a**2*(e*sin(c + d*x))**(m + 1)*sqrt(cos(c + d*x)**2)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*e*(m + 1)) + 2*a**2*(e*sin(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)) + a**2*(e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_136():
    f = (e*sin(c + d*x))**m*(a*sec(c + d*x) + a)
    F = a*(e*sin(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)) + a*(e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_137():
    f = (e*sin(c + d*x))**m/(a*sec(c + d*x) + a)
    F = -e*(e*sin(c + d*x))**(m - 1)/(a*d*(1 - m)) + e*(e*sin(c + d*x))**(m - 1)*cos(c + d*x)*hyper((sympy.S(-1)/2, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), sin(c + d*x)**2)/(a*d*(1 - m)*sqrt(cos(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_138():
    f = (e*sin(c + d*x))**m/(a*sec(c + d*x) + a)**2
    F = 2*e**3*(e*sin(c + d*x))**(m - 3)/(a**2*d*(3 - m)) - e**3*(e*sin(c + d*x))**(m - 3)*cos(c + d*x)*hyper((sympy.S(-3)/2, m/2 + sympy.S(-3)/2), (m/2 + sympy.S(-1)/2,), sin(c + d*x)**2)/(a**2*d*(3 - m)*sqrt(cos(c + d*x)**2)) - e**3*(e*sin(c + d*x))**(m - 3)*cos(c + d*x)*hyper((sympy.S(-1)/2, m/2 + sympy.S(-3)/2), (m/2 + sympy.S(-1)/2,), sin(c + d*x)**2)/(a**2*d*(3 - m)*sqrt(cos(c + d*x)**2)) - 2*e*(e*sin(c + d*x))**(m - 1)/(a**2*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_139():
    f = (e*sin(c + d*x))**m/(a*sec(c + d*x) + a)**3
    F = -4*e**5*(e*sin(c + d*x))**(m - 5)/(a**3*d*(5 - m)) + e**5*(e*sin(c + d*x))**(m - 5)*cos(c + d*x)*hyper((sympy.S(-5)/2, m/2 + sympy.S(-5)/2), (m/2 + sympy.S(-3)/2,), sin(c + d*x)**2)/(a**3*d*(5 - m)*sqrt(cos(c + d*x)**2)) + 3*e**5*(e*sin(c + d*x))**(m - 5)*cos(c + d*x)*hyper((sympy.S(-3)/2, m/2 + sympy.S(-5)/2), (m/2 + sympy.S(-3)/2,), sin(c + d*x)**2)/(a**3*d*(5 - m)*sqrt(cos(c + d*x)**2)) + 7*e**3*(e*sin(c + d*x))**(m - 3)/(a**3*d*(3 - m)) - 3*e*(e*sin(c + d*x))**(m - 1)/(a**3*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_140():
    f = (e*sin(c + d*x))**m*(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*a*e*(e*sin(c + d*x))**(m - 1)*(1 - cos(c + d*x))**(sympy.S.Half - m/2)*sqrt(a*sec(c + d*x) + a)*appellf1(sympy.S(-1)/2, sympy.S.Half - m/2, -m/2 - 1, sympy.S.Half, cos(c + d*x), -cos(c + d*x))/(d*(cos(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_141():
    f = (e*sin(c + d*x))**m*sqrt(a*sec(c + d*x) + a)
    F = -2*e*(e*sin(c + d*x))**(m - 1)*(1 - cos(c + d*x))**(sympy.S.Half - m/2)*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)*appellf1(sympy.S.Half, -m/2, sympy.S.Half - m/2, sympy.S(3)/2, -cos(c + d*x), cos(c + d*x))/(d*(cos(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_142():
    f = (e*sin(c + d*x))**m/sqrt(a*sec(c + d*x) + a)
    F = -2*e*(e*sin(c + d*x))**(m - 1)*(1 - cos(c + d*x))**(sympy.S.Half - m/2)*(cos(c + d*x) + 1)**(1 - m/2)*cos(c + d*x)*appellf1(sympy.S(3)/2, sympy.S.Half - m/2, 1 - m/2, sympy.S(5)/2, cos(c + d*x), -cos(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_143():
    f = (e*sin(c + d*x))**m/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*e*(e*sin(c + d*x))**(m - 1)*(1 - cos(c + d*x))**(sympy.S.Half - m/2)*(cos(c + d*x) + 1)**(1 - m/2)*cos(c + d*x)**2*appellf1(sympy.S(5)/2, sympy.S.Half - m/2, 2 - m/2, sympy.S(7)/2, cos(c + d*x), -cos(c + d*x))/(5*a*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_144():
    f = (e*sin(c + d*x))**m*(a*sec(c + d*x) + a)**n
    F = -e*(e*sin(c + d*x))**(m - 1)*(1 - cos(c + d*x))**(sympy.S.Half - m/2)*(a*sec(c + d*x) + a)**n*(cos(c + d*x) + 1)**(-m/2 - n + sympy.S.Half)*cos(c + d*x)*appellf1(1 - n, sympy.S.Half - m/2, -m/2 - n + sympy.S.Half, 2 - n, cos(c + d*x), -cos(c + d*x))/(d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_145():
    f = (a*sec(c + d*x) + a)**n*sin(c + d*x)**7
    F = -(1 - sec(c + d*x))**2*(a*sec(c + d*x) + a)**(n + 4)*cos(c + d*x)**7/(a**4*d*(1 - n)) - (3 - n)*(8 - n)*(16 - n)*(a*sec(c + d*x) + a)**(n + 4)*hyper((6, n + 4), (n + 5,), sec(c + d*x) + 1)/(42*a**4*d*(1 - n)*(n + 4)) + (a*sec(c + d*x) + a)**(n + 4)*(-6*n - (n**2 - 25*n + 108)*sec(c + d*x) + 48)*cos(c + d*x)**7/(42*a**4*d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_146():
    f = (a*sec(c + d*x) + a)**n*sin(c + d*x)**5
    F = (12 - n)*(a*sec(c + d*x) + a)**(n + 3)*cos(c + d*x)**4/(20*a**3*d) - (a*sec(c + d*x) + a)**(n + 3)*cos(c + d*x)**5/(5*a**3*d) + (a*sec(c + d*x) + a)**(n + 3)*(n**2 - 13*n + 32)*hyper((4, n + 3), (n + 4,), sec(c + d*x) + 1)/(20*a**3*d*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_147():
    f = (a*sec(c + d*x) + a)**n*sin(c + d*x)**3
    F = -(4 - n)*(a*sec(c + d*x) + a)**(n + 2)*hyper((3, n + 2), (n + 3,), sec(c + d*x) + 1)/(3*a**2*d*(n + 2)) + (a*sec(c + d*x) + a)**(n + 2)*cos(c + d*x)**3/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_148():
    f = (a*sec(c + d*x) + a)**n*sin(c + d*x)
    F = (a*sec(c + d*x) + a)**(n + 1)*hyper((2, n + 1), (n + 2,), sec(c + d*x) + 1)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_149():
    f = (a*sec(c + d*x) + a)**n*csc(c + d*x)
    F = -(a*sec(c + d*x) + a)**n*hyper((1, n), (n + 1,), sec(c + d*x)/2 + sympy.S.Half)/(2*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_150():
    f = (a*sec(c + d*x) + a)**n*csc(c + d*x)**3
    F = a*(a*sec(c + d*x) + a)**(n - 1)/(2*d*(1 - sec(c + d*x))) - a*(2 - n)*(a*sec(c + d*x) + a)**(n - 1)/(4*d*(1 - n)) - (n + 2)*(a*sec(c + d*x) + a)**n*hyper((1, n), (n + 1,), sec(c + d*x)/2 + sympy.S.Half)/(8*d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_151():
    f = (a*sec(c + d*x) + a)**n*csc(c + d*x)**5
    F = a**2*(a*sec(c + d*x) + a)**(n - 2)*(n**2 + 9*n + 12)*hyper((1, n - 2), (n - 1,), sec(c + d*x)/2 + sympy.S.Half)/(16*d*(2 - n)) - a**2*(a*sec(c + d*x) + a)**(n - 2)*(-n**3 - 7*n**2 + 4*n - (2 - 2*n)*(n + 6)*sec(c + d*x) + 12)/(8*d*(1 - sec(c + d*x))*(n**2 - 3*n + 2)) + a**2*(n + 3)*(a*sec(c + d*x) + a)**(n - 2)*sec(c + d*x)**2/(4*d*(1 - n)*(1 - sec(c + d*x))**2) - a**2*(a*sec(c + d*x) + a)**(n - 2)*sec(c + d*x)**3/(d*(1 - n)*(1 - sec(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_152():
    f = (a*sec(c + d*x) + a)**n*sin(c + d*x)**4
    F = 2**(n + sympy.S.Half)*(a*sec(c + d*x) + a)**n*(cos(c + d*x) + 1)**(-n + sympy.S(-1)/2)*sin(c + d*x)*cos(c + d*x)**n*appellf1(sympy.S.Half, sympy.S.Half - n, n - 4, sympy.S(3)/2, sympy.S.Half - cos(c + d*x)/2, 1 - cos(c + d*x))/d - (a*sec(c + d*x) + a)**n*sin(c + d*x)*cos(c + d*x)/d - (a*sec(c + d*x) + a)**n*(-n*cos(c + d*x) + n)*(cos(c + d*x) + 1)**(sympy.S.Half - n)*cot(c + d*x)*appellf1(1 - n, sympy.S(-1)/2, sympy.S.Half - n, 2 - n, cos(c + d*x), -cos(c + d*x))/(d*(1 - n)*sqrt(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_153():
    f = (a*sec(c + d*x) + a)**n*sin(c + d*x)**2
    F = -sqrt(1 - cos(c + d*x))*(a*sec(c + d*x) + a)**n*(cos(c + d*x) + 1)**(sympy.S.Half - n)*cot(c + d*x)*appellf1(1 - n, sympy.S(-1)/2, -n + sympy.S(-1)/2, 2 - n, cos(c + d*x), -cos(c + d*x))/(d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_154():
    f = (a*sec(c + d*x) + a)**n*csc(c + d*x)**2
    F = 2**(n + sympy.S(-1)/2)*n*(a*sec(c + d*x) + a)**n*(sec(c + d*x) + 1)**(-n + sympy.S(-1)/2)*tan(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/2 - n), (sympy.S(3)/2,), sympy.S.Half - sec(c + d*x)/2)/d - (a*sec(c + d*x) + a)**n*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_155():
    f = (a*sec(c + d*x) + a)**n*csc(c + d*x)**4
    F = -a**4*(a*sec(c + d*x) + a)**n*sin(c + d*x)*cos(c + d*x)/(d*(3 - 2*n)*(-a*cos(c + d*x) + a)**2*(a*cos(c + d*x) + a)**2) - a**3*(4 - n)*(a*sec(c + d*x) + a)**n*sin(c + d*x)*cos(c + d*x)/(d*(-a*cos(c + d*x) + a)**2*(a*cos(c + d*x) + a)*(4*n**2 - 8*n + 3)) + n*((cos(c + d*x) + 1)/(1 - cos(c + d*x)))**(-n + sympy.S(-1)/2)*(a*sec(c + d*x) + a)**n*(-n**2 - 3*n + 7)*sin(c + d*x)*cos(c + d*x)*hyper((1 - n, -n + sympy.S(-1)/2), (2 - n,), -2*cos(c + d*x)/(1 - cos(c + d*x)))/(d*(1 - 2*n)*(1 - n)*(1 - cos(c + d*x))**2*(3 - 2*n)*(2*n + 1)) + (a*sec(c + d*x) + a)**n*(n**2 - n + 2)*sin(c + d*x)*cos(c + d*x)/(d*(1 - 4*n**2)*(1 - cos(c + d*x))**2*(3 - 2*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_156():
    f = (a*sec(c + d*x) + a)**n*sin(c + d*x)**(sympy.S(3)/2)
    F = -(a*sec(c + d*x) + a)**n*(cos(c + d*x) + 1)**(-n + sympy.S(-1)/4)*sqrt(sin(c + d*x))*cos(c + d*x)*appellf1(1 - n, sympy.S(-1)/4, -n + sympy.S(-1)/4, 2 - n, cos(c + d*x), -cos(c + d*x))/(d*(1 - n)*(1 - cos(c + d*x))**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_157():
    f = (a*sec(c + d*x) + a)**n*sqrt(sin(c + d*x))
    F = -(1 - cos(c + d*x))**(sympy.S(1)/4)*(a*sec(c + d*x) + a)**n*(cos(c + d*x) + 1)**(sympy.S(1)/4 - n)*cos(c + d*x)*appellf1(1 - n, sympy.S(1)/4, sympy.S(1)/4 - n, 2 - n, cos(c + d*x), -cos(c + d*x))/(d*(1 - n)*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_158():
    f = (a*sec(c + d*x) + a)**n/sqrt(sin(c + d*x))
    F = -(1 - cos(c + d*x))**(sympy.S(3)/4)*(a*sec(c + d*x) + a)**n*(cos(c + d*x) + 1)**(sympy.S(3)/4 - n)*cos(c + d*x)*appellf1(1 - n, sympy.S(3)/4, sympy.S(3)/4 - n, 2 - n, cos(c + d*x), -cos(c + d*x))/(d*(1 - n)*sin(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_159():
    f = (a*sec(c + d*x) + a)**n/sin(c + d*x)**(sympy.S(3)/2)
    F = -(1 - cos(c + d*x))**(sympy.S(5)/4)*(a*sec(c + d*x) + a)**n*(cos(c + d*x) + 1)**(sympy.S(5)/4 - n)*cos(c + d*x)*appellf1(1 - n, sympy.S(5)/4, sympy.S(5)/4 - n, 2 - n, cos(c + d*x), -cos(c + d*x))/(d*(1 - n)*sin(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_160():
    f = (a + b*sec(c + d*x))*sin(c + d*x)**7
    F = a*cos(c + d*x)**7/(7*d) - 3*a*cos(c + d*x)**5/(5*d) + a*cos(c + d*x)**3/d - a*cos(c + d*x)/d - b*log(cos(c + d*x))/d + b*cos(c + d*x)**6/(6*d) - 3*b*cos(c + d*x)**4/(4*d) + 3*b*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_161():
    f = (a + b*sec(c + d*x))*sin(c + d*x)**5
    F = -a*cos(c + d*x)**5/(5*d) + 2*a*cos(c + d*x)**3/(3*d) - a*cos(c + d*x)/d - b*log(cos(c + d*x))/d - b*cos(c + d*x)**4/(4*d) + b*cos(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_162():
    f = (a + b*sec(c + d*x))*sin(c + d*x)**3
    F = a*cos(c + d*x)**3/(3*d) - a*cos(c + d*x)/d - b*log(cos(c + d*x))/d + b*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_163():
    f = (a + b*sec(c + d*x))*sin(c + d*x)
    F = -a*cos(c + d*x)/d - b*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_164():
    f = (a + b*sec(c + d*x))*csc(c + d*x)
    F = -a*atanh(cos(c + d*x))/d + b*log(tan(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_165():
    f = (a + b*sec(c + d*x))*csc(c + d*x)**3
    F = -a*cot(c + d*x)*csc(c + d*x)/(2*d) - a*atanh(cos(c + d*x))/(2*d) + b*log(tan(c + d*x))/d - b*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_166():
    f = (a + b*sec(c + d*x))*csc(c + d*x)**5
    F = -a*cot(c + d*x)*csc(c + d*x)**3/(4*d) - 3*a*cot(c + d*x)*csc(c + d*x)/(8*d) - 3*a*atanh(cos(c + d*x))/(8*d) + b*log(tan(c + d*x))/d - b*cot(c + d*x)**4/(4*d) - b*cot(c + d*x)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_167():
    f = (a + b*sec(c + d*x))*csc(c + d*x)**7
    F = -a*cot(c + d*x)*csc(c + d*x)**5/(6*d) - 5*a*cot(c + d*x)*csc(c + d*x)**3/(24*d) - 5*a*cot(c + d*x)*csc(c + d*x)/(16*d) - 5*a*atanh(cos(c + d*x))/(16*d) + b*log(tan(c + d*x))/d - b*cot(c + d*x)**6/(6*d) - 3*b*cot(c + d*x)**4/(4*d) - 3*b*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_168():
    f = (a + b*sec(c + d*x))*sin(c + d*x)**6
    F = 5*a*x/16 - a*sin(c + d*x)**5*cos(c + d*x)/(6*d) - 5*a*sin(c + d*x)**3*cos(c + d*x)/(24*d) - 5*a*sin(c + d*x)*cos(c + d*x)/(16*d) - b*sin(c + d*x)**5/(5*d) - b*sin(c + d*x)**3/(3*d) - b*sin(c + d*x)/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_169():
    f = (a + b*sec(c + d*x))*sin(c + d*x)**4
    F = 3*a*x/8 - a*sin(c + d*x)**3*cos(c + d*x)/(4*d) - 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - b*sin(c + d*x)**3/(3*d) - b*sin(c + d*x)/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_170():
    f = (a + b*sec(c + d*x))*sin(c + d*x)**2
    F = a*x/2 - a*sin(c + d*x)*cos(c + d*x)/(2*d) - b*sin(c + d*x)/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_171():
    f = (a + b*sec(c + d*x))*csc(c + d*x)**2
    F = -a*cot(c + d*x)/d + b*atanh(sin(c + d*x))/d - b*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_172():
    f = (a + b*sec(c + d*x))*csc(c + d*x)**4
    F = -a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + b*atanh(sin(c + d*x))/d - b*csc(c + d*x)**3/(3*d) - b*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_173():
    f = (a + b*sec(c + d*x))*csc(c + d*x)**6
    F = -a*cot(c + d*x)**5/(5*d) - 2*a*cot(c + d*x)**3/(3*d) - a*cot(c + d*x)/d + b*atanh(sin(c + d*x))/d - b*csc(c + d*x)**5/(5*d) - b*csc(c + d*x)**3/(3*d) - b*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_174():
    f = (a + b*sec(c + d*x))**2*sin(c + d*x)**5
    F = -a**2*cos(c + d*x)**5/(5*d) - 2*a*b*log(cos(c + d*x))/d - a*b*cos(c + d*x)**4/(2*d) + 2*a*b*cos(c + d*x)**2/d + b**2*sec(c + d*x)/d - (a**2 - 2*b**2)*cos(c + d*x)/d + (2*a**2 - b**2)*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_175():
    f = (a + b*sec(c + d*x))**2*sin(c + d*x)**3
    F = a**2*cos(c + d*x)**3/(3*d) - 2*a*b*log(cos(c + d*x))/d + a*b*cos(c + d*x)**2/d + b**2*sec(c + d*x)/d - (a**2 - b**2)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_176():
    f = (a + b*sec(c + d*x))**2*sin(c + d*x)
    F = -a**2*cos(c + d*x)/d - 2*a*b*log(cos(c + d*x))/d + b**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_177():
    f = (a + b*sec(c + d*x))**2*csc(c + d*x)
    F = -2*a*b*log(cos(c + d*x))/d + b**2*sec(c + d*x)/d - (a - b)**2*log(cos(c + d*x) + 1)/(2*d) + (a + b)**2*log(1 - cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_178():
    f = (a + b*sec(c + d*x))**2*csc(c + d*x)**3
    F = -2*a*b*log(cos(c + d*x))/d + b**2*sec(c + d*x)/d - (a - 3*b)*(a - b)*log(cos(c + d*x) + 1)/(4*d) + (a + b)*(a + 3*b)*log(1 - cos(c + d*x))/(4*d) - (2*a*b + (a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_179():
    f = (a + b*sec(c + d*x))**2*sin(c + d*x)**6
    F = -a**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) - 2*a*b*sin(c + d*x)**5/(5*d) - 2*a*b*sin(c + d*x)**3/(3*d) - 2*a*b*sin(c + d*x)/d + 2*a*b*atanh(sin(c + d*x))/d + b**2*tan(c + d*x)/d + x*(5*a**2/16 - 15*b**2/8) - (11*a**2 - 18*b**2)*sin(c + d*x)*cos(c + d*x)/(16*d) + (13*a**2 - 6*b**2)*sin(c + d*x)*cos(c + d*x)**3/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_180():
    f = (a + b*sec(c + d*x))**2*sin(c + d*x)**4
    F = 2*a*b*atanh(sin(c + d*x))/d + x*(3*a**2/8 - 3*b**2/2) - (39*a**2 + 2*b**2)*sin(c + d*x)*cos(c + d*x)/(24*d) + (a*cos(c + d*x) + b)**3*tan(c + d*x)/(b*d) - b*(28*a**2 + b**2)*sin(c + d*x)/(6*a*d) + (a*cos(c + d*x) + b)**3*sin(c + d*x)/(4*a*d) - (12*a**2 + b**2)*(a*cos(c + d*x) + b)**2*sin(c + d*x)/(12*a*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_181():
    f = (a + b*sec(c + d*x))**2*sin(c + d*x)**2
    F = a**2*x/2 - a**2*sin(c + d*x)*cos(c + d*x)/(2*d) - 2*a*b*sin(c + d*x)/d + 2*a*b*atanh(sin(c + d*x))/d - b**2*x + b**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_182():
    f = (a + b*sec(c + d*x))**2*csc(c + d*x)**2
    F = 2*a*b*atanh(sin(c + d*x))/d - 2*a*b*csc(c + d*x)/d + b**2*tan(c + d*x)/d - (a**2 + b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_183():
    f = (a + b*sec(c + d*x))**2*csc(c + d*x)**4
    F = 2*a*b*atanh(sin(c + d*x))/d - 2*a*b*csc(c + d*x)**3/(3*d) - 2*a*b*csc(c + d*x)/d + b**2*tan(c + d*x)/d - (a**2 + b**2)*cot(c + d*x)**3/(3*d) - (a**2 + 2*b**2)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_184():
    f = (a + b*sec(c + d*x))**2*csc(c + d*x)**6
    F = 2*a*b*atanh(sin(c + d*x))/d - 2*a*b*csc(c + d*x)**5/(5*d) - 2*a*b*csc(c + d*x)**3/(3*d) - 2*a*b*csc(c + d*x)/d + b**2*tan(c + d*x)/d - (a**2 + b**2)*cot(c + d*x)**5/(5*d) - (a**2 + 3*b**2)*cot(c + d*x)/d - (2*a**2 + 3*b**2)*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_185():
    f = (a + b*sec(c + d*x))**3*sin(c + d*x)**5
    F = -a**3*cos(c + d*x)**5/(5*d) - 3*a**2*b*cos(c + d*x)**4/(4*d) + 3*a*b**2*sec(c + d*x)/d - a*(a**2 - 6*b**2)*cos(c + d*x)/d + a*(2*a**2 - 3*b**2)*cos(c + d*x)**3/(3*d) + b**3*sec(c + d*x)**2/(2*d) - b*(3*a**2 - 2*b**2)*log(cos(c + d*x))/d + b*(6*a**2 - b**2)*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_186():
    f = (a + b*sec(c + d*x))**3*sin(c + d*x)**3
    F = a**3*cos(c + d*x)**3/(3*d) + 3*a**2*b*cos(c + d*x)**2/(2*d) + 3*a*b**2*sec(c + d*x)/d - a*(a**2 - 3*b**2)*cos(c + d*x)/d + b**3*sec(c + d*x)**2/(2*d) - b*(3*a**2 - b**2)*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_187():
    f = (a + b*sec(c + d*x))**3*sin(c + d*x)
    F = -a**3*cos(c + d*x)/d - 3*a**2*b*log(cos(c + d*x))/d + 3*a*b**2*sec(c + d*x)/d + b**3*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_188():
    f = (a + b*sec(c + d*x))**3*csc(c + d*x)
    F = 3*a*b**2*sec(c + d*x)/d + b**3*sec(c + d*x)**2/(2*d) - b*(3*a**2 + b**2)*log(cos(c + d*x))/d - (a - b)**3*log(cos(c + d*x) + 1)/(2*d) + (a + b)**3*log(1 - cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_189():
    f = (a + b*sec(c + d*x))**3*csc(c + d*x)**3
    F = -a**2*(a*(1 + 3*b**2/a**2)*cos(c + d*x) + b*(3 + b**2/a**2))*csc(c + d*x)**2/(2*d) + 3*a*b**2*sec(c + d*x)/d + b**3*sec(c + d*x)**2/(2*d) - b*(3*a**2 + 2*b**2)*log(cos(c + d*x))/d - (a - 4*b)*(a - b)**2*log(cos(c + d*x) + 1)/(4*d) + (a + b)**2*(a + 4*b)*log(1 - cos(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_190():
    f = (a + b*sec(c + d*x))**3*sin(c + d*x)**6
    F = 5*a**3*x/16 - a**3*sin(c + d*x)**5*cos(c + d*x)/(6*d) - 5*a**3*sin(c + d*x)**3*cos(c + d*x)/(24*d) - 5*a**3*sin(c + d*x)*cos(c + d*x)/(16*d) - 3*a**2*b*sin(c + d*x)**5/(5*d) - a**2*b*sin(c + d*x)**3/d - 3*a**2*b*sin(c + d*x)/d + 3*a**2*b*atanh(sin(c + d*x))/d - 45*a*b**2*x/8 - 3*a*b**2*sin(c + d*x)**4*tan(c + d*x)/(4*d) - 15*a*b**2*sin(c + d*x)**2*tan(c + d*x)/(8*d) + 45*a*b**2*tan(c + d*x)/(8*d) + b**3*sin(c + d*x)**3*tan(c + d*x)**2/(2*d) + 5*b**3*sin(c + d*x)**3/(6*d) + 5*b**3*sin(c + d*x)/(2*d) - 5*b**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_191():
    f = (a + b*sec(c + d*x))**3*sin(c + d*x)**4
    F = 3*a*x*(a**2 - 12*b**2)/8 - a*(21*a**2 - 2*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d) + a*(a*cos(c + d*x) + b)**4*tan(c + d*x)/(b**2*d) + 3*b*(2*a**2 - b**2)*atanh(sin(c + d*x))/(2*d) - b*(17*a**2 - b**2)*sin(c + d*x)/(2*d) - (6*a**2 - b**2)*(a*cos(c + d*x) + b)**2*sin(c + d*x)/(4*b*d) + (a*cos(c + d*x) + b)**4*tan(c + d*x)*sec(c + d*x)/(2*b*d) - (4*a**2 - b**2)*(a*cos(c + d*x) + b)**3*sin(c + d*x)/(4*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_192():
    f = (a + b*sec(c + d*x))**3*sin(c + d*x)**2
    F = -5*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) - 15*a**2*b*sin(c + d*x)/(2*d) + a*x*(a**2 - 6*b**2)/2 + 3*a*(a*cos(c + d*x) + b)**2*tan(c + d*x)/(2*d) + b*(6*a**2 - b**2)*atanh(sin(c + d*x))/(2*d) + (a*cos(c + d*x) + b)**3*tan(c + d*x)*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_193():
    f = (a + b*sec(c + d*x))**3*csc(c + d*x)**2
    F = -a**3*cot(c + d*x)/d + 3*a**2*b*atanh(sin(c + d*x))/d - 3*a**2*b*csc(c + d*x)/d + 3*a*b**2*tan(c + d*x)/d - 3*a*b**2*cot(c + d*x)/d + 3*b**3*atanh(sin(c + d*x))/(2*d) + b**3*csc(c + d*x)*sec(c + d*x)**2/(2*d) - 3*b**3*csc(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_194():
    f = (a + b*sec(c + d*x))**3*csc(c + d*x)**4
    F = -a**3*cot(c + d*x)**3/(3*d) - a**3*cot(c + d*x)/d + 3*a**2*b*atanh(sin(c + d*x))/d - a**2*b*csc(c + d*x)**3/d - 3*a**2*b*csc(c + d*x)/d + 3*a*b**2*tan(c + d*x)/d - a*b**2*cot(c + d*x)**3/d - 6*a*b**2*cot(c + d*x)/d + 5*b**3*atanh(sin(c + d*x))/(2*d) + b**3*csc(c + d*x)**3*sec(c + d*x)**2/(2*d) - 5*b**3*csc(c + d*x)**3/(6*d) - 5*b**3*csc(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_195():
    f = (a + b*sec(c + d*x))**3*csc(c + d*x)**6
    F = -a**3*cot(c + d*x)**5/(5*d) - 2*a**3*cot(c + d*x)**3/(3*d) - a**3*cot(c + d*x)/d + 3*a**2*b*atanh(sin(c + d*x))/d - 3*a**2*b*csc(c + d*x)**5/(5*d) - a**2*b*csc(c + d*x)**3/d - 3*a**2*b*csc(c + d*x)/d + 3*a*b**2*tan(c + d*x)/d - 3*a*b**2*cot(c + d*x)**5/(5*d) - 3*a*b**2*cot(c + d*x)**3/d - 9*a*b**2*cot(c + d*x)/d + 7*b**3*atanh(sin(c + d*x))/(2*d) + b**3*csc(c + d*x)**5*sec(c + d*x)**2/(2*d) - 7*b**3*csc(c + d*x)**5/(10*d) - 7*b**3*csc(c + d*x)**3/(6*d) - 7*b**3*csc(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_196():
    f = sin(c + d*x)**7/(a + b*sec(c + d*x))
    F = cos(c + d*x)**7/(7*a*d) - b*cos(c + d*x)**6/(6*a**2*d) - (3*a**2 - b**2)*cos(c + d*x)**5/(5*a**3*d) + b*(3*a**2 - b**2)*cos(c + d*x)**4/(4*a**4*d) + (3*a**4 - 3*a**2*b**2 + b**4)*cos(c + d*x)**3/(3*a**5*d) - b*(3*a**4 - 3*a**2*b**2 + b**4)*cos(c + d*x)**2/(2*a**6*d) - (a**2 - b**2)**3*cos(c + d*x)/(a**7*d) + b*(a**2 - b**2)**3*log(a*cos(c + d*x) + b)/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_197():
    f = sin(c + d*x)**5/(a + b*sec(c + d*x))
    F = -cos(c + d*x)**5/(5*a*d) + b*cos(c + d*x)**4/(4*a**2*d) + (2*a**2 - b**2)*cos(c + d*x)**3/(3*a**3*d) - b*(2*a**2 - b**2)*cos(c + d*x)**2/(2*a**4*d) - (a**2 - b**2)**2*cos(c + d*x)/(a**5*d) + b*(a**2 - b**2)**2*log(a*cos(c + d*x) + b)/(a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_198():
    f = sin(c + d*x)**3/(a + b*sec(c + d*x))
    F = cos(c + d*x)**3/(3*a*d) - b*cos(c + d*x)**2/(2*a**2*d) - (a**2 - b**2)*cos(c + d*x)/(a**3*d) + b*(a**2 - b**2)*log(a*cos(c + d*x) + b)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_199():
    f = sin(c + d*x)/(a + b*sec(c + d*x))
    F = -cos(c + d*x)/(a*d) + b*log(a*cos(c + d*x) + b)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_200():
    f = csc(c + d*x)/(a + b*sec(c + d*x))
    F = b*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)) + log(1 - cos(c + d*x))/(d*(2*a + 2*b)) - log(cos(c + d*x) + 1)/(d*(2*a - 2*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_201():
    f = csc(c + d*x)**3/(a + b*sec(c + d*x))
    F = a**2*b*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**2) + a*log(1 - cos(c + d*x))/(4*d*(a + b)**2) - a*log(cos(c + d*x) + 1)/(4*d*(a - b)**2) + (-a*cos(c + d*x) + b)*csc(c + d*x)**2/(d*(2*a**2 - 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_202():
    f = csc(c + d*x)**5/(a + b*sec(c + d*x))
    F = a**4*b*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**3) + a*(3*a + b)*log(1 - cos(c + d*x))/(16*d*(a + b)**3) - a*(3*a - b)*log(cos(c + d*x) + 1)/(16*d*(a - b)**3) + (-a*cos(c + d*x) + b)*csc(c + d*x)**4/(d*(4*a**2 - 4*b**2)) + (4*a**2*b - a*(3*a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**2/(8*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_203():
    f = sin(c + d*x)**6/(a + b*sec(c + d*x))
    F = (-5*a*cos(c + d*x) + 6*b)*sin(c + d*x)**5/(30*a**2*d) + (-a*(5*a**2 - 6*b**2)*cos(c + d*x) + 8*b*(a**2 - b**2))*sin(c + d*x)**3/(24*a**4*d) + (-a*(5*a**4 - 14*a**2*b**2 + 8*b**4)*cos(c + d*x) + 16*b*(a**2 - b**2)**2)*sin(c + d*x)/(16*a**6*d) - 2*b*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**7*d) + x*(5*a**6 - 30*a**4*b**2 + 40*a**2*b**4 - 16*b**6)/(16*a**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_204():
    f = sin(c + d*x)**4/(a + b*sec(c + d*x))
    F = (-3*a*cos(c + d*x) + 4*b)*sin(c + d*x)**3/(12*a**2*d) + (-a*(3*a**2 - 4*b**2)*cos(c + d*x) + 8*b*(a**2 - b**2))*sin(c + d*x)/(8*a**4*d) - 2*b*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d) + x*(3*a**4 - 12*a**2*b**2 + 8*b**4)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_205():
    f = sin(c + d*x)**2/(a + b*sec(c + d*x))
    F = (-a*cos(c + d*x) + 2*b)*sin(c + d*x)/(2*a**2*d) - 2*b*sqrt(a - b)*sqrt(a + b)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d) + x*(a**2 - 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_206():
    f = csc(c + d*x)**2/(a + b*sec(c + d*x))
    F = -2*a*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + (-a*cos(c + d*x) + b)*csc(c + d*x)/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_207():
    f = csc(c + d*x)**4/(a + b*sec(c + d*x))
    F = -2*a**3*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + (-a*cos(c + d*x) + b)*csc(c + d*x)**3/(d*(3*a**2 - 3*b**2)) + (3*a**2*b - a*(2*a**2 + b**2)*cos(c + d*x))*csc(c + d*x)/(3*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_208():
    f = csc(c + d*x)**6/(a + b*sec(c + d*x))
    F = -2*a**5*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (-a*cos(c + d*x) + b)*csc(c + d*x)**5/(d*(5*a**2 - 5*b**2)) + (5*a**2*b - a*(4*a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**3/(15*d*(a**2 - b**2)**2) + (15*a**4*b - a*(8*a**4 + 9*a**2*b**2 - 2*b**4)*cos(c + d*x))*csc(c + d*x)/(15*d*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_209():
    f = sin(c + d*x)**7/(a + b*sec(c + d*x))**2
    F = cos(c + d*x)**7/(7*a**2*d) - b*cos(c + d*x)**6/(3*a**3*d) - (3*a**2 - 3*b**2)*cos(c + d*x)**5/(5*a**4*d) + b*(3*a**2 - 2*b**2)*cos(c + d*x)**4/(2*a**5*d) + (3*a**4 - 9*a**2*b**2 + 5*b**4)*cos(c + d*x)**3/(3*a**6*d) - 3*b*(a**2 - b**2)**2*cos(c + d*x)**2/(a**7*d) - (a**2 - 7*b**2)*(a**2 - b**2)**2*cos(c + d*x)/(a**8*d) + b**2*(a**2 - b**2)**3/(a**9*d*(a*cos(c + d*x) + b)) + 2*b*(a**2 - 4*b**2)*(a**2 - b**2)**2*log(a*cos(c + d*x) + b)/(a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_210():
    f = sin(c + d*x)**5/(a + b*sec(c + d*x))**2
    F = -cos(c + d*x)**5/(5*a**2*d) + b*cos(c + d*x)**4/(2*a**3*d) + (2*a**2 - 3*b**2)*cos(c + d*x)**3/(3*a**4*d) - 2*b*(a**2 - b**2)*cos(c + d*x)**2/(a**5*d) - (a**4 - 6*a**2*b**2 + 5*b**4)*cos(c + d*x)/(a**6*d) + b**2*(a**2 - b**2)**2/(a**7*d*(a*cos(c + d*x) + b)) + 2*b*(a**4 - 4*a**2*b**2 + 3*b**4)*log(a*cos(c + d*x) + b)/(a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_211():
    f = sin(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = cos(c + d*x)**3/(3*a**2*d) - b*cos(c + d*x)**2/(a**3*d) - (a**2 - 3*b**2)*cos(c + d*x)/(a**4*d) + b**2*(a**2 - b**2)/(a**5*d*(a*cos(c + d*x) + b)) + 2*b*(a**2 - 2*b**2)*log(a*cos(c + d*x) + b)/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_212():
    f = sin(c + d*x)/(a + b*sec(c + d*x))**2
    F = -cos(c + d*x)/(a**2*d) + b**2/(a**3*d*(a*cos(c + d*x) + b)) + 2*b*log(a*cos(c + d*x) + b)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_213():
    f = csc(c + d*x)/(a + b*sec(c + d*x))**2
    F = 2*a*b*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**2) + log(1 - cos(c + d*x))/(2*d*(a + b)**2) - log(cos(c + d*x) + 1)/(2*d*(a - b)**2) + b**2/(a*d*(a**2 - b**2)*(a*cos(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_214():
    f = csc(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = a*b**2/(d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) + 2*a*b*(a**2 + b**2)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**3) + (a - b)*log(1 - cos(c + d*x))/(4*d*(a + b)**3) + (2*a*b - (a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**2/(2*d*(a**2 - b**2)**2) - (a + b)*log(cos(c + d*x) + 1)/(4*d*(a - b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_215():
    f = csc(c + d*x)**5/(a + b*sec(c + d*x))**2
    F = a**3*b**2/(d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) + 2*a**3*b*(a**2 + 2*b**2)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) + (2*a*b - (a**2 + b**2)*cos(c + d*x))*csc(c + d*x)**4/(4*d*(a**2 - b**2)**2) + (8*a*b*(a**2 + b**2) - (3*a**4 + 12*a**2*b**2 + b**4)*cos(c + d*x))*csc(c + d*x)**2/(8*d*(a**2 - b**2)**3) + (3*a**2 - 4*a*b - b**2)*log(1 - cos(c + d*x))/(16*d*(a + b)**4) - (3*a**2 + 4*a*b - b**2)*log(cos(c + d*x) + 1)/(16*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_216():
    f = sin(c + d*x)**6/(a + b*sec(c + d*x))**2
    F = a*sin(c + d*x)*cos(c + d*x)**4/(6*b**2*d*(a*cos(c + d*x) + b)) - sin(c + d*x)*cos(c + d*x)**3/(3*b*d*(a*cos(c + d*x) + b)) - sin(c + d*x)*cos(c + d*x)**6/(6*a*d*(a*cos(c + d*x) + b)) + 7*b*sin(c + d*x)*cos(c + d*x)**5/(30*a**2*d*(a*cos(c + d*x) + b)) + (5*a**4 - 20*a**2*b**2 + 14*b**4)*sin(c + d*x)*cos(c + d*x)**4/(10*a**3*b**2*d*(a*cos(c + d*x) + b)) - (16*a**4 - 61*a**2*b**2 + 42*b**4)*sin(c + d*x)*cos(c + d*x)**3/(24*a**4*b**2*d) + (15*a**4 - 52*a**2*b**2 + 35*b**4)*sin(c + d*x)*cos(c + d*x)**2/(15*a**5*b*d) - (27*a**4 - 86*a**2*b**2 + 56*b**4)*sin(c + d*x)*cos(c + d*x)/(16*a**6*d) + b*(61*a**4 - 170*a**2*b**2 + 105*b**4)*sin(c + d*x)/(15*a**7*d) - 2*b*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*(2*a**2 - 7*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**8*d) + x*(5*a**6 - 90*a**4*b**2 + 200*a**2*b**4 - 112*b**6)/(16*a**8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_217():
    f = sin(c + d*x)**4/(a + b*sec(c + d*x))**2
    F = sin(c + d*x)*cos(c + d*x)**3/(4*a**2*d) - (a**2 - b**2)*sin(c + d*x)*cos(c + d*x)**3/(a**2*b*d*(a*cos(c + d*x) + b)) + (3*a**2 - 5*b**2)*sin(c + d*x)*cos(c + d*x)**2/(3*a**3*b*d) - (13*a**2 - 20*b**2)*sin(c + d*x)*cos(c + d*x)/(8*a**4*d) + b*(11*a**2 - 15*b**2)*sin(c + d*x)/(3*a**5*d) - 2*b*sqrt(a - b)*sqrt(a + b)*(2*a**2 - 5*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**6*d) + x*(3*a**4 - 36*a**2*b**2 + 40*b**4)/(8*a**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_218():
    f = sin(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = sin(c + d*x)*cos(c + d*x)**2/(a*d*(a*cos(c + d*x) + b)) - 3*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) + 3*b*sin(c + d*x)/(a**3*d) - 2*b*(2*a**2 - 3*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*sqrt(a - b)*sqrt(a + b)) + x*(a**2 - 6*b**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_219():
    f = csc(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = -4*a**2*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a*b**2*sin(c + d*x)/(d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) - 2*b**3*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + sin(c + d*x)/(2*d*(a - b)**2*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_220():
    f = csc(c + d*x)**4/(a + b*sec(c + d*x))**2
    F = a**3*b**2*sin(c + d*x)/(d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) - 2*a**2*b**3*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - 4*a**2*b*(a**2 + b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + sin(c + d*x)/(12*d*(a - b)**2*(cos(c + d*x) + 1)) + sin(c + d*x)/(12*d*(a - b)**2*(cos(c + d*x) + 1)**2) + (a + b)*sin(c + d*x)/(4*d*(a - b)**3*(cos(c + d*x) + 1)) - (a - b)*sin(c + d*x)/(4*d*(1 - cos(c + d*x))*(a + b)**3) - sin(c + d*x)/(12*d*(1 - cos(c + d*x))*(a + b)**2) - sin(c + d*x)/(12*d*(1 - cos(c + d*x))**2*(a + b)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_221():
    f = sin(c + d*x)**7/(a + b*sec(c + d*x))**3
    F = cos(c + d*x)**7/(7*a**3*d) - b*cos(c + d*x)**6/(2*a**4*d) - (3*a**2 - 6*b**2)*cos(c + d*x)**5/(5*a**5*d) + b*(9*a**2 - 10*b**2)*cos(c + d*x)**4/(4*a**6*d) + (a**4 - 6*a**2*b**2 + 5*b**4)*cos(c + d*x)**3/(a**7*d) - 3*b*(3*a**4 - 10*a**2*b**2 + 7*b**4)*cos(c + d*x)**2/(2*a**8*d) - (a**6 - 18*a**4*b**2 + 45*a**2*b**4 - 28*b**6)*cos(c + d*x)/(a**9*d) - b**3*(a**2 - b**2)**3/(2*a**10*d*(a*cos(c + d*x) + b)**2) + 3*b**2*(a**2 - 3*b**2)*(a**2 - b**2)**2/(a**10*d*(a*cos(c + d*x) + b)) + 3*b*(a**2 - b**2)*(a**4 - 9*a**2*b**2 + 12*b**4)*log(a*cos(c + d*x) + b)/(a**10*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_222():
    f = sin(c + d*x)**5/(a + b*sec(c + d*x))**3
    F = -cos(c + d*x)**5/(5*a**3*d) + 3*b*cos(c + d*x)**4/(4*a**4*d) + (2*a**2 - 6*b**2)*cos(c + d*x)**3/(3*a**5*d) - b*(3*a**2 - 5*b**2)*cos(c + d*x)**2/(a**6*d) - (a**4 - 12*a**2*b**2 + 15*b**4)*cos(c + d*x)/(a**7*d) - b**3*(a**2 - b**2)**2/(2*a**8*d*(a*cos(c + d*x) + b)**2) + b**2*(3*a**4 - 10*a**2*b**2 + 7*b**4)/(a**8*d*(a*cos(c + d*x) + b)) + b*(3*a**4 - 20*a**2*b**2 + 21*b**4)*log(a*cos(c + d*x) + b)/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_223():
    f = sin(c + d*x)**3/(a + b*sec(c + d*x))**3
    F = cos(c + d*x)**3/(3*a**3*d) - 3*b*cos(c + d*x)**2/(2*a**4*d) - (a**2 - 6*b**2)*cos(c + d*x)/(a**5*d) - b**3*(a**2 - b**2)/(2*a**6*d*(a*cos(c + d*x) + b)**2) + b**2*(3*a**2 - 5*b**2)/(a**6*d*(a*cos(c + d*x) + b)) + b*(3*a**2 - 10*b**2)*log(a*cos(c + d*x) + b)/(a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_224():
    f = sin(c + d*x)/(a + b*sec(c + d*x))**3
    F = -cos(c + d*x)/(a**3*d) - b**3/(2*a**4*d*(a*cos(c + d*x) + b)**2) + 3*b**2/(a**4*d*(a*cos(c + d*x) + b)) + 3*b*log(a*cos(c + d*x) + b)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_225():
    f = csc(c + d*x)/(a + b*sec(c + d*x))**3
    F = b*(3*a**2 + b**2)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**3) + log(1 - cos(c + d*x))/(2*d*(a + b)**3) - log(cos(c + d*x) + 1)/(2*d*(a - b)**3) - b**3/(2*a**2*d*(a**2 - b**2)*(a*cos(c + d*x) + b)**2) + b**2*(3*a**2 - b**2)/(a**2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_226():
    f = csc(c + d*x)**3/(a + b*sec(c + d*x))**3
    F = -b**3/(2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) + b**2*(3*a**2 + b**2)/(d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) + b*(3*a**4 + 8*a**2*b**2 + b**4)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**4) + (a - 2*b)*log(1 - cos(c + d*x))/(4*d*(a + b)**4) + (-a*(a**2 + 3*b**2)*cos(c + d*x) + b*(3*a**2 + b**2))*csc(c + d*x)**2/(2*d*(a**2 - b**2)**3) - (a + 2*b)*log(cos(c + d*x) + 1)/(4*d*(a - b)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_227():
    f = csc(c + d*x)**5/(a + b*sec(c + d*x))**3
    F = -a**2*b**3/(2*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)**2) + 3*a**2*b**2*(a**2 + b**2)/(d*(a**2 - b**2)**4*(a*cos(c + d*x) + b)) + 3*a**2*b*(a**4 + 5*a**2*b**2 + 2*b**4)*log(a*cos(c + d*x) + b)/(d*(a**2 - b**2)**5) + 3*a*(a - 3*b)*log(1 - cos(c + d*x))/(16*d*(a + b)**5) - 3*a*(a + 3*b)*log(cos(c + d*x) + 1)/(16*d*(a - b)**5) + (-a*(a**2 + 3*b**2)*cos(c + d*x) + b*(3*a**2 + b**2))*csc(c + d*x)**4/(4*d*(a**2 - b**2)**3) + (-3*a*(a**4 + 10*a**2*b**2 + 5*b**4)*cos(c + d*x) + 4*b*(3*a**4 + 8*a**2*b**2 + b**4))*csc(c + d*x)**2/(8*d*(a**2 - b**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_228():
    f = sin(c + d*x)**6/(a + b*sec(c + d*x))**3
    F = a*sin(c + d*x)*cos(c + d*x)**5/(10*b**2*d*(a*cos(c + d*x) + b)**2) - sin(c + d*x)*cos(c + d*x)**4/(4*b*d*(a*cos(c + d*x) + b)**2) - sin(c + d*x)*cos(c + d*x)**7/(6*a*d*(a*cos(c + d*x) + b)**2) + 4*b*sin(c + d*x)*cos(c + d*x)**6/(15*a**2*d*(a*cos(c + d*x) + b)**2) + (9*a**4 - 60*a**2*b**2 + 56*b**4)*sin(c + d*x)*cos(c + d*x)**5/(60*a**3*b**2*d*(a*cos(c + d*x) + b)**2) + (15*a**4 - 110*a**2*b**2 + 112*b**4)*sin(c + d*x)*cos(c + d*x)**4/(20*a**4*b**2*d*(a*cos(c + d*x) + b)) - (24*a**4 - 169*a**2*b**2 + 168*b**4)*sin(c + d*x)*cos(c + d*x)**3/(24*a**5*b**2*d) + (45*a**4 - 291*a**2*b**2 + 280*b**4)*sin(c + d*x)*cos(c + d*x)**2/(30*a**6*b*d) - (43*a**4 - 244*a**2*b**2 + 224*b**4)*sin(c + d*x)*cos(c + d*x)/(16*a**7*d) + b*(213*a**4 - 985*a**2*b**2 + 840*b**4)*sin(c + d*x)/(30*a**8*d) - b*sqrt(a - b)*sqrt(a + b)*(6*a**4 - 47*a**2*b**2 + 56*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**9*d) + x*(5*a**6 - 180*a**4*b**2 + 600*a**2*b**4 - 448*b**6)/(16*a**9)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_229():
    f = sin(c + d*x)**4/(a + b*sec(c + d*x))**3
    F = -(a**2 - b**2)*sin(c + d*x)*cos(c + d*x)**4/(2*a**2*b*d*(a*cos(c + d*x) + b)**2) + (2*a**2 - 7*b**2)*sin(c + d*x)*cos(c + d*x)**4/(2*a**2*b**2*d*(a*cos(c + d*x) + b)) - (4*a**2 - 15*b**2)*sin(c + d*x)*cos(c + d*x)**3/(4*a**3*b**2*d) + (3*a**2 - 10*b**2)*sin(c + d*x)*cos(c + d*x)**2/(2*a**4*b*d) - (21*a**2 - 60*b**2)*sin(c + d*x)*cos(c + d*x)/(8*a**5*d) + b*(13*a**2 - 30*b**2)*sin(c + d*x)/(2*a**6*d) - 3*b*(2*a**4 - 11*a**2*b**2 + 10*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**7*d*sqrt(a - b)*sqrt(a + b)) + x*(3*a**4 - 72*a**2*b**2 + 120*b**4)/(8*a**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_230():
    f = sin(c + d*x)**2/(a + b*sec(c + d*x))**3
    F = sin(c + d*x)*cos(c + d*x)**3/(2*a*d*(a*cos(c + d*x) + b)**2) + (3*a**2 - 4*b**2)*sin(c + d*x)*cos(c + d*x)**2/(2*a**2*d*(a**2 - b**2)*(a*cos(c + d*x) + b)) - (5*a**2 - 6*b**2)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d*(a**2 - b**2)) + b*(11*a**2 - 12*b**2)*sin(c + d*x)/(2*a**4*d*(a**2 - b**2)) - b*(6*a**4 - 19*a**2*b**2 + 12*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x*(a**2 - 12*b**2)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_231():
    f = csc(c + d*x)**2/(a + b*sec(c + d*x))**3
    F = -2*a*b*(3*a**2 + b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + 3*b**4*sin(c + d*x)/(2*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) - b**3*sin(c + d*x)/(2*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)**2) + b**2*(3*a**2 - b**2)*sin(c + d*x)/(d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) + sin(c + d*x)/(2*d*(a - b)**3*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**3) - b**3*(a**2 + 2*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - 2*b**3*(3*a**2 - b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_232():
    f = csc(c + d*x)**4/(a + b*sec(c + d*x))**3
    F = 3*a**2*b**4*sin(c + d*x)/(2*d*(a**2 - b**2)**4*(a*cos(c + d*x) + b)) - a**2*b**3*sin(c + d*x)/(2*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)**2) + a**2*b**2*(3*a**2 + b**2)*sin(c + d*x)/(d*(a**2 - b**2)**4*(a*cos(c + d*x) + b)) - a*b**3*(a**2 + 2*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(9)/2)*(a + b)**(sympy.S(9)/2)) - 2*a*b**3*(3*a**2 + b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(9)/2)*(a + b)**(sympy.S(9)/2)) - 2*a*b*(3*a**4 + 8*a**2*b**2 + b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(9)/2)*(a + b)**(sympy.S(9)/2)) + sin(c + d*x)/(12*d*(a - b)**3*(cos(c + d*x) + 1)) + sin(c + d*x)/(12*d*(a - b)**3*(cos(c + d*x) + 1)**2) + (a + 2*b)*sin(c + d*x)/(4*d*(a - b)**4*(cos(c + d*x) + 1)) - (a - 2*b)*sin(c + d*x)/(4*d*(1 - cos(c + d*x))*(a + b)**4) - sin(c + d*x)/(12*d*(1 - cos(c + d*x))*(a + b)**3) - sin(c + d*x)/(12*d*(1 - cos(c + d*x))**2*(a + b)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_233():
    f = (e*sin(c + d*x))**(sympy.S(7)/2)/(a + b*sec(c + d*x))
    F = (Integer(-1) * ((Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Integer(5) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(28) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(21) * (Symbol('b'))**(Integer(4)))) * (Symbol('e'))**(Integer(4)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(21) * (Symbol('a'))**(Integer(5)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(5)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(5)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('e'))**(Integer(3)) * ((Integer(21) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) + (Integer(-1) * (Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(21) * (Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('e') * ((Integer(7) * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(35) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_234():
    f = (e*sin(c + d*x))**(sympy.S(5)/2)/(a + b*sec(c + d*x))
    F = ((Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * Symbol('e') * ((Integer(5) * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_235():
    f = (e*sin(c + d*x))**(sympy.S(3)/2)/(a + b*sec(c + d*x))
    F = (Integer(-1) * ((Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(2) * Symbol('e') * ((Integer(3) * Symbol('b')) + (Integer(-1) * (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_236():
    f = sqrt(e*sin(c + d*x))/(a + b*sec(c + d*x))
    F = ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('e') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * Symbol('e') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_237():
    f = 1/(sqrt(e*sin(c + d*x))*(a + b*sec(c + d*x)))
    F = (Integer(-1) * ((Symbol('b') * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(2) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_238():
    f = 1/((e*sin(c + d*x))**(sympy.S(3)/2)*(a + b*sec(c + d*x)))
    F = ((sympy.sqrt(Symbol('a')) * Symbol('b') * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('a')) * Symbol('b') * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b') + (Integer(-1) * (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_239():
    f = 1/((e*sin(c + d*x))**(sympy.S(5)/2)*(a + b*sec(c + d*x)))
    F = (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('b') * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('b') * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b') + (Integer(-1) * (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_240():
    f = 1/((e*sin(c + d*x))**(sympy.S(7)/2)*(a + b*sec(c + d*x)))
    F = (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('b') * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('b') * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b') + (Integer(-1) * (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(5) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b')) + (Integer(-1) * (Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * (Symbol('e'))**(Integer(3)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(4)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_241():
    f = (e*sin(c + d*x))**(sympy.S(9)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((Integer(7) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(7) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(13) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(7) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(5)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(7)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(5)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(7)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(7) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(5)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(7)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(5)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(7)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(14) * (Symbol('e'))**(Integer(4)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(6)) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(6)) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(14) * (Symbol('e'))**(Integer(3)) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(45) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * ((Integer(5) * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * (Symbol('e'))**(Integer(3)) * ((Integer(5) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) + (Integer(3) * Symbol('a') * Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('d')))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(7) * (Integer(2))**(Integer(-1))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_242():
    f = (e*sin(c + d*x))**(sympy.S(7)/2)/(a + b*sec(c + d*x))**2
    F = ((Integer(5) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(5) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(11) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(10) * (Symbol('e'))**(Integer(4)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(21) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(6)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(4)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(6)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('e'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(6)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(10) * (Symbol('e'))**(Integer(3)) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(21) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(3)) * ((Integer(3) * Symbol('b')) + (Integer(-1) * (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * (Symbol('e'))**(Integer(3)) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) + (Symbol('a') * Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('d')))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(5) * (Integer(2))**(Integer(-1))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_243():
    f = (e*sin(c + d*x))**(sympy.S(5)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(3)) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(9) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(4)) * (Symbol('e'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(5)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(5)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('b'))**(Integer(4)) * (Symbol('e'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(5)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(5)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(6) * (Symbol('e'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(4)) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(4) * Symbol('b') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_244():
    f = (e*sin(c + d*x))**(sympy.S(3)/2)/(a + b*sec(c + d*x))**2
    F = (((Symbol('b'))**(Integer(3)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(3)) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('e'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(2)) * (Symbol('e'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_245():
    f = sqrt(e*sin(c + d*x))/(a + b*sec(c + d*x))**2
    F = (((Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * sympy.sqrt(Symbol('e')) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt(Symbol('e')) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('e') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * Symbol('e') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('e') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(4)) * Symbol('e') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_246():
    f = 1/(sqrt(e*sin(c + d*x))*(a + b*sec(c + d*x))**2)
    F = (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((Symbol('a'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + ((Integer(2) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_247():
    f = 1/((e*sin(c + d*x))**(sympy.S(3)/2)*(a + b*sec(c + d*x))**2)
    F = ((Integer(5) * (Symbol('b'))**(Integer(3)) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('a')) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((sympy.sqrt(Symbol('a')) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(5) * Symbol('a') * Symbol('b')) + (Integer(-1) * (((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_248():
    f = 1/((e*sin(c + d*x))**(sympy.S(5)/2)*(a + b*sec(c + d*x))**2)
    F = (Integer(-1) * ((Integer(7) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**(Integer(3)) * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('a')) * Symbol('b') * sympy.atan(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * sympy.sqrt(Symbol('a')) * (Symbol('b'))**(Integer(3)) * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Symbol('a')) * Symbol('b') * sympy.atanh(((sympy.sqrt(Symbol('a')) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(4))**(Integer(-1))) * sympy.sqrt(Symbol('e'))))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b'))**(Integer(2)) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(4) * Symbol('b') * (Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(7) * Symbol('a') * Symbol('b')) + (Integer(-1) * (((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(4) * (Symbol('b'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(7) * (Symbol('b'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(7) * (Symbol('b'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(2))**(Integer(-1)))) + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))) + (Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_249():
    f = sqrt(a + b*sec(e + f*x))
    F = Integer(-1) * ((Integer(2) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('a') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_250():
    f = sqrt(a + b*sec(e + f*x))*csc(e + f*x)**2
    F = sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*sqrt(a + b)*cot(e + f*x)*elliptic_f(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/f - sqrt(a + b*sec(e + f*x))*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_251():
    f = (a + b*sec(e + f*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1)))) + ((Integer(2) * ((Integer(2) * Symbol('a')) + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_252():
    f = (a + b*sec(e + f*x))**(sympy.S(3)/2)*csc(e + f*x)**2
    F = -sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*sqrt(a + b)*(3*a - 3*b)*cot(e + f*x)*elliptic_e(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/f + sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*sqrt(a + b)*(3*a - 3*b)*cot(e + f*x)*elliptic_f(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/f - (a + b*sec(e + f*x))**(sympy.S(3)/2)*cot(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_253():
    f = 1/sqrt(a + b*sec(e + f*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_254():
    f = csc(e + f*x)**2/sqrt(a + b*sec(e + f*x))
    F = b**2*tan(e + f*x)/(f*sqrt(a + b*sec(e + f*x))*(a**2 - b**2)) + sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*cot(e + f*x)*elliptic_e(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/(f*sqrt(a + b)) - sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*cot(e + f*x)*elliptic_f(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/(f*sqrt(a + b)) - cot(e + f*x)/(f*sqrt(a + b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_255():
    f = (a + b*sec(e + f*x))**(sympy.S(-3)/2)
    F = ((Integer(2) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_256():
    f = csc(e + f*x)**2/(a + b*sec(e + f*x))**(sympy.S(3)/2)
    F = 4*a*b**2*tan(e + f*x)/(f*sqrt(a + b*sec(e + f*x))*(a**2 - b**2)**2) + 4*a*sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*cot(e + f*x)*elliptic_e(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/(f*(a - b)*(a + b)**(sympy.S(3)/2)) + b**2*tan(e + f*x)/(f*(a + b*sec(e + f*x))**(sympy.S(3)/2)*(a**2 - b**2)) - sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*(3*a - b)*cot(e + f*x)*elliptic_f(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/(f*(a - b)*(a + b)**(sympy.S(3)/2)) - cot(e + f*x)/(f*(a + b*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_257():
    f = (e*sin(c + d*x))**m*(a + b*sec(c + d*x))**3
    F = a**3*(e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)*sqrt(cos(c + d*x)**2)) + 3*a**2*b*(e*sin(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)) + 3*a*b**2*(e*sin(c + d*x))**(m + 1)*sqrt(cos(c + d*x)**2)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*e*(m + 1)) + b**3*(e*sin(c + d*x))**(m + 1)*hyper((2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_258():
    f = (e*sin(c + d*x))**m*(a + b*sec(c + d*x))**2
    F = a**2*(e*sin(c + d*x))**m*sin(c + d*x)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 1)*sqrt(cos(c + d*x)**2)) + 2*a*b*(e*sin(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)) + b**2*(e*sin(c + d*x))**m*sqrt(cos(c + d*x)**2)*tan(c + d*x)*hyper((sympy.S(3)/2, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_259():
    f = (e*sin(c + d*x))**m*(a + b*sec(c + d*x))
    F = a*(e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1)*sqrt(cos(c + d*x)**2)) + b*(e*sin(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_260():
    f = (e*sin(c + d*x))**m/(a + b*sec(c + d*x))
    F = (e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(a*d*e*(m + 1)*sqrt(cos(c + d*x)**2)) - b*e*(e*sin(c + d*x))**(m - 1)*(a*(cos(c + d*x) + 1)/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*(-a*(1 - cos(c + d*x))/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*appellf1(1 - m, sympy.S.Half - m/2, sympy.S.Half - m/2, 2 - m, -(a - b)/(a*cos(c + d*x) + b), (a + b)/(a*cos(c + d*x) + b))/(a**2*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_261():
    f = (e*sin(c + d*x))**m/(a + b*sec(c + d*x))**2
    F = (e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(a**2*d*e*(m + 1)*sqrt(cos(c + d*x)**2)) + b**2*e*(e*sin(c + d*x))**(m - 1)*(a*(cos(c + d*x) + 1)/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*(-a*(1 - cos(c + d*x))/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*appellf1(2 - m, sympy.S.Half - m/2, sympy.S.Half - m/2, 3 - m, -(a - b)/(a*cos(c + d*x) + b), (a + b)/(a*cos(c + d*x) + b))/(a**3*d*(2 - m)*(a*cos(c + d*x) + b)) - 2*b*e*(e*sin(c + d*x))**(m - 1)*(a*(cos(c + d*x) + 1)/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*(-a*(1 - cos(c + d*x))/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*appellf1(1 - m, sympy.S.Half - m/2, sympy.S.Half - m/2, 2 - m, -(a - b)/(a*cos(c + d*x) + b), (a + b)/(a*cos(c + d*x) + b))/(a**3*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_262():
    f = (e*sin(c + d*x))**m/(a + b*sec(c + d*x))**3
    F = (e*sin(c + d*x))**(m + 1)*cos(c + d*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)/(a**3*d*e*(m + 1)*sqrt(cos(c + d*x)**2)) - b**3*e*(e*sin(c + d*x))**(m - 1)*(a*(cos(c + d*x) + 1)/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*(-a*(1 - cos(c + d*x))/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*appellf1(3 - m, sympy.S.Half - m/2, sympy.S.Half - m/2, 4 - m, -(a - b)/(a*cos(c + d*x) + b), (a + b)/(a*cos(c + d*x) + b))/(a**4*d*(3 - m)*(a*cos(c + d*x) + b)**2) + 3*b**2*e*(e*sin(c + d*x))**(m - 1)*(a*(cos(c + d*x) + 1)/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*(-a*(1 - cos(c + d*x))/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*appellf1(2 - m, sympy.S.Half - m/2, sympy.S.Half - m/2, 3 - m, -(a - b)/(a*cos(c + d*x) + b), (a + b)/(a*cos(c + d*x) + b))/(a**4*d*(2 - m)*(a*cos(c + d*x) + b)) - 3*b*e*(e*sin(c + d*x))**(m - 1)*(a*(cos(c + d*x) + 1)/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*(-a*(1 - cos(c + d*x))/(a*cos(c + d*x) + b))**(sympy.S.Half - m/2)*appellf1(1 - m, sympy.S.Half - m/2, sympy.S.Half - m/2, 2 - m, -(a - b)/(a*cos(c + d*x) + b), (a + b)/(a*cos(c + d*x) + b))/(a**4*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_263():
    f = (e*sin(c + d*x))**m*(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_264():
    f = (e*sin(c + d*x))**m*sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_265():
    f = (e*sin(c + d*x))**m/sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m')) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_266():
    f = (e*sin(c + d*x))**m/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_267():
    f = (e*sin(c + d*x))**m*(a + b*sec(c + d*x))**n
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * ((Symbol('e') * sympy.sin((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_268():
    f = (a + b*sec(c + d*x))**n*sin(c + d*x)**5
    F = b*(a + b*sec(c + d*x))**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a**2*d*(n + 1)) - 2*b**3*(a + b*sec(c + d*x))**(n + 1)*hyper((4, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a**4*d*(n + 1)) + b**5*(a + b*sec(c + d*x))**(n + 1)*hyper((6, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a**6*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_269():
    f = (a + b*sec(c + d*x))**n*sin(c + d*x)**3
    F = (a + b*sec(c + d*x))**(n + 1)*(2*a - b*(2 - n)*sec(c + d*x))*cos(c + d*x)**3/(6*a**2*d) + b*(a + b*sec(c + d*x))**(n + 1)*(6*a**2 - b**2*(n**2 - 3*n + 2))*hyper((2, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(6*a**4*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_270():
    f = (a + b*sec(c + d*x))**n*sin(c + d*x)
    F = b*(a + b*sec(c + d*x))**(n + 1)*hyper((2, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a**2*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_271():
    f = (a + b*sec(c + d*x))**n*csc(c + d*x)
    F = -(a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a + b))/(d*(2*a + 2*b)*(n + 1)) + (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a - b))/(d*(2*a - 2*b)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_272():
    f = (a + b*sec(c + d*x))**n*csc(c + d*x)**3
    F = b*(a + b*sec(c + d*x))**(n + 1)*hyper((2, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a + b))/(4*d*(a + b)**2*(n + 1)) + b*(a + b*sec(c + d*x))**(n + 1)*hyper((2, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a - b))/(4*d*(a - b)**2*(n + 1)) - (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a + b))/(d*(4*a + 4*b)*(n + 1)) + (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a - b))/(d*(4*a - 4*b)*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_273():
    f = (a + b*sec(c + d*x))**n*sin(c + d*x)**4
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(4))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_274():
    f = (a + b*sec(c + d*x))**n*sin(c + d*x)**2
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**(Integer(2))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_275():
    f = (a + b*sec(c + d*x))**n*csc(c + d*x)**2
    F = sqrt(2)*b*n*(a + b*sec(c + d*x))**n*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, 1 - n, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(d*((a + b*sec(c + d*x))/(a + b))**n*(a + b)*sqrt(sec(c + d*x) + 1)) - (a + b*sec(c + d*x))**n*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_276():
    f = (a + b*sec(c + d*x))**n*csc(c + d*x)**4
    F = -sqrt(2)*(a + b*sec(c + d*x))**n*(sec(c + d*x) + 1)**(sympy.S(3)/2)*cot(c + d*x)**3*appellf1(sympy.S(-3)/2, sympy.S(5)/2, -n, sympy.S(-1)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(12*d*((a + b*sec(c + d*x))/(a + b))**n) - 3*sqrt(2)*(a + b*sec(c + d*x))**n*sqrt(sec(c + d*x) + 1)*cot(c + d*x)*appellf1(sympy.S(-1)/2, sympy.S(5)/2, -n, sympy.S.Half, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(4*d*((a + b*sec(c + d*x))/(a + b))**n) + sqrt(2)*(a + b*sec(c + d*x))**n*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(3)/2, -n, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(2*d*((a + b*sec(c + d*x))/(a + b))**n*sqrt(sec(c + d*x) + 1)) + sqrt(2)*(a + b*sec(c + d*x))**n*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(5)/2, -n, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(4*d*((a + b*sec(c + d*x))/(a + b))**n*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_277():
    f = (a + b*sec(c + d*x))**n*sin(c + d*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * (sympy.sin((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_278():
    f = (a + b*sec(c + d*x))**n*sqrt(sin(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_279():
    f = (a + b*sec(c + d*x))**n/sqrt(sin(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * (sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_280():
    f = (a + b*sec(c + d*x))**n/sin(c + d*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * ((sympy.sin((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_281():
    f = (e*csc(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)
    F = a*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atan(sqrt(sin(c + d*x)))/d + a*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atanh(sqrt(sin(c + d*x)))/d + 2*a*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d) - 2*a*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)/(3*d) - 2*a*e**2*sqrt(e*csc(c + d*x))*csc(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_282():
    f = (e*csc(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)
    F = -a*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atan(sqrt(sin(c + d*x)))/d + a*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atanh(sqrt(sin(c + d*x)))/d - 2*a*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/d - 2*a*e*sqrt(e*csc(c + d*x))*cos(c + d*x)/d - 2*a*e*sqrt(e*csc(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_283():
    f = sqrt(e*csc(c + d*x))*(a*sec(c + d*x) + a)
    F = a*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atan(sqrt(sin(c + d*x)))/d + a*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atanh(sqrt(sin(c + d*x)))/d + 2*a*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_284():
    f = (a*sec(c + d*x) + a)/sqrt(e*csc(c + d*x))
    F = -a*atan(sqrt(sin(c + d*x)))/(d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + a*atanh(sqrt(sin(c + d*x)))/(d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + 2*a*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_285():
    f = (a*sec(c + d*x) + a)/(e*csc(c + d*x))**(sympy.S(3)/2)
    F = -2*a*cos(c + d*x)/(3*d*e*sqrt(e*csc(c + d*x))) - 2*a/(d*e*sqrt(e*csc(c + d*x))) + a*atan(sqrt(sin(c + d*x)))/(d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + a*atanh(sqrt(sin(c + d*x)))/(d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + 2*a*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_286():
    f = (a*sec(c + d*x) + a)/(e*csc(c + d*x))**(sympy.S(5)/2)
    F = -2*a*sin(c + d*x)*cos(c + d*x)/(5*d*e**2*sqrt(e*csc(c + d*x))) - 2*a*sin(c + d*x)/(3*d*e**2*sqrt(e*csc(c + d*x))) - a*atan(sqrt(sin(c + d*x)))/(d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + a*atanh(sqrt(sin(c + d*x)))/(d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + 6*a*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_287():
    f = (e*csc(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)**2
    F = 2*a**2*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atan(sqrt(sin(c + d*x)))/d + 2*a**2*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atanh(sqrt(sin(c + d*x)))/d + 7*a**2*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d) + 5*a**2*e**2*sqrt(e*csc(c + d*x))*tan(c + d*x)/(3*d) - 2*a**2*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)/(3*d) - 2*a**2*e**2*sqrt(e*csc(c + d*x))*csc(c + d*x)*sec(c + d*x)/(3*d) - 4*a**2*e**2*sqrt(e*csc(c + d*x))*csc(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_288():
    f = (e*csc(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)**2
    F = -2*a**2*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atan(sqrt(sin(c + d*x)))/d + 2*a**2*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atanh(sqrt(sin(c + d*x)))/d - 5*a**2*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/d + 3*a**2*e*sqrt(e*csc(c + d*x))*sin(c + d*x)*tan(c + d*x)/d - 2*a**2*e*sqrt(e*csc(c + d*x))*cos(c + d*x)/d - 2*a**2*e*sqrt(e*csc(c + d*x))*sec(c + d*x)/d - 4*a**2*e*sqrt(e*csc(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_289():
    f = sqrt(e*csc(c + d*x))*(a*sec(c + d*x) + a)**2
    F = 2*a**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atan(sqrt(sin(c + d*x)))/d + 2*a**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*atanh(sqrt(sin(c + d*x)))/d + 3*a**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/d + a**2*sqrt(e*csc(c + d*x))*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_290():
    f = (a*sec(c + d*x) + a)**2/sqrt(e*csc(c + d*x))
    F = a**2*tan(c + d*x)/(d*sqrt(e*csc(c + d*x))) - 2*a**2*atan(sqrt(sin(c + d*x)))/(d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + 2*a**2*atanh(sqrt(sin(c + d*x)))/(d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + a**2*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_291():
    f = (a*sec(c + d*x) + a)**2/(e*csc(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*cos(c + d*x)/(3*d*e*sqrt(e*csc(c + d*x))) + a**2*sec(c + d*x)/(d*e*sqrt(e*csc(c + d*x))) - 4*a**2/(d*e*sqrt(e*csc(c + d*x))) + 2*a**2*atan(sqrt(sin(c + d*x)))/(d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + 2*a**2*atanh(sqrt(sin(c + d*x)))/(d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) - a**2*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_292():
    f = (a*sec(c + d*x) + a)**2/(e*csc(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*sin(c + d*x)*cos(c + d*x)/(5*d*e**2*sqrt(e*csc(c + d*x))) - 4*a**2*sin(c + d*x)/(3*d*e**2*sqrt(e*csc(c + d*x))) + a**2*tan(c + d*x)/(d*e**2*sqrt(e*csc(c + d*x))) - 2*a**2*atan(sqrt(sin(c + d*x)))/(d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) + 2*a**2*atanh(sqrt(sin(c + d*x)))/(d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))) - 9*a**2*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_293():
    f = (e*csc(c + d*x))**(sympy.S(5)/2)/(a*sec(c + d*x) + a)
    F = 4*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a*d) + 2*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)*csc(c + d*x)**2/(7*a*d) - 4*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)/(21*a*d) - 2*e**2*sqrt(e*csc(c + d*x))*csc(c + d*x)**3/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_294():
    f = (e*csc(c + d*x))**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = -4*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a*d) - 4*e*sqrt(e*csc(c + d*x))*cos(c + d*x)/(5*a*d) + 2*e*sqrt(e*csc(c + d*x))*cot(c + d*x)*csc(c + d*x)/(5*a*d) - 2*e*sqrt(e*csc(c + d*x))*csc(c + d*x)**2/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_295():
    f = sqrt(e*csc(c + d*x))/(a*sec(c + d*x) + a)
    F = 4*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*a*d) + 2*sqrt(e*csc(c + d*x))*cot(c + d*x)/(3*a*d) - 2*sqrt(e*csc(c + d*x))*csc(c + d*x)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_296():
    f = 1/(sqrt(e*csc(c + d*x))*(a*sec(c + d*x) + a))
    F = 2*cot(c + d*x)/(a*d*sqrt(e*csc(c + d*x))) - 2*csc(c + d*x)/(a*d*sqrt(e*csc(c + d*x))) + 4*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(a*d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_297():
    f = 1/((e*csc(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a))
    F = -2*cos(c + d*x)/(3*a*d*e*sqrt(e*csc(c + d*x))) + 2/(a*d*e*sqrt(e*csc(c + d*x))) - 4*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(3*a*d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_298():
    f = 1/((e*csc(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a))
    F = -2*sin(c + d*x)*cos(c + d*x)/(5*a*d*e**2*sqrt(e*csc(c + d*x))) + 2*sin(c + d*x)/(3*a*d*e**2*sqrt(e*csc(c + d*x))) - 4*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a*d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_299():
    f = 1/((e*csc(c + d*x))**(sympy.S(7)/2)*(a*sec(c + d*x) + a))
    F = 2*sin(c + d*x)**2/(5*a*d*e**3*sqrt(e*csc(c + d*x))) + 2*cos(c + d*x)**3/(7*a*d*e**3*sqrt(e*csc(c + d*x))) - 2*cos(c + d*x)/(21*a*d*e**3*sqrt(e*csc(c + d*x))) - 4*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a*d*e**3*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_300():
    f = (e*csc(c + d*x))**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**2
    F = 4*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(231*a**2*d) - 2*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)**3*csc(c + d*x)**2/(11*a**2*d) - 2*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)*csc(c + d*x)**4/(11*a**2*d) + 16*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)*csc(c + d*x)**2/(77*a**2*d) - 4*e**2*sqrt(e*csc(c + d*x))*cot(c + d*x)/(231*a**2*d) + 4*e**2*sqrt(e*csc(c + d*x))*csc(c + d*x)**5/(11*a**2*d) - 4*e**2*sqrt(e*csc(c + d*x))*csc(c + d*x)**3/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_301():
    f = (e*csc(c + d*x))**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**2
    F = -4*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(15*a**2*d) - 4*e*sqrt(e*csc(c + d*x))*cos(c + d*x)/(15*a**2*d) - 2*e*sqrt(e*csc(c + d*x))*cot(c + d*x)**3*csc(c + d*x)/(9*a**2*d) - 2*e*sqrt(e*csc(c + d*x))*cot(c + d*x)*csc(c + d*x)**3/(9*a**2*d) + 16*e*sqrt(e*csc(c + d*x))*cot(c + d*x)*csc(c + d*x)/(45*a**2*d) + 4*e*sqrt(e*csc(c + d*x))*csc(c + d*x)**4/(9*a**2*d) - 4*e*sqrt(e*csc(c + d*x))*csc(c + d*x)**2/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_302():
    f = sqrt(e*csc(c + d*x))/(a*sec(c + d*x) + a)**2
    F = 20*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x))*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a**2*d) - 2*sqrt(e*csc(c + d*x))*cot(c + d*x)**3/(7*a**2*d) - 2*sqrt(e*csc(c + d*x))*cot(c + d*x)*csc(c + d*x)**2/(7*a**2*d) + 16*sqrt(e*csc(c + d*x))*cot(c + d*x)/(21*a**2*d) + 4*sqrt(e*csc(c + d*x))*csc(c + d*x)**3/(7*a**2*d) - 4*sqrt(e*csc(c + d*x))*csc(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_303():
    f = 1/(sqrt(e*csc(c + d*x))*(a*sec(c + d*x) + a)**2)
    F = -2*cot(c + d*x)**3/(5*a**2*d*sqrt(e*csc(c + d*x))) - 2*cot(c + d*x)*csc(c + d*x)**2/(5*a**2*d*sqrt(e*csc(c + d*x))) + 16*cot(c + d*x)/(5*a**2*d*sqrt(e*csc(c + d*x))) + 4*csc(c + d*x)**3/(5*a**2*d*sqrt(e*csc(c + d*x))) - 4*csc(c + d*x)/(a**2*d*sqrt(e*csc(c + d*x))) + 28*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a**2*d*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_304():
    f = 1/((e*csc(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)**2)
    F = -2*cos(c + d*x)*cot(c + d*x)**2/(3*a**2*d*e*sqrt(e*csc(c + d*x))) - 4*cos(c + d*x)/(3*a**2*d*e*sqrt(e*csc(c + d*x))) - 2*cot(c + d*x)*csc(c + d*x)/(3*a**2*d*e*sqrt(e*csc(c + d*x))) + 4*csc(c + d*x)**2/(3*a**2*d*e*sqrt(e*csc(c + d*x))) + 4/(a**2*d*e*sqrt(e*csc(c + d*x))) - 4*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(a**2*d*e*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_305():
    f = 1/((e*csc(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)**2)
    F = -12*sin(c + d*x)*cos(c + d*x)/(5*a**2*d*e**2*sqrt(e*csc(c + d*x))) + 4*sin(c + d*x)/(3*a**2*d*e**2*sqrt(e*csc(c + d*x))) - 2*cos(c + d*x)**2*cot(c + d*x)/(a**2*d*e**2*sqrt(e*csc(c + d*x))) - 2*cot(c + d*x)/(a**2*d*e**2*sqrt(e*csc(c + d*x))) + 4*csc(c + d*x)/(a**2*d*e**2*sqrt(e*csc(c + d*x))) - 44*elliptic_e(c/2 + d*x/2 - pi/4, 2)/(5*a**2*d*e**2*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_3_d_sin_pow_n_a_plus_b_sec_pow_m_306():
    f = 1/((e*csc(c + d*x))**(sympy.S(7)/2)*(a*sec(c + d*x) + a)**2)
    F = 4*sin(c + d*x)**2/(5*a**2*d*e**3*sqrt(e*csc(c + d*x))) + 2*cos(c + d*x)**3/(7*a**2*d*e**3*sqrt(e*csc(c + d*x))) + 26*cos(c + d*x)/(21*a**2*d*e**3*sqrt(e*csc(c + d*x))) - 4/(a**2*d*e**3*sqrt(e*csc(c + d*x))) + 52*elliptic_f(c/2 + d*x/2 - pi/4, 2)/(21*a**2*d*e**3*sqrt(e*csc(c + d*x))*sqrt(sin(c + d*x)))
    assert integrate(f, x) == F

