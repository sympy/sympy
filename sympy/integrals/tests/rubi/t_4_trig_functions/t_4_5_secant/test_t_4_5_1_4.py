"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.1.4 (d tan)^n (a+b sec)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_1():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**9
    F = -a*log(cos(c + d*x))/d + a*sec(c + d*x)**9/(9*d) + a*sec(c + d*x)**8/(8*d) - 4*a*sec(c + d*x)**7/(7*d) - 2*a*sec(c + d*x)**6/(3*d) + 6*a*sec(c + d*x)**5/(5*d) + 3*a*sec(c + d*x)**4/(2*d) - 4*a*sec(c + d*x)**3/(3*d) - 2*a*sec(c + d*x)**2/d + a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_2():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**7
    F = a*log(cos(c + d*x))/d + a*sec(c + d*x)**7/(7*d) + a*sec(c + d*x)**6/(6*d) - 3*a*sec(c + d*x)**5/(5*d) - 3*a*sec(c + d*x)**4/(4*d) + a*sec(c + d*x)**3/d + 3*a*sec(c + d*x)**2/(2*d) - a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_3():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**5
    F = -a*log(cos(c + d*x))/d + a*sec(c + d*x)**5/(5*d) + a*sec(c + d*x)**4/(4*d) - 2*a*sec(c + d*x)**3/(3*d) - a*sec(c + d*x)**2/d + a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_4():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**3
    F = a*log(cos(c + d*x))/d + a*sec(c + d*x)**3/(3*d) + a*sec(c + d*x)**2/(2*d) - a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_5():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)
    F = -a*log(cos(c + d*x))/d + a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_6():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)
    F = a*log(1 - cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_7():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**3
    F = -3*a*log(1 - cos(c + d*x))/(4*d) - a*log(cos(c + d*x) + 1)/(4*d) - a/(2*d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_8():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**5
    F = 11*a*log(1 - cos(c + d*x))/(16*d) + 5*a*log(cos(c + d*x) + 1)/(16*d) + a/(8*d*(cos(c + d*x) + 1)) + 3*a/(4*d*(1 - cos(c + d*x))) - a/(8*d*(1 - cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_9():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**7
    F = -21*a*log(1 - cos(c + d*x))/(32*d) - 11*a*log(cos(c + d*x) + 1)/(32*d) - a/(4*d*(cos(c + d*x) + 1)) + a/(32*d*(cos(c + d*x) + 1)**2) - 15*a/(16*d*(1 - cos(c + d*x))) + 9*a/(32*d*(1 - cos(c + d*x))**2) - a/(24*d*(1 - cos(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_10():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**8
    F = a*x + 35*a*atanh(sin(c + d*x))/(128*d) + (7*a*sec(c + d*x) + 8*a)*tan(c + d*x)**7/(56*d) - (35*a*sec(c + d*x) + 48*a)*tan(c + d*x)**5/(240*d) + (35*a*sec(c + d*x) + 64*a)*tan(c + d*x)**3/(192*d) - (35*a*sec(c + d*x) + 128*a)*tan(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_11():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**6
    F = -a*x - 5*a*atanh(sin(c + d*x))/(16*d) + (5*a*sec(c + d*x) + 6*a)*tan(c + d*x)**5/(30*d) - (5*a*sec(c + d*x) + 8*a)*tan(c + d*x)**3/(24*d) + (5*a*sec(c + d*x) + 16*a)*tan(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_12():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**4
    F = a*x + 3*a*atanh(sin(c + d*x))/(8*d) + (3*a*sec(c + d*x) + 4*a)*tan(c + d*x)**3/(12*d) - (3*a*sec(c + d*x) + 8*a)*tan(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_13():
    f = (a*sec(c + d*x) + a)*tan(c + d*x)**2
    F = -a*x - a*atanh(sin(c + d*x))/(2*d) + (a*sec(c + d*x) + 2*a)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_14():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**2
    F = -a*x - (a*sec(c + d*x) + a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_15():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**4
    F = a*x - (a*sec(c + d*x) + a)*cot(c + d*x)**3/(3*d) + (2*a*sec(c + d*x) + 3*a)*cot(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_16():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**6
    F = -a*x - (a*sec(c + d*x) + a)*cot(c + d*x)**5/(5*d) + (4*a*sec(c + d*x) + 5*a)*cot(c + d*x)**3/(15*d) - (8*a*sec(c + d*x) + 15*a)*cot(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_17():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**8
    F = a*x - (a*sec(c + d*x) + a)*cot(c + d*x)**7/(7*d) + (6*a*sec(c + d*x) + 7*a)*cot(c + d*x)**5/(35*d) + (16*a*sec(c + d*x) + 35*a)*cot(c + d*x)/(35*d) - (24*a*sec(c + d*x) + 35*a)*cot(c + d*x)**3/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_18():
    f = (a*sec(c + d*x) + a)*cot(c + d*x)**10
    F = -a*x - (a*sec(c + d*x) + a)*cot(c + d*x)**9/(9*d) + (8*a*sec(c + d*x) + 9*a)*cot(c + d*x)**7/(63*d) - (16*a*sec(c + d*x) + 21*a)*cot(c + d*x)**5/(105*d) + (64*a*sec(c + d*x) + 105*a)*cot(c + d*x)**3/(315*d) - (128*a*sec(c + d*x) + 315*a)*cot(c + d*x)/(315*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_19():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)**9
    F = -a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)**10/(10*d) + 2*a**2*sec(c + d*x)**9/(9*d) - 3*a**2*sec(c + d*x)**8/(8*d) - 8*a**2*sec(c + d*x)**7/(7*d) + a**2*sec(c + d*x)**6/(3*d) + 12*a**2*sec(c + d*x)**5/(5*d) + a**2*sec(c + d*x)**4/(2*d) - 8*a**2*sec(c + d*x)**3/(3*d) - 3*a**2*sec(c + d*x)**2/(2*d) + 2*a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_20():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)**7
    F = a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)**8/(8*d) + 2*a**2*sec(c + d*x)**7/(7*d) - a**2*sec(c + d*x)**6/(3*d) - 6*a**2*sec(c + d*x)**5/(5*d) + 2*a**2*sec(c + d*x)**3/d + a**2*sec(c + d*x)**2/d - 2*a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_21():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)**5
    F = -a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)**6/(6*d) + 2*a**2*sec(c + d*x)**5/(5*d) - a**2*sec(c + d*x)**4/(4*d) - 4*a**2*sec(c + d*x)**3/(3*d) - a**2*sec(c + d*x)**2/(2*d) + 2*a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_22():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)**3
    F = a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)**4/(4*d) + 2*a**2*sec(c + d*x)**3/(3*d) - 2*a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_23():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)
    F = -a**2*log(cos(c + d*x))/d + a**2*sec(c + d*x)**2/(2*d) + 2*a**2*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_24():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)
    F = 2*a**2*log(1 - cos(c + d*x))/d - a**2*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_25():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**3
    F = -a**2*log(1 - cos(c + d*x))/d - a**2/(d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_26():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**5
    F = 7*a**2*log(1 - cos(c + d*x))/(8*d) + a**2*log(cos(c + d*x) + 1)/(8*d) + 5*a**2/(4*d*(1 - cos(c + d*x))) - a**2/(4*d*(1 - cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_27():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**7
    F = -13*a**2*log(1 - cos(c + d*x))/(16*d) - 3*a**2*log(cos(c + d*x) + 1)/(16*d) - a**2/(16*d*(cos(c + d*x) + 1)) - 23*a**2/(16*d*(1 - cos(c + d*x))) + a**2/(2*d*(1 - cos(c + d*x))**2) - a**2/(12*d*(1 - cos(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_28():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**9
    F = 99*a**2*log(1 - cos(c + d*x))/(128*d) + 29*a**2*log(cos(c + d*x) + 1)/(128*d) + 9*a**2/(64*d*(cos(c + d*x) + 1)) - a**2/(64*d*(cos(c + d*x) + 1)**2) + 51*a**2/(32*d*(1 - cos(c + d*x))) - 3*a**2/(4*d*(1 - cos(c + d*x))**2) + 11*a**2/(48*d*(1 - cos(c + d*x))**3) - a**2/(32*d*(1 - cos(c + d*x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_29():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)**6
    F = -a**2*x + a**2*tan(c + d*x)**7/(7*d) + a**2*tan(c + d*x)**5*sec(c + d*x)/(3*d) + a**2*tan(c + d*x)**5/(5*d) - 5*a**2*tan(c + d*x)**3*sec(c + d*x)/(12*d) - a**2*tan(c + d*x)**3/(3*d) + 5*a**2*tan(c + d*x)*sec(c + d*x)/(8*d) + a**2*tan(c + d*x)/d - 5*a**2*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_30():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)**4
    F = a**2*x + a**2*tan(c + d*x)**5/(5*d) + a**2*tan(c + d*x)**3*sec(c + d*x)/(2*d) + a**2*tan(c + d*x)**3/(3*d) - 3*a**2*tan(c + d*x)*sec(c + d*x)/(4*d) - a**2*tan(c + d*x)/d + 3*a**2*atanh(sin(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_31():
    f = (a*sec(c + d*x) + a)**2*tan(c + d*x)**2
    F = -a**2*x + a**2*tan(c + d*x)**3/(3*d) + a**2*tan(c + d*x)*sec(c + d*x)/d + a**2*tan(c + d*x)/d - a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_32():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**2
    F = -a**2*x - 2*a**2*cot(c + d*x)/d - 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_33():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**4
    F = a**2*x - 2*a**2*cot(c + d*x)**3/(3*d) + a**2*cot(c + d*x)/d - 2*a**2*csc(c + d*x)**3/(3*d) + 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_34():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**6
    F = -a**2*x - 2*a**2*cot(c + d*x)**5/(5*d) + a**2*cot(c + d*x)**3/(3*d) - a**2*cot(c + d*x)/d - 2*a**2*csc(c + d*x)**5/(5*d) + 4*a**2*csc(c + d*x)**3/(3*d) - 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_35():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**8
    F = a**2*x - 2*a**2*cot(c + d*x)**7/(7*d) + a**2*cot(c + d*x)**5/(5*d) - a**2*cot(c + d*x)**3/(3*d) + a**2*cot(c + d*x)/d - 2*a**2*csc(c + d*x)**7/(7*d) + 6*a**2*csc(c + d*x)**5/(5*d) - 2*a**2*csc(c + d*x)**3/d + 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_36():
    f = (a*sec(c + d*x) + a)**2*cot(c + d*x)**10
    F = -a**2*x - 2*a**2*cot(c + d*x)**9/(9*d) + a**2*cot(c + d*x)**7/(7*d) - a**2*cot(c + d*x)**5/(5*d) + a**2*cot(c + d*x)**3/(3*d) - a**2*cot(c + d*x)/d - 2*a**2*csc(c + d*x)**9/(9*d) + 8*a**2*csc(c + d*x)**7/(7*d) - 12*a**2*csc(c + d*x)**5/(5*d) + 8*a**2*csc(c + d*x)**3/(3*d) - 2*a**2*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_37():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)**9
    F = -a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**11/(11*d) + 3*a**3*sec(c + d*x)**10/(10*d) - a**3*sec(c + d*x)**9/(9*d) - 11*a**3*sec(c + d*x)**8/(8*d) - 6*a**3*sec(c + d*x)**7/(7*d) + 7*a**3*sec(c + d*x)**6/(3*d) + 14*a**3*sec(c + d*x)**5/(5*d) - 3*a**3*sec(c + d*x)**4/(2*d) - 11*a**3*sec(c + d*x)**3/(3*d) - a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_38():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)**7
    F = a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**9/(9*d) + 3*a**3*sec(c + d*x)**8/(8*d) - 4*a**3*sec(c + d*x)**6/(3*d) - 6*a**3*sec(c + d*x)**5/(5*d) + 3*a**3*sec(c + d*x)**4/(2*d) + 8*a**3*sec(c + d*x)**3/(3*d) - 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_39():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)**5
    F = -a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**7/(7*d) + a**3*sec(c + d*x)**6/(2*d) + a**3*sec(c + d*x)**5/(5*d) - 5*a**3*sec(c + d*x)**4/(4*d) - 5*a**3*sec(c + d*x)**3/(3*d) + a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_40():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)**3
    F = a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**5/(5*d) + 3*a**3*sec(c + d*x)**4/(4*d) + 2*a**3*sec(c + d*x)**3/(3*d) - a**3*sec(c + d*x)**2/d - 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_41():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)
    F = -a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)**3/(3*d) + 3*a**3*sec(c + d*x)**2/(2*d) + 3*a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_42():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)
    F = 4*a**3*log(1 - cos(c + d*x))/d - 3*a**3*log(cos(c + d*x))/d + a**3*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_43():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**3
    F = -a**3*log(1 - cos(c + d*x))/d - 2*a**3/(d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_44():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**5
    F = a**3*log(1 - cos(c + d*x))/d + 2*a**3/(d*(1 - cos(c + d*x))) - a**3/(2*d*(1 - cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_45():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**7
    F = -15*a**3*log(1 - cos(c + d*x))/(16*d) - a**3*log(cos(c + d*x) + 1)/(16*d) - 17*a**3/(8*d*(1 - cos(c + d*x))) + 7*a**3/(8*d*(1 - cos(c + d*x))**2) - a**3/(6*d*(1 - cos(c + d*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_46():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**9
    F = 57*a**3*log(1 - cos(c + d*x))/(64*d) + 7*a**3*log(cos(c + d*x) + 1)/(64*d) + a**3/(32*d*(cos(c + d*x) + 1)) + 9*a**3/(4*d*(1 - cos(c + d*x))) - 39*a**3/(32*d*(1 - cos(c + d*x))**2) + 5*a**3/(12*d*(1 - cos(c + d*x))**3) - a**3/(16*d*(1 - cos(c + d*x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_47():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)**6
    F = -a**3*x + 3*a**3*tan(c + d*x)**7/(7*d) + a**3*tan(c + d*x)**5*sec(c + d*x)**3/(8*d) + a**3*tan(c + d*x)**5*sec(c + d*x)/(2*d) + a**3*tan(c + d*x)**5/(5*d) - 5*a**3*tan(c + d*x)**3*sec(c + d*x)**3/(48*d) - 5*a**3*tan(c + d*x)**3*sec(c + d*x)/(8*d) - a**3*tan(c + d*x)**3/(3*d) + 5*a**3*tan(c + d*x)*sec(c + d*x)**3/(64*d) + 115*a**3*tan(c + d*x)*sec(c + d*x)/(128*d) + a**3*tan(c + d*x)/d - 125*a**3*atanh(sin(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_48():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)**4
    F = a**3*x + 3*a**3*tan(c + d*x)**5/(5*d) + a**3*tan(c + d*x)**3*sec(c + d*x)**3/(6*d) + 3*a**3*tan(c + d*x)**3*sec(c + d*x)/(4*d) + a**3*tan(c + d*x)**3/(3*d) - a**3*tan(c + d*x)*sec(c + d*x)**3/(8*d) - 17*a**3*tan(c + d*x)*sec(c + d*x)/(16*d) - a**3*tan(c + d*x)/d + 19*a**3*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_49():
    f = (a*sec(c + d*x) + a)**3*tan(c + d*x)**2
    F = -a**3*x + a**3*tan(c + d*x)**3/d + a**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 11*a**3*tan(c + d*x)*sec(c + d*x)/(8*d) + a**3*tan(c + d*x)/d - 13*a**3*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_50():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**2
    F = -a**3*x - 4*a**3*cot(c + d*x)/d + a**3*atanh(sin(c + d*x))/d - 4*a**3*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_51():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**4
    F = a**3*x - 4*a**3*cot(c + d*x)**3/(3*d) + a**3*cot(c + d*x)/d - 4*a**3*csc(c + d*x)**3/(3*d) + 3*a**3*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_52():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**6
    F = -a**3*x - 4*a**3*cot(c + d*x)**5/(5*d) + a**3*cot(c + d*x)**3/(3*d) - a**3*cot(c + d*x)/d - 4*a**3*csc(c + d*x)**5/(5*d) + 7*a**3*csc(c + d*x)**3/(3*d) - 3*a**3*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_53():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**8
    F = a**3*x - 4*a**3*cot(c + d*x)**7/(7*d) + a**3*cot(c + d*x)**5/(5*d) - a**3*cot(c + d*x)**3/(3*d) + a**3*cot(c + d*x)/d - 4*a**3*csc(c + d*x)**7/(7*d) + 11*a**3*csc(c + d*x)**5/(5*d) - 10*a**3*csc(c + d*x)**3/(3*d) + 3*a**3*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_54():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**10
    F = -a**3*x - 4*a**3*cot(c + d*x)**9/(9*d) + a**3*cot(c + d*x)**7/(7*d) - a**3*cot(c + d*x)**5/(5*d) + a**3*cot(c + d*x)**3/(3*d) - a**3*cot(c + d*x)/d - 4*a**3*csc(c + d*x)**9/(9*d) + 15*a**3*csc(c + d*x)**7/(7*d) - 21*a**3*csc(c + d*x)**5/(5*d) + 13*a**3*csc(c + d*x)**3/(3*d) - 3*a**3*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_55():
    f = (a*sec(c + d*x) + a)**3*cot(c + d*x)**12
    F = a**3*x - 4*a**3*cot(c + d*x)**11/(11*d) + a**3*cot(c + d*x)**9/(9*d) - a**3*cot(c + d*x)**7/(7*d) + a**3*cot(c + d*x)**5/(5*d) - a**3*cot(c + d*x)**3/(3*d) + a**3*cot(c + d*x)/d - 4*a**3*csc(c + d*x)**11/(11*d) + 19*a**3*csc(c + d*x)**9/(9*d) - 36*a**3*csc(c + d*x)**7/(7*d) + 34*a**3*csc(c + d*x)**5/(5*d) - 16*a**3*csc(c + d*x)**3/(3*d) + 3*a**3*csc(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_56():
    f = tan(c + d*x)**9/(a*sec(c + d*x) + a)
    F = -log(cos(c + d*x))/(a*d) + sec(c + d*x)**7/(7*a*d) - sec(c + d*x)**6/(6*a*d) - 3*sec(c + d*x)**5/(5*a*d) + 3*sec(c + d*x)**4/(4*a*d) + sec(c + d*x)**3/(a*d) - 3*sec(c + d*x)**2/(2*a*d) - sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_57():
    f = tan(c + d*x)**7/(a*sec(c + d*x) + a)
    F = log(cos(c + d*x))/(a*d) + sec(c + d*x)**5/(5*a*d) - sec(c + d*x)**4/(4*a*d) - 2*sec(c + d*x)**3/(3*a*d) + sec(c + d*x)**2/(a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_58():
    f = tan(c + d*x)**5/(a*sec(c + d*x) + a)
    F = -log(cos(c + d*x))/(a*d) + sec(c + d*x)**3/(3*a*d) - sec(c + d*x)**2/(2*a*d) - sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_59():
    f = tan(c + d*x)**3/(a*sec(c + d*x) + a)
    F = log(cos(c + d*x))/(a*d) + sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_60():
    f = tan(c + d*x)/(a*sec(c + d*x) + a)
    F = -log(cos(c + d*x) + 1)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_61():
    f = cot(c + d*x)/(a*sec(c + d*x) + a)
    F = log(1 - cos(c + d*x))/(4*a*d) + 3*log(cos(c + d*x) + 1)/(4*a*d) + 1/(2*a*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_62():
    f = cot(c + d*x)**3/(a*sec(c + d*x) + a)
    F = -5*log(1 - cos(c + d*x))/(16*a*d) - 11*log(cos(c + d*x) + 1)/(16*a*d) - 3/(4*a*d*(cos(c + d*x) + 1)) + 1/(8*a*d*(cos(c + d*x) + 1)**2) - 1/(8*a*d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_63():
    f = cot(c + d*x)**5/(a*sec(c + d*x) + a)
    F = 11*log(1 - cos(c + d*x))/(32*a*d) + 21*log(cos(c + d*x) + 1)/(32*a*d) + 15/(16*a*d*(cos(c + d*x) + 1)) - 9/(32*a*d*(cos(c + d*x) + 1)**2) + 1/(24*a*d*(cos(c + d*x) + 1)**3) + 1/(4*a*d*(1 - cos(c + d*x))) - 1/(32*a*d*(1 - cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_64():
    f = tan(c + d*x)**8/(a*sec(c + d*x) + a)
    F = x/a - (6 - 5*sec(c + d*x))*tan(c + d*x)**5/(30*a*d) + (8 - 5*sec(c + d*x))*tan(c + d*x)**3/(24*a*d) - (16 - 5*sec(c + d*x))*tan(c + d*x)/(16*a*d) - 5*atanh(sin(c + d*x))/(16*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_65():
    f = tan(c + d*x)**6/(a*sec(c + d*x) + a)
    F = -x/a - (4 - 3*sec(c + d*x))*tan(c + d*x)**3/(12*a*d) + (8 - 3*sec(c + d*x))*tan(c + d*x)/(8*a*d) + 3*atanh(sin(c + d*x))/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_66():
    f = tan(c + d*x)**4/(a*sec(c + d*x) + a)
    F = x/a - (2 - sec(c + d*x))*tan(c + d*x)/(2*a*d) - atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_67():
    f = tan(c + d*x)**2/(a*sec(c + d*x) + a)
    F = -x/a + atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_68():
    f = cot(c + d*x)**2/(a*sec(c + d*x) + a)
    F = -x/a + (1 - sec(c + d*x))*cot(c + d*x)**3/(3*a*d) - (3 - 2*sec(c + d*x))*cot(c + d*x)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_69():
    f = cot(c + d*x)**4/(a*sec(c + d*x) + a)
    F = x/a + (1 - sec(c + d*x))*cot(c + d*x)**5/(5*a*d) - (5 - 4*sec(c + d*x))*cot(c + d*x)**3/(15*a*d) + (15 - 8*sec(c + d*x))*cot(c + d*x)/(15*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_70():
    f = cot(c + d*x)**6/(a*sec(c + d*x) + a)
    F = -x/a + (1 - sec(c + d*x))*cot(c + d*x)**7/(7*a*d) - (7 - 6*sec(c + d*x))*cot(c + d*x)**5/(35*a*d) + (35 - 24*sec(c + d*x))*cot(c + d*x)**3/(105*a*d) - (35 - 16*sec(c + d*x))*cot(c + d*x)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_71():
    f = tan(c + d*x)**9/(a*sec(c + d*x) + a)**2
    F = -log(cos(c + d*x))/(a**2*d) + sec(c + d*x)**6/(6*a**2*d) - 2*sec(c + d*x)**5/(5*a**2*d) - sec(c + d*x)**4/(4*a**2*d) + 4*sec(c + d*x)**3/(3*a**2*d) - sec(c + d*x)**2/(2*a**2*d) - 2*sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_72():
    f = tan(c + d*x)**7/(a*sec(c + d*x) + a)**2
    F = log(cos(c + d*x))/(a**2*d) + sec(c + d*x)**4/(4*a**2*d) - 2*sec(c + d*x)**3/(3*a**2*d) + 2*sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_73():
    f = tan(c + d*x)**5/(a*sec(c + d*x) + a)**2
    F = -log(cos(c + d*x))/(a**2*d) + sec(c + d*x)**2/(2*a**2*d) - 2*sec(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_74():
    f = tan(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = 2*log(cos(c + d*x) + 1)/(a**2*d) - log(cos(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_75():
    f = tan(c + d*x)/(a*sec(c + d*x) + a)**2
    F = -log(cos(c + d*x) + 1)/(a**2*d) - 1/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_76():
    f = cot(c + d*x)/(a*sec(c + d*x) + a)**2
    F = log(1 - cos(c + d*x))/(8*a**2*d) + 7*log(cos(c + d*x) + 1)/(8*a**2*d) + 5/(4*a**2*d*(cos(c + d*x) + 1)) - 1/(4*a**2*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_77():
    f = cot(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = -3*log(1 - cos(c + d*x))/(16*a**2*d) - 13*log(cos(c + d*x) + 1)/(16*a**2*d) - 23/(16*a**2*d*(cos(c + d*x) + 1)) + 1/(2*a**2*d*(cos(c + d*x) + 1)**2) - 1/(12*a**2*d*(cos(c + d*x) + 1)**3) - 1/(16*a**2*d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_78():
    f = cot(c + d*x)**5/(a*sec(c + d*x) + a)**2
    F = 29*log(1 - cos(c + d*x))/(128*a**2*d) + 99*log(cos(c + d*x) + 1)/(128*a**2*d) + 51/(32*a**2*d*(cos(c + d*x) + 1)) - 3/(4*a**2*d*(cos(c + d*x) + 1)**2) + 11/(48*a**2*d*(cos(c + d*x) + 1)**3) - 1/(32*a**2*d*(cos(c + d*x) + 1)**4) + 9/(64*a**2*d*(1 - cos(c + d*x))) - 1/(64*a**2*d*(1 - cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_79():
    f = tan(c + d*x)**8/(a*sec(c + d*x) + a)**2
    F = x/a**2 + tan(c + d*x)**5/(5*a**2*d) - tan(c + d*x)**3*sec(c + d*x)/(2*a**2*d) + tan(c + d*x)**3/(3*a**2*d) + 3*tan(c + d*x)*sec(c + d*x)/(4*a**2*d) - tan(c + d*x)/(a**2*d) - 3*atanh(sin(c + d*x))/(4*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_80():
    f = tan(c + d*x)**6/(a*sec(c + d*x) + a)**2
    F = -x/a**2 + tan(c + d*x)**3/(3*a**2*d) - tan(c + d*x)*sec(c + d*x)/(a**2*d) + tan(c + d*x)/(a**2*d) + atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_81():
    f = tan(c + d*x)**4/(a*sec(c + d*x) + a)**2
    F = x/a**2 + tan(c + d*x)/(a**2*d) - 2*atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_82():
    f = cot(c + d*x)**2/(a*sec(c + d*x) + a)**2
    F = -x/a**2 - 2*cot(c + d*x)**5/(5*a**2*d) + cot(c + d*x)**3/(3*a**2*d) - cot(c + d*x)/(a**2*d) + 2*csc(c + d*x)**5/(5*a**2*d) - 4*csc(c + d*x)**3/(3*a**2*d) + 2*csc(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_83():
    f = cot(c + d*x)**4/(a*sec(c + d*x) + a)**2
    F = x/a**2 - 2*cot(c + d*x)**7/(7*a**2*d) + cot(c + d*x)**5/(5*a**2*d) - cot(c + d*x)**3/(3*a**2*d) + cot(c + d*x)/(a**2*d) + 2*csc(c + d*x)**7/(7*a**2*d) - 6*csc(c + d*x)**5/(5*a**2*d) + 2*csc(c + d*x)**3/(a**2*d) - 2*csc(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_84():
    f = cot(c + d*x)**6/(a*sec(c + d*x) + a)**2
    F = -x/a**2 - 2*cot(c + d*x)**9/(9*a**2*d) + cot(c + d*x)**7/(7*a**2*d) - cot(c + d*x)**5/(5*a**2*d) + cot(c + d*x)**3/(3*a**2*d) - cot(c + d*x)/(a**2*d) + 2*csc(c + d*x)**9/(9*a**2*d) - 8*csc(c + d*x)**7/(7*a**2*d) + 12*csc(c + d*x)**5/(5*a**2*d) - 8*csc(c + d*x)**3/(3*a**2*d) + 2*csc(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_85():
    f = tan(c + d*x)**11/(a*sec(c + d*x) + a)**3
    F = log(cos(c + d*x))/(a**3*d) + sec(c + d*x)**7/(7*a**3*d) - sec(c + d*x)**6/(2*a**3*d) + sec(c + d*x)**5/(5*a**3*d) + 5*sec(c + d*x)**4/(4*a**3*d) - 5*sec(c + d*x)**3/(3*a**3*d) - sec(c + d*x)**2/(2*a**3*d) + 3*sec(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_86():
    f = tan(c + d*x)**9/(a*sec(c + d*x) + a)**3
    F = -log(cos(c + d*x))/(a**3*d) + sec(c + d*x)**5/(5*a**3*d) - 3*sec(c + d*x)**4/(4*a**3*d) + 2*sec(c + d*x)**3/(3*a**3*d) + sec(c + d*x)**2/(a**3*d) - 3*sec(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_87():
    f = tan(c + d*x)**7/(a*sec(c + d*x) + a)**3
    F = log(cos(c + d*x))/(a**3*d) + sec(c + d*x)**3/(3*a**3*d) - 3*sec(c + d*x)**2/(2*a**3*d) + 3*sec(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_88():
    f = tan(c + d*x)**5/(a*sec(c + d*x) + a)**3
    F = -4*log(cos(c + d*x) + 1)/(a**3*d) + 3*log(cos(c + d*x))/(a**3*d) + sec(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_89():
    f = tan(c + d*x)**3/(a*sec(c + d*x) + a)**3
    F = log(cos(c + d*x) + 1)/(a**3*d) + 2/(a**3*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_90():
    f = tan(c + d*x)/(a*sec(c + d*x) + a)**3
    F = -log(cos(c + d*x) + 1)/(a**3*d) - 2/(a**3*d*(cos(c + d*x) + 1)) + 1/(2*a**3*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_91():
    f = cot(c + d*x)/(a*sec(c + d*x) + a)**3
    F = log(1 - cos(c + d*x))/(16*a**3*d) + 15*log(cos(c + d*x) + 1)/(16*a**3*d) + 17/(8*a**3*d*(cos(c + d*x) + 1)) - 7/(8*a**3*d*(cos(c + d*x) + 1)**2) + 1/(6*a**3*d*(cos(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_92():
    f = cot(c + d*x)**3/(a*sec(c + d*x) + a)**3
    F = -7*log(1 - cos(c + d*x))/(64*a**3*d) - 57*log(cos(c + d*x) + 1)/(64*a**3*d) - 9/(4*a**3*d*(cos(c + d*x) + 1)) + 39/(32*a**3*d*(cos(c + d*x) + 1)**2) - 5/(12*a**3*d*(cos(c + d*x) + 1)**3) + 1/(16*a**3*d*(cos(c + d*x) + 1)**4) - 1/(32*a**3*d*(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_93():
    f = cot(c + d*x)**5/(a*sec(c + d*x) + a)**3
    F = 37*log(1 - cos(c + d*x))/(256*a**3*d) + 219*log(cos(c + d*x) + 1)/(256*a**3*d) + 303/(128*a**3*d*(cos(c + d*x) + 1)) - 99/(64*a**3*d*(cos(c + d*x) + 1)**2) + 35/(48*a**3*d*(cos(c + d*x) + 1)**3) - 13/(64*a**3*d*(cos(c + d*x) + 1)**4) + 1/(40*a**3*d*(cos(c + d*x) + 1)**5) + 5/(64*a**3*d*(1 - cos(c + d*x))) - 1/(128*a**3*d*(1 - cos(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_94():
    f = tan(c + d*x)**12/(a*sec(c + d*x) + a)**3
    F = x/a**3 - 3*tan(c + d*x)**7/(7*a**3*d) + tan(c + d*x)**5*sec(c + d*x)**3/(8*a**3*d) + tan(c + d*x)**5*sec(c + d*x)/(2*a**3*d) - tan(c + d*x)**5/(5*a**3*d) - 5*tan(c + d*x)**3*sec(c + d*x)**3/(48*a**3*d) - 5*tan(c + d*x)**3*sec(c + d*x)/(8*a**3*d) + tan(c + d*x)**3/(3*a**3*d) + 5*tan(c + d*x)*sec(c + d*x)**3/(64*a**3*d) + 115*tan(c + d*x)*sec(c + d*x)/(128*a**3*d) - tan(c + d*x)/(a**3*d) - 125*atanh(sin(c + d*x))/(128*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_95():
    f = tan(c + d*x)**10/(a*sec(c + d*x) + a)**3
    F = -x/a**3 - 3*tan(c + d*x)**5/(5*a**3*d) + tan(c + d*x)**3*sec(c + d*x)**3/(6*a**3*d) + 3*tan(c + d*x)**3*sec(c + d*x)/(4*a**3*d) - tan(c + d*x)**3/(3*a**3*d) - tan(c + d*x)*sec(c + d*x)**3/(8*a**3*d) - 17*tan(c + d*x)*sec(c + d*x)/(16*a**3*d) + tan(c + d*x)/(a**3*d) + 19*atanh(sin(c + d*x))/(16*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_96():
    f = tan(c + d*x)**8/(a*sec(c + d*x) + a)**3
    F = x/a**3 - tan(c + d*x)**3/(a**3*d) + tan(c + d*x)*sec(c + d*x)**3/(4*a**3*d) + 11*tan(c + d*x)*sec(c + d*x)/(8*a**3*d) - tan(c + d*x)/(a**3*d) - 13*atanh(sin(c + d*x))/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_97():
    f = tan(c + d*x)**6/(a*sec(c + d*x) + a)**3
    F = -x/a**3 - (1 - sec(c + d*x))*tan(c + d*x)/(2*a**3*d) - 5*tan(c + d*x)/(2*a**3*d) + 7*atanh(sin(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_98():
    f = cot(c + d*x)**2/(a*sec(c + d*x) + a)**3
    F = -x/a**3 + 4*cot(c + d*x)**7/(7*a**3*d) - cot(c + d*x)**5/(5*a**3*d) + cot(c + d*x)**3/(3*a**3*d) - cot(c + d*x)/(a**3*d) - 4*csc(c + d*x)**7/(7*a**3*d) + 11*csc(c + d*x)**5/(5*a**3*d) - 10*csc(c + d*x)**3/(3*a**3*d) + 3*csc(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_99():
    f = cot(c + d*x)**4/(a*sec(c + d*x) + a)**3
    F = x/a**3 + 4*cot(c + d*x)**9/(9*a**3*d) - cot(c + d*x)**7/(7*a**3*d) + cot(c + d*x)**5/(5*a**3*d) - cot(c + d*x)**3/(3*a**3*d) + cot(c + d*x)/(a**3*d) - 4*csc(c + d*x)**9/(9*a**3*d) + 15*csc(c + d*x)**7/(7*a**3*d) - 21*csc(c + d*x)**5/(5*a**3*d) + 13*csc(c + d*x)**3/(3*a**3*d) - 3*csc(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_100():
    f = cot(c + d*x)**6/(a*sec(c + d*x) + a)**3
    F = -x/a**3 + 4*cot(c + d*x)**11/(11*a**3*d) - cot(c + d*x)**9/(9*a**3*d) + cot(c + d*x)**7/(7*a**3*d) - cot(c + d*x)**5/(5*a**3*d) + cot(c + d*x)**3/(3*a**3*d) - cot(c + d*x)/(a**3*d) - 4*csc(c + d*x)**11/(11*a**3*d) + 19*csc(c + d*x)**9/(9*a**3*d) - 36*csc(c + d*x)**7/(7*a**3*d) + 34*csc(c + d*x)**5/(5*a**3*d) - 16*csc(c + d*x)**3/(3*a**3*d) + 3*csc(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_101():
    f = (e*tan(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)
    F = -sqrt(2)*a*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) + sqrt(2)*a*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) + sqrt(2)*a*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*a*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) + 6*a*e**2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*d*sqrt(sin(2*c + 2*d*x))) - 6*a*e*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*d) + 2*e*(e*tan(c + d*x))**(sympy.S(3)/2)*(3*a*sec(c + d*x) + 5*a)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_102():
    f = (e*tan(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)
    F = sqrt(2)*a*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) - sqrt(2)*a*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) + sqrt(2)*a*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*a*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - a*e**2*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*d*sqrt(e*tan(c + d*x))) + 2*e*sqrt(e*tan(c + d*x))*(a*sec(c + d*x) + 3*a)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_103():
    f = sqrt(e*tan(c + d*x))*(a*sec(c + d*x) + a)
    F = sqrt(2)*a*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) - sqrt(2)*a*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) - sqrt(2)*a*sqrt(e)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) + sqrt(2)*a*sqrt(e)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - 2*a*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*sqrt(sin(2*c + 2*d*x))) + 2*a*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_104():
    f = (a*sec(c + d*x) + a)/sqrt(e*tan(c + d*x))
    F = a*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(d*sqrt(e*tan(c + d*x))) - sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*sqrt(e)) + sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*sqrt(e)) - sqrt(2)*a*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e)) + sqrt(2)*a*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_105():
    f = (a*sec(c + d*x) + a)/(e*tan(c + d*x))**(sympy.S(3)/2)
    F = -2*a*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*e**2*sqrt(sin(2*c + 2*d*x))) + 2*a*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(d*e**3) - sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(3)/2)) + sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(3)/2)) + sqrt(2)*a*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) - sqrt(2)*a*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) - (2*a*sec(c + d*x) + 2*a)/(d*e*sqrt(e*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_106():
    f = (a*sec(c + d*x) + a)/(e*tan(c + d*x))**(sympy.S(5)/2)
    F = -a*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*d*e**2*sqrt(e*tan(c + d*x))) + sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(5)/2)) - sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(5)/2)) + sqrt(2)*a*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)) - sqrt(2)*a*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)) - (2*a*sec(c + d*x) + 2*a)/(3*d*e*(e*tan(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_107():
    f = (a*sec(c + d*x) + a)/(e*tan(c + d*x))**(sympy.S(7)/2)
    F = 6*a*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*d*e**4*sqrt(sin(2*c + 2*d*x))) - 6*a*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*d*e**5) + sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(7)/2)) - sqrt(2)*a*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(7)/2)) - sqrt(2)*a*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2)) + sqrt(2)*a*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2)) - (2*a*sec(c + d*x) + 2*a)/(5*d*e*(e*tan(c + d*x))**(sympy.S(5)/2)) + (6*a*sec(c + d*x) + 10*a)/(5*d*e**3*sqrt(e*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_108():
    f = (e*tan(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)**2
    F = -sqrt(2)*a**2*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) + sqrt(2)*a**2*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) + sqrt(2)*a**2*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*a**2*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) + 12*a**2*e**2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*d*sqrt(sin(2*c + 2*d*x))) - 12*a**2*e*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*d) + 4*a**2*e*(e*tan(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)/(5*d) + 2*a**2*e*(e*tan(c + d*x))**(sympy.S(3)/2)/(3*d) + 2*a**2*(e*tan(c + d*x))**(sympy.S(7)/2)/(7*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_109():
    f = (e*tan(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)**2
    F = sqrt(2)*a**2*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) - sqrt(2)*a**2*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) + sqrt(2)*a**2*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - sqrt(2)*a**2*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - 2*a**2*e**2*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*d*sqrt(e*tan(c + d*x))) + 4*a**2*e*sqrt(e*tan(c + d*x))*sec(c + d*x)/(3*d) + 2*a**2*e*sqrt(e*tan(c + d*x))/d + 2*a**2*(e*tan(c + d*x))**(sympy.S(5)/2)/(5*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_110():
    f = sqrt(e*tan(c + d*x))*(a*sec(c + d*x) + a)**2
    F = sqrt(2)*a**2*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) - sqrt(2)*a**2*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d) - sqrt(2)*a**2*sqrt(e)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) + sqrt(2)*a**2*sqrt(e)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d) - 4*a**2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*sqrt(sin(2*c + 2*d*x))) + 4*a**2*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(d*e) + 2*a**2*(e*tan(c + d*x))**(sympy.S(3)/2)/(3*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_111():
    f = (a*sec(c + d*x) + a)**2/sqrt(e*tan(c + d*x))
    F = 2*a**2*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(d*sqrt(e*tan(c + d*x))) + 2*a**2*sqrt(e*tan(c + d*x))/(d*e) - sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*sqrt(e)) + sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*sqrt(e)) - sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e)) + sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_112():
    f = (a*sec(c + d*x) + a)**2/(e*tan(c + d*x))**(sympy.S(3)/2)
    F = -4*a**2*cos(c + d*x)/(d*e*sqrt(e*tan(c + d*x))) - 4*a**2/(d*e*sqrt(e*tan(c + d*x))) - 4*a**2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*e**2*sqrt(sin(2*c + 2*d*x))) - sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(3)/2)) + sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(3)/2)) + sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2)) - sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_113():
    f = (a*sec(c + d*x) + a)**2/(e*tan(c + d*x))**(sympy.S(5)/2)
    F = -4*a**2*sec(c + d*x)/(3*d*e*(e*tan(c + d*x))**(sympy.S(3)/2)) - 4*a**2/(3*d*e*(e*tan(c + d*x))**(sympy.S(3)/2)) - 2*a**2*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*d*e**2*sqrt(e*tan(c + d*x))) + sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(5)/2)) - sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(5)/2)) + sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2)) - sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_114():
    f = (a*sec(c + d*x) + a)**2/(e*tan(c + d*x))**(sympy.S(7)/2)
    F = -4*a**2*sec(c + d*x)/(5*d*e*(e*tan(c + d*x))**(sympy.S(5)/2)) - 4*a**2/(5*d*e*(e*tan(c + d*x))**(sympy.S(5)/2)) + 12*a**2*cos(c + d*x)/(5*d*e**3*sqrt(e*tan(c + d*x))) + 2*a**2/(d*e**3*sqrt(e*tan(c + d*x))) + 12*a**2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*d*e**4*sqrt(sin(2*c + 2*d*x))) + sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(7)/2)) - sqrt(2)*a**2*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*d*e**(sympy.S(7)/2)) - sqrt(2)*a**2*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2)) + sqrt(2)*a**2*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*d*e**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_115():
    f = (e*tan(c + d*x))**(sympy.S(11)/2)/(a*sec(c + d*x) + a)
    F = sqrt(2)*e**(sympy.S(11)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) - sqrt(2)*e**(sympy.S(11)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) + sqrt(2)*e**(sympy.S(11)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) - sqrt(2)*e**(sympy.S(11)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) + 5*e**6*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(21*a*d*sqrt(e*tan(c + d*x))) + 2*e**5*sqrt(e*tan(c + d*x))*(21 - 5*sec(c + d*x))/(21*a*d) - 2*e**3*(e*tan(c + d*x))**(sympy.S(5)/2)*(7 - 5*sec(c + d*x))/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_116():
    f = (e*tan(c + d*x))**(sympy.S(9)/2)/(a*sec(c + d*x) + a)
    F = sqrt(2)*e**(sympy.S(9)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) - sqrt(2)*e**(sympy.S(9)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) - sqrt(2)*e**(sympy.S(9)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) + sqrt(2)*e**(sympy.S(9)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) + 6*e**4*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*a*d*sqrt(sin(2*c + 2*d*x))) - 2*e**3*(e*tan(c + d*x))**(sympy.S(3)/2)*(5 - 3*sec(c + d*x))/(15*a*d) - 6*e**3*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_117():
    f = (e*tan(c + d*x))**(sympy.S(7)/2)/(a*sec(c + d*x) + a)
    F = -sqrt(2)*e**(sympy.S(7)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) + sqrt(2)*e**(sympy.S(7)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) - sqrt(2)*e**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) + sqrt(2)*e**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) - e**4*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*a*d*sqrt(e*tan(c + d*x))) - 2*e**3*sqrt(e*tan(c + d*x))*(3 - sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_118():
    f = (e*tan(c + d*x))**(sympy.S(5)/2)/(a*sec(c + d*x) + a)
    F = -sqrt(2)*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) + sqrt(2)*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) + sqrt(2)*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) - sqrt(2)*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) - 2*e**2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(a*d*sqrt(sin(2*c + 2*d*x))) + 2*e*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_119():
    f = (e*tan(c + d*x))**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = sqrt(2)*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) - sqrt(2)*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) + sqrt(2)*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) - sqrt(2)*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) + e**2*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(a*d*sqrt(e*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_120():
    f = sqrt(e*tan(c + d*x))/(a*sec(c + d*x) + a)
    F = sqrt(2)*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) - sqrt(2)*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d) - sqrt(2)*sqrt(e)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) + sqrt(2)*sqrt(e)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d) + 2*e*(1 - sec(c + d*x))/(a*d*sqrt(e*tan(c + d*x))) - 2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(a*d*sqrt(sin(2*c + 2*d*x))) + 2*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(a*d*e)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_121():
    f = 1/(sqrt(e*tan(c + d*x))*(a*sec(c + d*x) + a))
    F = 2*e*(1 - sec(c + d*x))/(3*a*d*(e*tan(c + d*x))**(sympy.S(3)/2)) - sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*a*d*sqrt(e*tan(c + d*x))) - sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d*sqrt(e)) + sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d*sqrt(e)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d*sqrt(e)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_122():
    f = 1/((e*tan(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a))
    F = 2*e*(1 - sec(c + d*x))/(5*a*d*(e*tan(c + d*x))**(sympy.S(5)/2)) - (10 - 6*sec(c + d*x))/(5*a*d*e*sqrt(e*tan(c + d*x))) + 6*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*a*d*e**2*sqrt(sin(2*c + 2*d*x))) - 6*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*a*d*e**3) - sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d*e**(sympy.S(3)/2)) + sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d*e**(sympy.S(3)/2)) + sqrt(2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d*e**(sympy.S(3)/2)) - sqrt(2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d*e**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_123():
    f = 1/((e*tan(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a))
    F = 2*e*(1 - sec(c + d*x))/(7*a*d*(e*tan(c + d*x))**(sympy.S(7)/2)) - (14 - 10*sec(c + d*x))/(21*a*d*e*(e*tan(c + d*x))**(sympy.S(3)/2)) + 5*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(21*a*d*e**2*sqrt(e*tan(c + d*x))) + sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d*e**(sympy.S(5)/2)) - sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a*d*e**(sympy.S(5)/2)) + sqrt(2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d*e**(sympy.S(5)/2)) - sqrt(2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a*d*e**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_124():
    f = (e*tan(c + d*x))**(sympy.S(13)/2)/(a*sec(c + d*x) + a)**2
    F = -sqrt(2)*e**(sympy.S(13)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) + sqrt(2)*e**(sympy.S(13)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) + sqrt(2)*e**(sympy.S(13)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) - sqrt(2)*e**(sympy.S(13)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) - 12*e**6*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*a**2*d*sqrt(sin(2*c + 2*d*x))) + 12*e**5*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(5*a**2*d) - 4*e**5*(e*tan(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)/(5*a**2*d) + 2*e**5*(e*tan(c + d*x))**(sympy.S(3)/2)/(3*a**2*d) + 2*e**3*(e*tan(c + d*x))**(sympy.S(7)/2)/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_125():
    f = (e*tan(c + d*x))**(sympy.S(11)/2)/(a*sec(c + d*x) + a)**2
    F = sqrt(2)*e**(sympy.S(11)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) - sqrt(2)*e**(sympy.S(11)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) + sqrt(2)*e**(sympy.S(11)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) - sqrt(2)*e**(sympy.S(11)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + 2*e**6*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*a**2*d*sqrt(e*tan(c + d*x))) - 4*e**5*sqrt(e*tan(c + d*x))*sec(c + d*x)/(3*a**2*d) + 2*e**5*sqrt(e*tan(c + d*x))/(a**2*d) + 2*e**3*(e*tan(c + d*x))**(sympy.S(5)/2)/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_126():
    f = (e*tan(c + d*x))**(sympy.S(9)/2)/(a*sec(c + d*x) + a)**2
    F = sqrt(2)*e**(sympy.S(9)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) - sqrt(2)*e**(sympy.S(9)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) - sqrt(2)*e**(sympy.S(9)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + sqrt(2)*e**(sympy.S(9)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + 4*e**4*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(a**2*d*sqrt(sin(2*c + 2*d*x))) - 4*e**3*(e*tan(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)/(a**2*d) + 2*e**3*(e*tan(c + d*x))**(sympy.S(3)/2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_127():
    f = (e*tan(c + d*x))**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**2
    F = -sqrt(2)*e**(sympy.S(7)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) + sqrt(2)*e**(sympy.S(7)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) - sqrt(2)*e**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + sqrt(2)*e**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) - 2*e**4*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(a**2*d*sqrt(e*tan(c + d*x))) + 2*e**3*sqrt(e*tan(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_128():
    f = (e*tan(c + d*x))**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**2
    F = -sqrt(2)*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) + sqrt(2)*e**(sympy.S(5)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) + sqrt(2)*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) - sqrt(2)*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + 4*e**3*cos(c + d*x)/(a**2*d*sqrt(e*tan(c + d*x))) - 4*e**3/(a**2*d*sqrt(e*tan(c + d*x))) + 4*e**2*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(a**2*d*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_129():
    f = (e*tan(c + d*x))**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**2
    F = sqrt(2)*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) - sqrt(2)*e**(sympy.S(3)/2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) + sqrt(2)*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) - sqrt(2)*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + 4*e**3*sec(c + d*x)/(3*a**2*d*(e*tan(c + d*x))**(sympy.S(3)/2)) - 4*e**3/(3*a**2*d*(e*tan(c + d*x))**(sympy.S(3)/2)) + 2*e**2*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*a**2*d*sqrt(e*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_130():
    f = sqrt(e*tan(c + d*x))/(a*sec(c + d*x) + a)**2
    F = sqrt(2)*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) - sqrt(2)*sqrt(e)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d) - sqrt(2)*sqrt(e)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + sqrt(2)*sqrt(e)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d) + 4*e**3*sec(c + d*x)/(5*a**2*d*(e*tan(c + d*x))**(sympy.S(5)/2)) - 4*e**3/(5*a**2*d*(e*tan(c + d*x))**(sympy.S(5)/2)) - 12*e*cos(c + d*x)/(5*a**2*d*sqrt(e*tan(c + d*x))) + 2*e/(a**2*d*sqrt(e*tan(c + d*x))) - 12*sqrt(e*tan(c + d*x))*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*a**2*d*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_131():
    f = 1/(sqrt(e*tan(c + d*x))*(a*sec(c + d*x) + a)**2)
    F = 4*e**3*sec(c + d*x)/(7*a**2*d*(e*tan(c + d*x))**(sympy.S(7)/2)) - 4*e**3/(7*a**2*d*(e*tan(c + d*x))**(sympy.S(7)/2)) - 20*e*sec(c + d*x)/(21*a**2*d*(e*tan(c + d*x))**(sympy.S(3)/2)) + 2*e/(3*a**2*d*(e*tan(c + d*x))**(sympy.S(3)/2)) - 10*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(21*a**2*d*sqrt(e*tan(c + d*x))) - sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) - sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d*sqrt(e)) + sqrt(2)*log(sqrt(e)*tan(c + d*x) + sqrt(e) + sqrt(2)*sqrt(e*tan(c + d*x)))/(4*a**2*d*sqrt(e)) - sqrt(2)*atan(1 - sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d*sqrt(e)) + sqrt(2)*atan(1 + sqrt(2)*sqrt(e*tan(c + d*x))/sqrt(e))/(2*a**2*d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_132():
    f = sqrt(a*sec(c + d*x) + a)*tan(c + d*x)**5
    F = -2*sqrt(a)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 2*sqrt(a*sec(c + d*x) + a)/d + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*a**2*d) - 6*(a*sec(c + d*x) + a)**(sympy.S(7)/2)/(7*a**3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(9)/2)/(9*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_133():
    f = sqrt(a*sec(c + d*x) + a)*tan(c + d*x)**3
    F = 2*sqrt(a)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 2*sqrt(a*sec(c + d*x) + a)/d - 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_134():
    f = sqrt(a*sec(c + d*x) + a)*tan(c + d*x)
    F = -2*sqrt(a)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 2*sqrt(a*sec(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_135():
    f = sqrt(a*sec(c + d*x) + a)*cot(c + d*x)
    F = 2*sqrt(a)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_136():
    f = sqrt(a*sec(c + d*x) + a)*cot(c + d*x)**3
    F = -2*sqrt(a)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 7*sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(8*d) + a/(4*d*sqrt(a*sec(c + d*x) + a)) + a/(2*d*(1 - sec(c + d*x))*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_137():
    f = sqrt(a*sec(c + d*x) + a)*cot(c + d*x)**5
    F = 2*sqrt(a)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 107*sqrt(2)*sqrt(a)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(128*d) + 43*a**2/(96*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 15*a**2/(16*d*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - a**2/(4*d*(1 - sec(c + d*x))**2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 21*a/(64*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_138():
    f = sqrt(a*sec(c + d*x) + a)*tan(c + d*x)**6
    F = -2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**6*tan(c + d*x)**11/(11*d*(a*sec(c + d*x) + a)**(sympy.S(11)/2)) + 10*a**5*tan(c + d*x)**9/(9*d*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) + 2*a**4*tan(c + d*x)**7/(d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 2*a**3*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 2*a**2*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*a*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_139():
    f = sqrt(a*sec(c + d*x) + a)*tan(c + d*x)**4
    F = 2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**4*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 6*a**3*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 2*a**2*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 2*a*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_140():
    f = sqrt(a*sec(c + d*x) + a)*tan(c + d*x)**2
    F = -2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**2*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*a*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_141():
    f = sqrt(a*sec(c + d*x) + a)*cot(c + d*x)**2
    F = -2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + sqrt(2)*sqrt(a)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(2*d) - sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_142():
    f = sqrt(a*sec(c + d*x) + a)*cot(c + d*x)**4
    F = 2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d - 9*sqrt(2)*sqrt(a)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(16*d) + 7*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(8*d) - (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)*cot(c + d*x)**3*sec(c/2 + d*x/2)**2/(4*a*d) + (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(12*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_143():
    f = sqrt(a*sec(c + d*x) + a)*cot(c + d*x)**6
    F = -2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 151*sqrt(2)*sqrt(a)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(256*d) - 105*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(128*d) - 23*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(192*a*d) - (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2*cot(c + d*x)**5*sec(c/2 + d*x/2)**4/(16*a**2*d) - 17*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)*cot(c + d*x)**5*sec(c/2 + d*x/2)**2/(32*a**2*d) + 87*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5/(160*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_144():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**5
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 2*a*sqrt(a*sec(c + d*x) + a)/d + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(7)/2)/(7*a**2*d) - 2*(a*sec(c + d*x) + a)**(sympy.S(9)/2)/(3*a**3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(11)/2)/(11*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_145():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**3
    F = 2*a**(sympy.S(3)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 2*a*sqrt(a*sec(c + d*x) + a)/d - 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*d) - 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(7)/2)/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_146():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 2*a*sqrt(a*sec(c + d*x) + a)/d + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_147():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)
    F = 2*a**(sympy.S(3)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 2*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_148():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3
    F = -2*a**(sympy.S(3)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 5*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(4*d) + a*sqrt(a*sec(c + d*x) + a)/(2*d*(1 - sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_149():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**5
    F = 2*a**(sympy.S(3)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 71*sqrt(2)*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(64*d) + 7*a**2/(32*d*sqrt(a*sec(c + d*x) + a)) - 13*a**2/(16*d*(1 - sec(c + d*x))*sqrt(a*sec(c + d*x) + a)) - a**2/(4*d*(1 - sec(c + d*x))**2*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_150():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**6
    F = -2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**8*tan(c + d*x)**13/(13*d*(a*sec(c + d*x) + a)**(sympy.S(13)/2)) + 14*a**7*tan(c + d*x)**11/(11*d*(a*sec(c + d*x) + a)**(sympy.S(11)/2)) + 34*a**6*tan(c + d*x)**9/(9*d*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) + 30*a**5*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 2*a**4*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 2*a**3*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*a**2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_151():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**4
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**6*tan(c + d*x)**9/(9*d*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) + 10*a**5*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 14*a**4*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 2*a**3*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 2*a**2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_152():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)**2
    F = -2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**4*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 2*a**3*tan(c + d*x)**3/(d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*a**2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_153():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**2
    F = -2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d - 2*a*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_154():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**4
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d - sqrt(2)*a**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*d) + 3*a*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(2*d) - (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_155():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**6
    F = -2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 11*sqrt(2)*a**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*d) - 21*a*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(16*d) + 5*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(24*d) - (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)*cot(c + d*x)**5*sec(c/2 + d*x/2)**2/(4*a*d) + 3*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5/(20*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_156():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**5
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 2*a**2*sqrt(a*sec(c + d*x) + a)/d + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(7)/2)/(7*a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(9)/2)/(9*a**2*d) - 6*(a*sec(c + d*x) + a)**(sympy.S(11)/2)/(11*a**3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(13)/2)/(13*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_157():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**3
    F = 2*a**(sympy.S(5)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 2*a**2*sqrt(a*sec(c + d*x) + a)/d - 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*d) - 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*d) - 2*(a*sec(c + d*x) + a)**(sympy.S(7)/2)/(7*a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(9)/2)/(9*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_158():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 2*a**2*sqrt(a*sec(c + d*x) + a)/d + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_159():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)
    F = 2*a**(sympy.S(5)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 4*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/d + 2*a**2*sqrt(a*sec(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_160():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**3
    F = -2*a**(sympy.S(5)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d + 3*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(2*d) + a**2*sqrt(a*sec(c + d*x) + a)/(d*(1 - sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_161():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5
    F = 2*a**(sympy.S(5)/2)*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/d - 43*sqrt(2)*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(32*d) - 11*a**2*sqrt(a*sec(c + d*x) + a)/(16*d*(1 - sec(c + d*x))) - a**2*sqrt(a*sec(c + d*x) + a)/(4*d*(1 - sec(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_162():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**6
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**10*tan(c + d*x)**15/(15*d*(a*sec(c + d*x) + a)**(sympy.S(15)/2)) + 18*a**9*tan(c + d*x)**13/(13*d*(a*sec(c + d*x) + a)**(sympy.S(13)/2)) + 62*a**8*tan(c + d*x)**11/(11*d*(a*sec(c + d*x) + a)**(sympy.S(11)/2)) + 98*a**7*tan(c + d*x)**9/(9*d*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) + 62*a**6*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 2*a**5*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 2*a**4*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*a**3*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_163():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**4
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**8*tan(c + d*x)**11/(11*d*(a*sec(c + d*x) + a)**(sympy.S(11)/2)) + 14*a**7*tan(c + d*x)**9/(9*d*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) + 34*a**6*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 6*a**5*tan(c + d*x)**5/(d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 2*a**4*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 2*a**3*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_164():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)**2
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**6*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 2*a**5*tan(c + d*x)**5/(d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 14*a**4*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*a**3*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_165():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**2
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d - 4*a**2*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_166():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**4
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**2*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/d - 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_167():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**6
    F = -2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + sqrt(2)*a**(sympy.S(5)/2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(8*d) - 7*a**2*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(4*d) + a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(2*d) - (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_168():
    f = tan(c + d*x)**5/sqrt(a*sec(c + d*x) + a)
    F = 2*sqrt(a*sec(c + d*x) + a)/(a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*a**2*d) - 6*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*a**3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(7)/2)/(7*a**4*d) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_169():
    f = tan(c + d*x)**3/sqrt(a*sec(c + d*x) + a)
    F = -2*sqrt(a*sec(c + d*x) + a)/(a*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*a**2*d) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_170():
    f = tan(c + d*x)/sqrt(a*sec(c + d*x) + a)
    F = -2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_171():
    f = cot(c + d*x)/sqrt(a*sec(c + d*x) + a)
    F = -1/(d*sqrt(a*sec(c + d*x) + a)) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_172():
    f = cot(c + d*x)**3/sqrt(a*sec(c + d*x) + a)
    F = -a/(12*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + a/(2*d*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 7/(8*d*sqrt(a*sec(c + d*x) + a)) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(sqrt(a)*d) + 9*sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(16*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_173():
    f = cot(c + d*x)**5/sqrt(a*sec(c + d*x) + a)
    F = 87*a**2/(160*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 17*a**2/(16*d*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - a**2/(4*d*(1 - sec(c + d*x))**2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 23*a/(192*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 105/(128*d*sqrt(a*sec(c + d*x) + a)) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(sqrt(a)*d) - 151*sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(256*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_174():
    f = tan(c + d*x)**6/sqrt(a*sec(c + d*x) + a)
    F = 2*a**4*tan(c + d*x)**9/(9*d*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) + 6*a**3*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 2*a**2*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 2*a*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_175():
    f = tan(c + d*x)**4/sqrt(a*sec(c + d*x) + a)
    F = 2*a**2*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 2*a*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_176():
    f = tan(c + d*x)**2/sqrt(a*sec(c + d*x) + a)
    F = 2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_177():
    f = cot(c + d*x)**2/sqrt(a*sec(c + d*x) + a)
    F = -sqrt(a*sec(c + d*x) + a)*cos(c + d*x)*cot(c + d*x)*sec(c/2 + d*x/2)**2/(4*a*d) - sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(4*a*d) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) + 7*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(8*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_178():
    f = cot(c + d*x)**4/sqrt(a*sec(c + d*x) + a)
    F = 21*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(64*a*d) - (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2*cot(c + d*x)**3*sec(c/2 + d*x/2)**4/(16*a**2*d) - 15*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)*cot(c + d*x)**3*sec(c/2 + d*x/2)**2/(32*a**2*d) + 43*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(96*a**2*d) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) - 107*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(128*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_179():
    f = cot(c + d*x)**6/sqrt(a*sec(c + d*x) + a)
    F = -189*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(512*a*d) - 323*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(768*a**2*d) - (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3*cot(c + d*x)**5*sec(c/2 + d*x/2)**6/(48*a**3*d) - 23*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2*cot(c + d*x)**5*sec(c/2 + d*x/2)**4/(192*a**3*d) - 101*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)*cot(c + d*x)**5*sec(c/2 + d*x/2)**2/(128*a**3*d) + 579*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5/(640*a**3*d) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) + 835*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(1024*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_180():
    f = tan(c + d*x)**5/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*sqrt(a*sec(c + d*x) + a)/(a**2*d) - 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(a**3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)/(5*a**4*d) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_181():
    f = tan(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*sqrt(a*sec(c + d*x) + a)/(a**2*d) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_182():
    f = tan(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2/(a*d*sqrt(a*sec(c + d*x) + a)) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_183():
    f = cot(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -1/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 3/(2*a*d*sqrt(a*sec(c + d*x) + a)) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_184():
    f = cot(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -3*a/(20*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + a/(2*d*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 5/(24*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 21/(16*a*d*sqrt(a*sec(c + d*x) + a)) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d) + 11*sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_185():
    f = cot(c + d*x)**5/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 139*a**2/(224*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) - 19*a**2/(16*d*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) - a**2/(4*d*(1 - sec(c + d*x))**2*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 15*a/(64*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 53/(384*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 309/(256*a*d*sqrt(a*sec(c + d*x) + a)) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*d) - 203*sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(512*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_186():
    f = tan(c + d*x)**6/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*a**2*tan(c + d*x)**7/(7*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 2*a*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 2*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*tan(c + d*x)/(a*d*sqrt(a*sec(c + d*x) + a)) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_187():
    f = tan(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*tan(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 2*tan(c + d*x)/(a*d*sqrt(a*sec(c + d*x) + a)) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_188():
    f = tan(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 2*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_189():
    f = cot(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**2*cot(c + d*x)*sec(c/2 + d*x/2)**4/(16*a**2*d) - 13*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)*cot(c + d*x)*sec(c/2 + d*x/2)**2/(32*a**2*d) + 7*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(32*a**2*d) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 71*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(64*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_190():
    f = cot(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -21*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(256*a**2*d) - (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3*cot(c + d*x)**3*sec(c/2 + d*x/2)**6/(48*a**3*d) - 7*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2*cot(c + d*x)**3*sec(c/2 + d*x/2)**4/(64*a**3*d) - 81*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)*cot(c + d*x)**3*sec(c/2 + d*x/2)**2/(128*a**3*d) + 277*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(384*a**3*d) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 533*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(512*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_191():
    f = cot(c + d*x)**6/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -21*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(8192*a**2*d) - 8171*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(12288*a**3*d) - (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**4*cot(c + d*x)**5*sec(c/2 + d*x/2)**8/(128*a**4*d) - 29*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3*cot(c + d*x)**5*sec(c/2 + d*x/2)**6/(768*a**4*d) - 511*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2*cot(c + d*x)**5*sec(c/2 + d*x/2)**4/(3072*a**4*d) - 2045*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)*cot(c + d*x)**5*sec(c/2 + d*x/2)**2/(2048*a**4*d) + 12267*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5/(10240*a**4*d) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 16363*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(16384*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_192():
    f = tan(c + d*x)**5/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -6*sqrt(a*sec(c + d*x) + a)/(a**3*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)/(3*a**4*d) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_193():
    f = tan(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -4/(a**2*d*sqrt(a*sec(c + d*x) + a)) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_194():
    f = tan(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 2/(3*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2/(a**2*d*sqrt(a*sec(c + d*x) + a)) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_195():
    f = cot(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -1/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 1/(2*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 7/(4*a**2*d*sqrt(a*sec(c + d*x) + a)) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(8*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_196():
    f = cot(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -5*a/(28*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + a/(2*d*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 3/(40*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 19/(48*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 51/(32*a**2*d*sqrt(a*sec(c + d*x) + a)) - 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d) + 13*sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(64*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_197():
    f = cot(c + d*x)**5/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 199*a**2/(288*d*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) - 21*a**2/(16*d*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) - a**2/(4*d*(1 - sec(c + d*x))**2*(a*sec(c + d*x) + a)**(sympy.S(9)/2)) + 135*a/(448*d*(a*sec(c + d*x) + a)**(sympy.S(7)/2)) + 7/(640*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 83/(256*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 761/(512*a**2*d*sqrt(a*sec(c + d*x) + a)) + 2*atanh(sqrt(a*sec(c + d*x) + a)/sqrt(a))/(a**(sympy.S(5)/2)*d) - 263*sqrt(2)*atanh(sqrt(2)*sqrt(a*sec(c + d*x) + a)/(2*sqrt(a)))/(1024*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_198():
    f = tan(c + d*x)**6/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*tan(c + d*x)**5/(5*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 2*tan(c + d*x)**3/(3*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*tan(c + d*x)/(a**2*d*sqrt(a*sec(c + d*x) + a)) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_199():
    f = tan(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*tan(c + d*x)/(a**2*d*sqrt(a*sec(c + d*x) + a)) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 4*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_200():
    f = tan(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = sin(c + d*x)*sec(c/2 + d*x/2)**2/(2*a**2*d*sqrt(a*sec(c + d*x) + a)) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 3*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_201():
    f = cot(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**3*cot(c + d*x)*sec(c/2 + d*x/2)**6/(48*a**3*d) - 19*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**2*cot(c + d*x)*sec(c/2 + d*x/2)**4/(192*a**3*d) - 191*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)*cot(c + d*x)*sec(c/2 + d*x/2)**2/(384*a**3*d) + 63*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(128*a**3*d) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 319*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(256*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_202():
    f = cot(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -1491*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(4096*a**3*d) - (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**4*cot(c + d*x)**3*sec(c/2 + d*x/2)**8/(128*a**4*d) - 9*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3*cot(c + d*x)**3*sec(c/2 + d*x/2)**6/(256*a**4*d) - 145*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2*cot(c + d*x)**3*sec(c/2 + d*x/2)**4/(1024*a**4*d) - 1527*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)*cot(c + d*x)**3*sec(c/2 + d*x/2)**2/(2048*a**4*d) + 5587*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(6144*a**4*d) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 9683*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(8192*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_203():
    f = cot(c + d*x)**6/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 8925*sqrt(a*sec(c + d*x) + a)*cot(c + d*x)/(32768*a**3*d) - 41693*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cot(c + d*x)**3/(49152*a**4*d) - (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5*cot(c + d*x)**5*sec(c/2 + d*x/2)**10/(320*a**5*d) - 7*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**4*cot(c + d*x)**5*sec(c/2 + d*x/2)**8/(512*a**5*d) - 155*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3*cot(c + d*x)**5*sec(c/2 + d*x/2)**6/(3072*a**5*d) - 2473*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2*cot(c + d*x)**5*sec(c/2 + d*x/2)**4/(12288*a**5*d) - 9467*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)*cot(c + d*x)**5*sec(c/2 + d*x/2)**2/(8192*a**5*d) + 58077*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cot(c + d*x)**5/(40960*a**5*d) - 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 74461*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(65536*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_204():
    f = tan(e + f*x)**2/(a*sec(e + f*x) + a)**(sympy.S(9)/2)
    F = tan(e + f*x)/(3*a*f*(a*sec(e + f*x) + a)**(sympy.S(7)/2)) + 11*tan(e + f*x)/(24*a**2*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 27*tan(e + f*x)/(32*a**3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(9)/2)*f) + 91*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(64*a**(sympy.S(9)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_205():
    f = (e*tan(c + d*x))**m*(a*sec(c + d*x) + a)**n
    F = 2**(m + n + 1)*(e*tan(c + d*x))**(m + 1)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(m + n + 1)*appellf1(m/2 + sympy.S.Half, 1, m + n, m/2 + sympy.S(3)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_206():
    f = (e*tan(c + d*x))**m*(a*sec(c + d*x) + a)**3
    F = 3*a**3*(e*tan(c + d*x))**(m + 1)*(cos(c + d*x)**2)**(m/2 + 1)*hyper((m/2 + sympy.S.Half, m/2 + 1), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*e*(m + 1)) + a**3*(e*tan(c + d*x))**(m + 1)*(cos(c + d*x)**2)**(m/2 + 2)*hyper((m/2 + sympy.S.Half, m/2 + 2), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)*sec(c + d*x)**3/(d*e*(m + 1)) + a**3*(e*tan(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*e*(m + 1)) + 3*a**3*(e*tan(c + d*x))**(m + 1)/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_207():
    f = (e*tan(c + d*x))**m*(a*sec(c + d*x) + a)**2
    F = 2*a**2*(e*tan(c + d*x))**(m + 1)*(cos(c + d*x)**2)**(m/2 + 1)*hyper((m/2 + sympy.S.Half, m/2 + 1), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*e*(m + 1)) + a**2*(e*tan(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*e*(m + 1)) + a**2*(e*tan(c + d*x))**(m + 1)/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_208():
    f = (e*tan(c + d*x))**m*(a*sec(c + d*x) + a)
    F = a*(e*tan(c + d*x))**(m + 1)*(cos(c + d*x)**2)**(m/2 + 1)*hyper((m/2 + sympy.S.Half, m/2 + 1), (m/2 + sympy.S(3)/2,), sin(c + d*x)**2)*sec(c + d*x)/(d*e*(m + 1)) + a*(e*tan(c + d*x))**(m + 1)*hyper((1, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), -tan(c + d*x)**2)/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_209():
    f = (e*tan(c + d*x))**m/(a*sec(c + d*x) + a)
    F = -e*(e*tan(c + d*x))**(m - 1)*(cos(c + d*x)**2)**(m/2)*hyper((m/2, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), sin(c + d*x)**2)*sec(c + d*x)/(a*d*(1 - m)) + e*(e*tan(c + d*x))**(m - 1)*hyper((1, m/2 + sympy.S(-1)/2), (m/2 + sympy.S.Half,), -tan(c + d*x)**2)/(a*d*(1 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_210():
    f = (e*tan(c + d*x))**m/(a*sec(c + d*x) + a)**2
    F = 2*e**3*(e*tan(c + d*x))**(m - 3)*(cos(c + d*x)**2)**(m/2 - 1)*hyper((m/2 + sympy.S(-3)/2, m/2 - 1), (m/2 + sympy.S(-1)/2,), sin(c + d*x)**2)*sec(c + d*x)/(a**2*d*(3 - m)) - e**3*(e*tan(c + d*x))**(m - 3)*hyper((1, m/2 + sympy.S(-3)/2), (m/2 + sympy.S(-1)/2,), -tan(c + d*x)**2)/(a**2*d*(3 - m)) - e**3*(e*tan(c + d*x))**(m - 3)/(a**2*d*(3 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_211():
    f = (e*tan(c + d*x))**m/(a*sec(c + d*x) + a)**3
    F = -3*e**5*(e*tan(c + d*x))**(m - 5)*(cos(c + d*x)**2)**(m/2 - 2)*hyper((m/2 + sympy.S(-5)/2, m/2 - 2), (m/2 + sympy.S(-3)/2,), sin(c + d*x)**2)*sec(c + d*x)/(a**3*d*(5 - m)) - e**5*(e*tan(c + d*x))**(m - 5)*(cos(c + d*x)**2)**(m/2 - 1)*hyper((m/2 + sympy.S(-5)/2, m/2 - 1), (m/2 + sympy.S(-3)/2,), sin(c + d*x)**2)*sec(c + d*x)**3/(a**3*d*(5 - m)) + e**5*(e*tan(c + d*x))**(m - 5)*hyper((1, m/2 + sympy.S(-5)/2), (m/2 + sympy.S(-3)/2,), -tan(c + d*x)**2)/(a**3*d*(5 - m)) + 3*e**5*(e*tan(c + d*x))**(m - 5)/(a**3*d*(5 - m))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_212():
    f = (e*tan(c + d*x))**m*(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2**(m + sympy.S(5)/2)*(e*tan(c + d*x))**(m + 1)*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*(1/(sec(c + d*x) + 1))**(m + sympy.S(5)/2)*appellf1(m/2 + sympy.S.Half, 1, m + sympy.S(3)/2, m/2 + sympy.S(3)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_213():
    f = (e*tan(c + d*x))**m*sqrt(a*sec(c + d*x) + a)
    F = 2**(m + sympy.S(3)/2)*(e*tan(c + d*x))**(m + 1)*sqrt(a*sec(c + d*x) + a)*(1/(sec(c + d*x) + 1))**(m + sympy.S(3)/2)*appellf1(m/2 + sympy.S.Half, 1, m + sympy.S.Half, m/2 + sympy.S(3)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(d*e*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_214():
    f = (e*tan(c + d*x))**m/sqrt(a*sec(c + d*x) + a)
    F = 2**(m + sympy.S.Half)*(e*tan(c + d*x))**(m + 1)*(1/(sec(c + d*x) + 1))**(m + sympy.S.Half)*appellf1(m/2 + sympy.S.Half, 1, m + sympy.S(-1)/2, m/2 + sympy.S(3)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(d*e*(m + 1)*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_215():
    f = (e*tan(c + d*x))**m/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2**(m + sympy.S(-1)/2)*(e*tan(c + d*x))**(m + 1)*(1/(sec(c + d*x) + 1))**(m + sympy.S(-1)/2)*appellf1(m/2 + sympy.S.Half, 1, m + sympy.S(-3)/2, m/2 + sympy.S(3)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(d*e*(m + 1)*(a*sec(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_216():
    f = (a*sec(c + d*x) + a)**n*tan(c + d*x)**7
    F = (a*sec(c + d*x) + a)**(n + 4)*hyper((1, n + 4), (n + 5,), sec(c + d*x) + 1)/(a**4*d*(n + 4)) + 7*(a*sec(c + d*x) + a)**(n + 4)/(a**4*d*(n + 4)) - 5*(a*sec(c + d*x) + a)**(n + 5)/(a**5*d*(n + 5)) + (a*sec(c + d*x) + a)**(n + 6)/(a**6*d*(n + 6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_217():
    f = (a*sec(c + d*x) + a)**n*tan(c + d*x)**5
    F = -(a*sec(c + d*x) + a)**(n + 3)*hyper((1, n + 3), (n + 4,), sec(c + d*x) + 1)/(a**3*d*(n + 3)) - 3*(a*sec(c + d*x) + a)**(n + 3)/(a**3*d*(n + 3)) + (a*sec(c + d*x) + a)**(n + 4)/(a**4*d*(n + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_218():
    f = (a*sec(c + d*x) + a)**n*tan(c + d*x)**3
    F = (a*sec(c + d*x) + a)**(n + 2)*hyper((1, n + 2), (n + 3,), sec(c + d*x) + 1)/(a**2*d*(n + 2)) + (a*sec(c + d*x) + a)**(n + 2)/(a**2*d*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_219():
    f = (a*sec(c + d*x) + a)**n*tan(c + d*x)
    F = -(a*sec(c + d*x) + a)**(n + 1)*hyper((1, n + 1), (n + 2,), sec(c + d*x) + 1)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_220():
    f = (a*sec(c + d*x) + a)**n*cot(c + d*x)
    F = -(a*sec(c + d*x) + a)**n*hyper((1, n), (n + 1,), sec(c + d*x)/2 + sympy.S.Half)/(2*d*n) + (a*sec(c + d*x) + a)**n*hyper((1, n), (n + 1,), sec(c + d*x) + 1)/(d*n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_221():
    f = (a*sec(c + d*x) + a)**n*cot(c + d*x)**3
    F = a*(a*sec(c + d*x) + a)**(n - 1)/(2*d*(1 - sec(c + d*x))) - a*(4 - n)*(a*sec(c + d*x) + a)**(n - 1)*hyper((1, n - 1), (n,), sec(c + d*x)/2 + sympy.S.Half)/(4*d*(1 - n)) + a*(a*sec(c + d*x) + a)**(n - 1)*hyper((1, n - 1), (n,), sec(c + d*x) + 1)/(d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_222():
    f = (a*sec(c + d*x) + a)**n*tan(c + d*x)**4
    F = 2**(n + 5)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n + 5)*tan(c + d*x)**5*appellf1(sympy.S(5)/2, 1, n + 4, sympy.S(7)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_223():
    f = (a*sec(c + d*x) + a)**n*tan(c + d*x)**2
    F = 2**(n + 3)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n + 3)*tan(c + d*x)**3*appellf1(sympy.S(3)/2, 1, n + 2, sympy.S(5)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_224():
    f = (a*sec(c + d*x) + a)**n*cot(c + d*x)**2
    F = -2**(n - 1)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n - 1)*cot(c + d*x)*appellf1(sympy.S(-1)/2, 1, n - 2, sympy.S.Half, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_225():
    f = (a*sec(c + d*x) + a)**n*cot(c + d*x)**4
    F = -2**(n - 3)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n - 3)*cot(c + d*x)**3*appellf1(sympy.S(-3)/2, 1, n - 4, sympy.S(-1)/2, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_226():
    f = (a*sec(c + d*x) + a)**n*tan(c + d*x)**(sympy.S(3)/2)
    F = 2**(n + sympy.S(7)/2)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n + sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)*appellf1(sympy.S(5)/4, 1, n + sympy.S(3)/2, sympy.S(9)/4, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_227():
    f = (a*sec(c + d*x) + a)**n*sqrt(tan(c + d*x))
    F = 2**(n + sympy.S(5)/2)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n + sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)*appellf1(sympy.S(3)/4, 1, n + sympy.S.Half, sympy.S(7)/4, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_228():
    f = (a*sec(c + d*x) + a)**n/sqrt(tan(c + d*x))
    F = 2**(n + sympy.S(3)/2)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n + sympy.S.Half)*sqrt(tan(c + d*x))*appellf1(sympy.S(1)/4, 1, n + sympy.S(-1)/2, sympy.S(5)/4, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_229():
    f = (a*sec(c + d*x) + a)**n/tan(c + d*x)**(sympy.S(3)/2)
    F = -2**(n + sympy.S.Half)*(a*sec(c + d*x) + a)**n*(1/(sec(c + d*x) + 1))**(n + sympy.S(-1)/2)*appellf1(sympy.S(-1)/4, 1, n + sympy.S(-3)/2, sympy.S(3)/4, (-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a), -(-a*sec(c + d*x) + a)/(a*sec(c + d*x) + a))/(d*sqrt(tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_230():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)
    F = sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(5)/2)/(4*d) - sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(5)/2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(5)/2)/(4*d) - a*(e*cot(c + d*x))**(sympy.S(5)/2)*sqrt(sin(2*c + 2*d*x))*tan(c + d*x)**2*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*d) - sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - 2*(e*cot(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_231():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)
    F = -sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*d) + sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(3)/2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*d) + 2*a*(e*cot(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*tan(c + d*x)**2/d - 2*a*(e*cot(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*tan(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*sqrt(sin(2*c + 2*d*x))) - sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*a*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - 2*(e*cot(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_232():
    f = sqrt(e*cot(c + d*x))*(a*sec(c + d*x) + a)
    F = -sqrt(2)*a*sqrt(e*cot(c + d*x))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*sqrt(tan(c + d*x))/(4*d) + sqrt(2)*a*sqrt(e*cot(c + d*x))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*sqrt(tan(c + d*x))/(4*d) + a*sqrt(e*cot(c + d*x))*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/d + sqrt(2)*a*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*a*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_233():
    f = (a*sec(c + d*x) + a)/sqrt(e*cot(c + d*x))
    F = sqrt(2)*a*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - sqrt(2)*a*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) + 2*a*sin(c + d*x)/(d*sqrt(e*cot(c + d*x))) + sqrt(2)*a*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) + sqrt(2)*a*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - 2*a*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*sqrt(e*cot(c + d*x))*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_234():
    f = (a*sec(c + d*x) + a)/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = sqrt(2)*a*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*a*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - a*sqrt(sin(2*c + 2*d*x))*cot(c + d*x)*csc(c + d*x)*elliptic_f(c + d*x - pi/4, 2)/(3*d*(e*cot(c + d*x))**(sympy.S(3)/2)) - sqrt(2)*a*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*a*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + 2*(a*sec(c + d*x) + 3*a)*cot(c + d*x)/(3*d*(e*cot(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_235():
    f = (e*cot(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)**2
    F = sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(5)/2)/(4*d) - sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(5)/2)/(4*d) - 2*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)*sqrt(sin(2*c + 2*d*x))*tan(c + d*x)**2*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*d) - sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - 4*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)*sec(c + d*x)/(3*d) - 4*a**2*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_236():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)**2
    F = -sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*d) + sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*d) - 4*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/d - 4*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*tan(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*sqrt(sin(2*c + 2*d*x))) - sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) - sqrt(2)*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) - 4*a**2*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_237():
    f = sqrt(e*cot(c + d*x))*(a*sec(c + d*x) + a)**2
    F = -sqrt(2)*a**2*sqrt(e*cot(c + d*x))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*sqrt(tan(c + d*x))/(4*d) + sqrt(2)*a**2*sqrt(e*cot(c + d*x))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*sqrt(tan(c + d*x))/(4*d) + 2*a**2*sqrt(e*cot(c + d*x))*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/d + sqrt(2)*a**2*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d) + sqrt(2)*a**2*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d) + 2*a**2*sqrt(e*cot(c + d*x))*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_238():
    f = (a*sec(c + d*x) + a)**2/sqrt(e*cot(c + d*x))
    F = sqrt(2)*a**2*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - sqrt(2)*a**2*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) + 4*a**2*sin(c + d*x)/(d*sqrt(e*cot(c + d*x))) + 2*a**2*tan(c + d*x)/(3*d*sqrt(e*cot(c + d*x))) + sqrt(2)*a**2*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) + sqrt(2)*a**2*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - 4*a**2*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(d*sqrt(e*cot(c + d*x))*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_239():
    f = (a*sec(c + d*x) + a)**2/(e*cot(c + d*x))**(sympy.S(3)/2)
    F = sqrt(2)*a**2*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*a**2*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - 2*a**2*sqrt(sin(2*c + 2*d*x))*cot(c + d*x)*csc(c + d*x)*elliptic_f(c + d*x - pi/4, 2)/(3*d*(e*cot(c + d*x))**(sympy.S(3)/2)) + 2*a**2*tan(c + d*x)/(5*d*(e*cot(c + d*x))**(sympy.S(3)/2)) + 2*a**2*cot(c + d*x)/(d*(e*cot(c + d*x))**(sympy.S(3)/2)) + 4*a**2*csc(c + d*x)/(3*d*(e*cot(c + d*x))**(sympy.S(3)/2)) - sqrt(2)*a**2*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*a**2*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_240():
    f = (e*cot(c + d*x))**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = 2*(e*cot(c + d*x))**(sympy.S(3)/2)*(1 - sec(c + d*x))*cot(c + d*x)/(5*a*d) - 2*(e*cot(c + d*x))**(sympy.S(3)/2)*(5 - 3*sec(c + d*x))*tan(c + d*x)/(5*a*d) - sqrt(2)*(e*cot(c + d*x))**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*a*d) + sqrt(2)*(e*cot(c + d*x))**(sympy.S(3)/2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*tan(c + d*x)**(sympy.S(3)/2)/(4*a*d) - 6*(e*cot(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*tan(c + d*x)**2/(5*a*d) + 6*(e*cot(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*tan(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*a*d*sqrt(sin(2*c + 2*d*x))) - sqrt(2)*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d) - sqrt(2)*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_241():
    f = sqrt(e*cot(c + d*x))/(a*sec(c + d*x) + a)
    F = 2*sqrt(e*cot(c + d*x))*(1 - sec(c + d*x))*cot(c + d*x)/(3*a*d) - sqrt(2)*sqrt(e*cot(c + d*x))*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*sqrt(tan(c + d*x))/(4*a*d) + sqrt(2)*sqrt(e*cot(c + d*x))*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)*sqrt(tan(c + d*x))/(4*a*d) - sqrt(e*cot(c + d*x))*sqrt(sin(2*c + 2*d*x))*elliptic_f(c + d*x - pi/4, 2)*sec(c + d*x)/(3*a*d) + sqrt(2)*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d) + sqrt(2)*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_242():
    f = 1/(sqrt(e*cot(c + d*x))*(a*sec(c + d*x) + a))
    F = 2*(1 - sec(c + d*x))*cot(c + d*x)/(a*d*sqrt(e*cot(c + d*x))) + sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) + 2*sin(c + d*x)/(a*d*sqrt(e*cot(c + d*x))) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - 2*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(a*d*sqrt(e*cot(c + d*x))*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_243():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a))
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + sqrt(sin(2*c + 2*d*x))*cot(c + d*x)*csc(c + d*x)*elliptic_f(c + d*x - pi/4, 2)/(a*d*(e*cot(c + d*x))**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_244():
    f = 1/((e*cot(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a))
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) + 2*cos(c + d*x)*cot(c + d*x)/(a*d*(e*cot(c + d*x))**(sympy.S(5)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) - 2*cos(c + d*x)*cot(c + d*x)**2*elliptic_e(c + d*x - pi/4, 2)/(a*d*(e*cot(c + d*x))**(sympy.S(5)/2)*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_245():
    f = 1/((e*cot(c + d*x))**(sympy.S(7)/2)*(a*sec(c + d*x) + a))
    F = -2*(3 - sec(c + d*x))*cot(c + d*x)**3/(3*a*d*(e*cot(c + d*x))**(sympy.S(7)/2)) - sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2)) + sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2)) - sqrt(sin(2*c + 2*d*x))*cot(c + d*x)**3*csc(c + d*x)*elliptic_f(c + d*x - pi/4, 2)/(3*a*d*(e*cot(c + d*x))**(sympy.S(7)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_246():
    f = 1/((e*cot(c + d*x))**(sympy.S(9)/2)*(a*sec(c + d*x) + a))
    F = -2*(5 - 3*sec(c + d*x))*cot(c + d*x)**3/(15*a*d*(e*cot(c + d*x))**(sympy.S(9)/2)) + sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) - 6*cos(c + d*x)*cot(c + d*x)**3/(5*a*d*(e*cot(c + d*x))**(sympy.S(9)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) + 6*cos(c + d*x)*cot(c + d*x)**4*elliptic_e(c + d*x - pi/4, 2)/(5*a*d*(e*cot(c + d*x))**(sympy.S(9)/2)*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_247():
    f = 1/(sqrt(e*cot(c + d*x))*(a*sec(c + d*x) + a)**2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - 12*cos(c + d*x)*cot(c + d*x)/(5*a**2*d*sqrt(e*cot(c + d*x))) - 4*cot(c + d*x)**3/(5*a**2*d*sqrt(e*cot(c + d*x))) + 4*cot(c + d*x)**2*csc(c + d*x)/(5*a**2*d*sqrt(e*cot(c + d*x))) + 2*cot(c + d*x)/(a**2*d*sqrt(e*cot(c + d*x))) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d*sqrt(e*cot(c + d*x))*sqrt(tan(c + d*x))) - 12*cos(c + d*x)*elliptic_e(c + d*x - pi/4, 2)/(5*a**2*d*sqrt(e*cot(c + d*x))*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_248():
    f = 1/((e*cot(c + d*x))**(sympy.S(3)/2)*(a*sec(c + d*x) + a)**2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) + 2*sqrt(sin(2*c + 2*d*x))*cot(c + d*x)*csc(c + d*x)*elliptic_f(c + d*x - pi/4, 2)/(3*a**2*d*(e*cot(c + d*x))**(sympy.S(3)/2)) - 4*cot(c + d*x)**3/(3*a**2*d*(e*cot(c + d*x))**(sympy.S(3)/2)) + 4*cot(c + d*x)**2*csc(c + d*x)/(3*a**2*d*(e*cot(c + d*x))**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_249():
    f = 1/((e*cot(c + d*x))**(sympy.S(5)/2)*(a*sec(c + d*x) + a)**2)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) + sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) + 4*cos(c + d*x)*cot(c + d*x)**3/(a**2*d*(e*cot(c + d*x))**(sympy.S(5)/2)) - 4*cot(c + d*x)**3/(a**2*d*(e*cot(c + d*x))**(sympy.S(5)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)**(sympy.S(5)/2)) + 4*cos(c + d*x)*cot(c + d*x)**2*elliptic_e(c + d*x - pi/4, 2)/(a**2*d*(e*cot(c + d*x))**(sympy.S(5)/2)*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_250():
    f = 1/((e*cot(c + d*x))**(sympy.S(7)/2)*(a*sec(c + d*x) + a)**2)
    F = -sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2)) + sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2)) - 2*sqrt(sin(2*c + 2*d*x))*cot(c + d*x)**3*csc(c + d*x)*elliptic_f(c + d*x - pi/4, 2)/(a**2*d*(e*cot(c + d*x))**(sympy.S(7)/2)) + 2*cot(c + d*x)**3/(a**2*d*(e*cot(c + d*x))**(sympy.S(7)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_251():
    f = 1/((e*cot(c + d*x))**(sympy.S(9)/2)*(a*sec(c + d*x) + a)**2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) - 4*cos(c + d*x)*cot(c + d*x)**3/(a**2*d*(e*cot(c + d*x))**(sympy.S(9)/2)) + 2*cot(c + d*x)**3/(3*a**2*d*(e*cot(c + d*x))**(sympy.S(9)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(9)/2)*tan(c + d*x)**(sympy.S(9)/2)) + 4*cos(c + d*x)*cot(c + d*x)**4*elliptic_e(c + d*x - pi/4, 2)/(a**2*d*(e*cot(c + d*x))**(sympy.S(9)/2)*sqrt(sin(2*c + 2*d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_252():
    f = 1/((e*cot(c + d*x))**(sympy.S(11)/2)*(a*sec(c + d*x) + a)**2)
    F = sqrt(2)*log(-sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)*tan(c + d*x)**(sympy.S(11)/2)) - sqrt(2)*log(sqrt(2)*sqrt(tan(c + d*x)) + tan(c + d*x) + 1)/(4*a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)*tan(c + d*x)**(sympy.S(11)/2)) + 2*sqrt(sin(2*c + 2*d*x))*cot(c + d*x)**5*csc(c + d*x)*elliptic_f(c + d*x - pi/4, 2)/(3*a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)) + 2*cot(c + d*x)**5/(a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)) - 4*cot(c + d*x)**4*csc(c + d*x)/(3*a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)) + 2*cot(c + d*x)**3/(5*a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) - 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)*tan(c + d*x)**(sympy.S(11)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(tan(c + d*x)) + 1)/(2*a**2*d*(e*cot(c + d*x))**(sympy.S(11)/2)*tan(c + d*x)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_253():
    f = (a + b*sec(c + d*x))*tan(c + d*x)**7
    F = a*log(cos(c + d*x))/d - 16*b*sec(c + d*x)/(35*d) + (7*a + 6*b*sec(c + d*x))*tan(c + d*x)**6/(42*d) + (35*a + 16*b*sec(c + d*x))*tan(c + d*x)**2/(70*d) - (35*a + 24*b*sec(c + d*x))*tan(c + d*x)**4/(140*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_254():
    f = (a + b*sec(c + d*x))*tan(c + d*x)**5
    F = -a*log(cos(c + d*x))/d + 8*b*sec(c + d*x)/(15*d) + (5*a + 4*b*sec(c + d*x))*tan(c + d*x)**4/(20*d) - (15*a + 8*b*sec(c + d*x))*tan(c + d*x)**2/(30*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_255():
    f = (a + b*sec(c + d*x))*tan(c + d*x)**3
    F = a*log(cos(c + d*x))/d - 2*b*sec(c + d*x)/(3*d) + (3*a + 2*b*sec(c + d*x))*tan(c + d*x)**2/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_256():
    f = (a + b*sec(c + d*x))*tan(c + d*x)
    F = -a*log(cos(c + d*x))/d + b*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_257():
    f = (a + b*sec(c + d*x))*cot(c + d*x)
    F = (a - b)*log(cos(c + d*x) + 1)/(2*d) + (a + b)*log(1 - cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_258():
    f = (a + b*sec(c + d*x))*cot(c + d*x)**3
    F = -(a + b*sec(c + d*x))*cot(c + d*x)**2/(2*d) - (2*a - b)*log(cos(c + d*x) + 1)/(4*d) - (2*a + b)*log(1 - cos(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_259():
    f = (a + b*sec(c + d*x))*cot(c + d*x)**5
    F = -(a + b*sec(c + d*x))*cot(c + d*x)**4/(4*d) + (4*a + 3*b*sec(c + d*x))*cot(c + d*x)**2/(8*d) + (8*a - 3*b)*log(cos(c + d*x) + 1)/(16*d) + (8*a + 3*b)*log(1 - cos(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_260():
    f = (a + b*sec(c + d*x))*cot(c + d*x)**7
    F = -(a + b*sec(c + d*x))*cot(c + d*x)**6/(6*d) + (6*a + 5*b*sec(c + d*x))*cot(c + d*x)**4/(24*d) - (8*a + 5*b*sec(c + d*x))*cot(c + d*x)**2/(16*d) - (16*a - 5*b)*log(cos(c + d*x) + 1)/(32*d) - (16*a + 5*b)*log(1 - cos(c + d*x))/(32*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_261():
    f = (a + b*sec(c + d*x))*tan(c + d*x)**6
    F = -a*x - 5*b*atanh(sin(c + d*x))/(16*d) + (6*a + 5*b*sec(c + d*x))*tan(c + d*x)**5/(30*d) - (8*a + 5*b*sec(c + d*x))*tan(c + d*x)**3/(24*d) + (16*a + 5*b*sec(c + d*x))*tan(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_262():
    f = (a + b*sec(c + d*x))*tan(c + d*x)**4
    F = a*x + 3*b*atanh(sin(c + d*x))/(8*d) + (4*a + 3*b*sec(c + d*x))*tan(c + d*x)**3/(12*d) - (8*a + 3*b*sec(c + d*x))*tan(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_263():
    f = (a + b*sec(c + d*x))*tan(c + d*x)**2
    F = -a*x - b*atanh(sin(c + d*x))/(2*d) + (2*a + b*sec(c + d*x))*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_264():
    f = (a + b*sec(c + d*x))*cot(c + d*x)**2
    F = -a*x - (a + b*sec(c + d*x))*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_265():
    f = (a + b*sec(c + d*x))*cot(c + d*x)**4
    F = a*x - (a + b*sec(c + d*x))*cot(c + d*x)**3/(3*d) + (3*a + 2*b*sec(c + d*x))*cot(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_266():
    f = (a + b*sec(c + d*x))*cot(c + d*x)**6
    F = -a*x - (a + b*sec(c + d*x))*cot(c + d*x)**5/(5*d) + (5*a + 4*b*sec(c + d*x))*cot(c + d*x)**3/(15*d) - (15*a + 8*b*sec(c + d*x))*cot(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_267():
    f = (a + b*sec(c + d*x))*cot(c + d*x)**8
    F = a*x - (a + b*sec(c + d*x))*cot(c + d*x)**7/(7*d) + (7*a + 6*b*sec(c + d*x))*cot(c + d*x)**5/(35*d) + (35*a + 16*b*sec(c + d*x))*cot(c + d*x)/(35*d) - (35*a + 24*b*sec(c + d*x))*cot(c + d*x)**3/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_268():
    f = (a + b*sec(c + d*x))**2*tan(c + d*x)**3
    F = a**2*log(cos(c + d*x))/d + 2*a*b*sec(c + d*x)**3/(3*d) - 2*a*b*sec(c + d*x)/d + b**2*sec(c + d*x)**4/(4*d) + (a**2 - b**2)*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_269():
    f = (a + b*sec(c + d*x))**2*tan(c + d*x)
    F = -a**2*log(cos(c + d*x))/d + 2*a*b*sec(c + d*x)/d + b**2*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_270():
    f = (a + b*sec(c + d*x))**2*cot(c + d*x)
    F = a**2*log(cos(c + d*x))/d + (a - b)**2*log(sec(c + d*x) + 1)/(2*d) + (a + b)**2*log(1 - sec(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_271():
    f = (a + b*sec(c + d*x))**2*cot(c + d*x)**3
    F = -a**2*log(cos(c + d*x))/d - a*(a - b)*log(sec(c + d*x) + 1)/(2*d) - a*(a + b)*log(1 - sec(c + d*x))/(2*d) - (a**2 + 2*a*b*sec(c + d*x) + b**2)*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_272():
    f = (a + b*sec(c + d*x))**2*cot(c + d*x)**5
    F = a**2*log(cos(c + d*x))/d + a*(2*a + 3*b*sec(c + d*x))*cot(c + d*x)**2/(4*d) + a*(4*a - 3*b)*log(sec(c + d*x) + 1)/(8*d) + a*(4*a + 3*b)*log(1 - sec(c + d*x))/(8*d) - (a**2 + 2*a*b*sec(c + d*x) + b**2)*cot(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_273():
    f = (a + b*sec(c + d*x))**2*tan(c + d*x)**6
    F = -a**2*x + a**2*tan(c + d*x)**5/(5*d) - a**2*tan(c + d*x)**3/(3*d) + a**2*tan(c + d*x)/d + a*b*tan(c + d*x)**5*sec(c + d*x)/(3*d) - 5*a*b*tan(c + d*x)**3*sec(c + d*x)/(12*d) + 5*a*b*tan(c + d*x)*sec(c + d*x)/(8*d) - 5*a*b*atanh(sin(c + d*x))/(8*d) + b**2*tan(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_274():
    f = (a + b*sec(c + d*x))**2*tan(c + d*x)**4
    F = a**2*x + a**2*tan(c + d*x)**3/(3*d) - a**2*tan(c + d*x)/d + a*b*tan(c + d*x)**3*sec(c + d*x)/(2*d) - 3*a*b*tan(c + d*x)*sec(c + d*x)/(4*d) + 3*a*b*atanh(sin(c + d*x))/(4*d) + b**2*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_275():
    f = (a + b*sec(c + d*x))**2*tan(c + d*x)**2
    F = -a**2*x + a**2*tan(c + d*x)/d + a*b*tan(c + d*x)*sec(c + d*x)/d - a*b*atanh(sin(c + d*x))/d + b**2*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_276():
    f = (a + b*sec(c + d*x))**2*cot(c + d*x)**2
    F = -a**2*x - a**2*cot(c + d*x)/d - 2*a*b*csc(c + d*x)/d - b**2*cot(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_277():
    f = (a + b*sec(c + d*x))**2*cot(c + d*x)**4
    F = a**2*x - a**2*cot(c + d*x)**3/(3*d) + a**2*cot(c + d*x)/d - 2*a*b*csc(c + d*x)**3/(3*d) + 2*a*b*csc(c + d*x)/d - b**2*cot(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_278():
    f = (a + b*sec(c + d*x))**2*cot(c + d*x)**6
    F = -a**2*x - a**2*cot(c + d*x)**5/(5*d) + a**2*cot(c + d*x)**3/(3*d) - a**2*cot(c + d*x)/d - 2*a*b*csc(c + d*x)**5/(5*d) + 4*a*b*csc(c + d*x)**3/(3*d) - 2*a*b*csc(c + d*x)/d - b**2*cot(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_279():
    f = (a + b*sec(c + d*x))**2*cot(c + d*x)**8
    F = a**2*x - a**2*cot(c + d*x)**7/(7*d) + a**2*cot(c + d*x)**5/(5*d) - a**2*cot(c + d*x)**3/(3*d) + a**2*cot(c + d*x)/d - 2*a*b*csc(c + d*x)**7/(7*d) + 6*a*b*csc(c + d*x)**5/(5*d) - 2*a*b*csc(c + d*x)**3/d + 2*a*b*csc(c + d*x)/d - b**2*cot(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_280():
    f = tan(c + d*x)**9/(a + b*sec(c + d*x))
    F = -a*sec(c + d*x)**6/(6*b**2*d) - a*(a**2 - 4*b**2)*sec(c + d*x)**4/(4*b**4*d) - a*(a**4 - 4*a**2*b**2 + 6*b**4)*sec(c + d*x)**2/(2*b**6*d) + sec(c + d*x)**7/(7*b*d) + (a**2 - 4*b**2)*sec(c + d*x)**5/(5*b**3*d) + (a**4 - 4*a**2*b**2 + 6*b**4)*sec(c + d*x)**3/(3*b**5*d) + (a**6 - 4*a**4*b**2 + 6*a**2*b**4 - 4*b**6)*sec(c + d*x)/(b**7*d) - log(cos(c + d*x))/(a*d) - (a**2 - b**2)**4*log(a + b*sec(c + d*x))/(a*b**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_281():
    f = tan(c + d*x)**7/(a + b*sec(c + d*x))
    F = -a*sec(c + d*x)**4/(4*b**2*d) - a*(a**2 - 3*b**2)*sec(c + d*x)**2/(2*b**4*d) + sec(c + d*x)**5/(5*b*d) + (a**2 - 3*b**2)*sec(c + d*x)**3/(3*b**3*d) + (a**4 - 3*a**2*b**2 + 3*b**4)*sec(c + d*x)/(b**5*d) + log(cos(c + d*x))/(a*d) - (a**2 - b**2)**3*log(a + b*sec(c + d*x))/(a*b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_282():
    f = tan(c + d*x)**5/(a + b*sec(c + d*x))
    F = -a*sec(c + d*x)**2/(2*b**2*d) + sec(c + d*x)**3/(3*b*d) + (a**2 - 2*b**2)*sec(c + d*x)/(b**3*d) - log(cos(c + d*x))/(a*d) - (a**2 - b**2)**2*log(a + b*sec(c + d*x))/(a*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_283():
    f = tan(c + d*x)**3/(a + b*sec(c + d*x))
    F = sec(c + d*x)/(b*d) + log(cos(c + d*x))/(a*d) - (a**2 - b**2)*log(a + b*sec(c + d*x))/(a*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_284():
    f = tan(c + d*x)/(a + b*sec(c + d*x))
    F = -log(a + b*sec(c + d*x))/(a*d) - log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_285():
    f = cot(c + d*x)/(a + b*sec(c + d*x))
    F = log(1 - sec(c + d*x))/(d*(2*a + 2*b)) + log(sec(c + d*x) + 1)/(d*(2*a - 2*b)) - b**2*log(a + b*sec(c + d*x))/(a*d*(a**2 - b**2)) + log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_286():
    f = cot(c + d*x)**3/(a + b*sec(c + d*x))
    F = 1/(d*(4*a - 4*b)*(sec(c + d*x) + 1)) - (2*a + 3*b)*log(1 - sec(c + d*x))/(4*d*(a + b)**2) - (2*a - 3*b)*log(sec(c + d*x) + 1)/(4*d*(a - b)**2) + 1/(d*(1 - sec(c + d*x))*(4*a + 4*b)) - b**4*log(a + b*sec(c + d*x))/(a*d*(a**2 - b**2)**2) - log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_287():
    f = cot(c + d*x)**5/(a + b*sec(c + d*x))
    F = -1/(d*(16*a - 16*b)*(sec(c + d*x) + 1)**2) + (8*a**2 + 21*a*b + 15*b**2)*log(1 - sec(c + d*x))/(16*d*(a + b)**3) - (5*a - 7*b)/(16*d*(a - b)**2*(sec(c + d*x) + 1)) + (8*a**2 - 21*a*b + 15*b**2)*log(sec(c + d*x) + 1)/(16*d*(a - b)**3) - (5*a + 7*b)/(16*d*(1 - sec(c + d*x))*(a + b)**2) - 1/(d*(1 - sec(c + d*x))**2*(16*a + 16*b)) - b**6*log(a + b*sec(c + d*x))/(a*d*(a**2 - b**2)**3) + log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_288():
    f = tan(c + d*x)**4/(a + b*sec(c + d*x))
    F = -a*tan(c + d*x)/(b**2*d) + tan(c + d*x)*sec(c + d*x)/(2*b*d) + (2*a**2 - 3*b**2)*atanh(sin(c + d*x))/(2*b**3*d) + x/a - 2*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_289():
    f = tan(c + d*x)**2/(a + b*sec(c + d*x))
    F = atanh(sin(c + d*x))/(b*d) - x/a - 2*sqrt(a - b)*sqrt(a + b)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_290():
    f = tan(c + d*x)**9/(a + b*sec(c + d*x))**2
    F = -2*a*sec(c + d*x)**5/(5*b**3*d) - 4*a*(a**2 - 2*b**2)*sec(c + d*x)**3/(3*b**5*d) - 2*a*(3*a**4 - 8*a**2*b**2 + 6*b**4)*sec(c + d*x)/(b**7*d) + sec(c + d*x)**6/(6*b**2*d) + (3*a**2 - 4*b**2)*sec(c + d*x)**4/(4*b**4*d) + (5*a**4 - 12*a**2*b**2 + 6*b**4)*sec(c + d*x)**2/(2*b**6*d) + (a**2 - b**2)**4/(a*b**8*d*(a + b*sec(c + d*x))) - log(cos(c + d*x))/(a**2*d) + (a**2 - b**2)**3*(7*a**2 + b**2)*log(a + b*sec(c + d*x))/(a**2*b**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_291():
    f = tan(c + d*x)**7/(a + b*sec(c + d*x))**2
    F = -2*a*sec(c + d*x)**3/(3*b**3*d) - 2*a*(2*a**2 - 3*b**2)*sec(c + d*x)/(b**5*d) + sec(c + d*x)**4/(4*b**2*d) + (3*a**2 - 3*b**2)*sec(c + d*x)**2/(2*b**4*d) + (a**2 - b**2)**3/(a*b**6*d*(a + b*sec(c + d*x))) + log(cos(c + d*x))/(a**2*d) + (a**2 - b**2)**2*(5*a**2 + b**2)*log(a + b*sec(c + d*x))/(a**2*b**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_292():
    f = tan(c + d*x)**5/(a + b*sec(c + d*x))**2
    F = -2*a*sec(c + d*x)/(b**3*d) + sec(c + d*x)**2/(2*b**2*d) + (a**2 - b**2)**2/(a*b**4*d*(a + b*sec(c + d*x))) - log(cos(c + d*x))/(a**2*d) + (a**2 - b**2)*(3*a**2 + b**2)*log(a + b*sec(c + d*x))/(a**2*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_293():
    f = tan(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = (a**2 - b**2)/(a*b**2*d*(a + b*sec(c + d*x))) + log(cos(c + d*x))/(a**2*d) + (a**2 + b**2)*log(a + b*sec(c + d*x))/(a**2*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_294():
    f = tan(c + d*x)/(a + b*sec(c + d*x))**2
    F = 1/(a*d*(a + b*sec(c + d*x))) - log(a + b*sec(c + d*x))/(a**2*d) - log(cos(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_295():
    f = cot(c + d*x)/(a + b*sec(c + d*x))**2
    F = log(1 - sec(c + d*x))/(2*d*(a + b)**2) + log(sec(c + d*x) + 1)/(2*d*(a - b)**2) + b**2/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - b**2*(3*a**2 - b**2)*log(a + b*sec(c + d*x))/(a**2*d*(a**2 - b**2)**2) + log(cos(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_296():
    f = cot(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = -(a - 2*b)*log(sec(c + d*x) + 1)/(2*d*(a - b)**3) - (a + 2*b)*log(1 - sec(c + d*x))/(2*d*(a + b)**3) + 1/(4*d*(a - b)**2*(sec(c + d*x) + 1)) + 1/(4*d*(1 - sec(c + d*x))*(a + b)**2) + b**4/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - b**4*(5*a**2 - b**2)*log(a + b*sec(c + d*x))/(a**2*d*(a**2 - b**2)**3) - log(cos(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_297():
    f = cot(c + d*x)**5/(a + b*sec(c + d*x))**2
    F = (4*a**2 + 13*a*b + 12*b**2)*log(1 - sec(c + d*x))/(8*d*(a + b)**4) - 1/(16*d*(a - b)**2*(sec(c + d*x) + 1)**2) - (5*a - 9*b)/(16*d*(a - b)**3*(sec(c + d*x) + 1)) + (4*a**2 - 13*a*b + 12*b**2)*log(sec(c + d*x) + 1)/(8*d*(a - b)**4) - (5*a + 9*b)/(16*d*(1 - sec(c + d*x))*(a + b)**3) - 1/(16*d*(1 - sec(c + d*x))**2*(a + b)**2) + b**6/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - b**6*(7*a**2 - b**2)*log(a + b*sec(c + d*x))/(a**2*d*(a**2 - b**2)**4) + log(cos(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_298():
    f = tan(c + d*x)**4/(a + b*sec(c + d*x))**2
    F = -2*a*atanh(sin(c + d*x))/(b**3*d) + tan(c + d*x)/(b*d*(a*cos(c + d*x) + b)) + (2*a**2 - b**2)*sin(c + d*x)/(a*b**2*d*(a*cos(c + d*x) + b)) + x/a**2 + 2*sqrt(a - b)*sqrt(a + b)*(2*a**2 + b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_299():
    f = tan(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = tan(c + d*x)/(a*d*(a + b*sec(c + d*x))) + 2*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*sqrt(a - b)*sqrt(a + b)) - x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_300():
    f = cot(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = sin(c + d*x)/(2*d*(a - b)**2*(cos(c + d*x) + 1)) - sin(c + d*x)/(2*d*(1 - cos(c + d*x))*(a + b)**2) + b**4*sin(c + d*x)/(a*d*(a**2 - b**2)**2*(a*cos(c + d*x) + b)) - 2*b**5*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - 4*b**3*(2*a**2 - b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_301():
    f = cot(c + d*x)**4/(a + b*sec(c + d*x))**2
    F = sin(c + d*x)/(12*d*(a - b)**2*(cos(c + d*x) + 1)) + sin(c + d*x)/(12*d*(a - b)**2*(cos(c + d*x) + 1)**2) - (3*a - 5*b)*sin(c + d*x)/(4*d*(a - b)**3*(cos(c + d*x) + 1)) - sin(c + d*x)/(12*d*(1 - cos(c + d*x))*(a + b)**2) + (3*a + 5*b)*sin(c + d*x)/(4*d*(1 - cos(c + d*x))*(a + b)**3) - sin(c + d*x)/(12*d*(1 - cos(c + d*x))**2*(a + b)**2) + b**6*sin(c + d*x)/(a*d*(a**2 - b**2)**3*(a*cos(c + d*x) + b)) - 2*b**7*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - 4*b**5*(3*a**2 - b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_302():
    f = (e*tan(c + d*x))**(sympy.S(5)/2)/(a + b*sec(c + d*x))
    F = ((Symbol('a') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('e'))**(Integer(2)) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e((Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(4))**(Integer(-1)))) + (Symbol('d') * x)), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sin(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('e') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * ((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_303():
    f = (e*tan(c + d*x))**(sympy.S(3)/2)/(a + b*sec(c + d*x))
    F = ((Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Integer(1) + sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * (Symbol('e'))**(Integer(2)) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Integer(1) + sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + (((Symbol('e'))**(Integer(2)) * sympy.elliptic_f((Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(4))**(Integer(-1)))) + (Symbol('d') * x)), Integer(2)) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt(sympy.sin(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_304():
    f = sqrt(e*tan(c + d*x))/(a + b*sec(c + d*x))
    F = (Integer(-1) * ((sympy.sqrt(Symbol('e')) * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt(Symbol('e')) * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt(Symbol('e')) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Symbol('e')) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_305():
    f = 1/(sqrt(e*tan(c + d*x))*(a + b*sec(c + d*x)))
    F = (Integer(-1) * (sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * (sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1)))) + (sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * Symbol('d') * sympy.sqrt(Symbol('e'))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('b') * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Integer(1) + sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('b') * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Integer(1) + sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))) * Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_306():
    f = 1/((e*tan(c + d*x))**(sympy.S(3)/2)*(a + b*sec(c + d*x)))
    F = ((Symbol('a') * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((Integer(-1) * (sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')((sympy.sqrt((Symbol('a') + (Integer(-1) * Symbol('b')))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1))), sympy.asin((sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (sympy.sqrt((Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e((Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(4))**(Integer(-1)))) + (Symbol('d') * x)), Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt(sympy.sin(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * ((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(3))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_307():
    f = 1/((e*tan(c + d*x))**(sympy.S(5)/2)*(a + b*sec(c + d*x)))
    F = ((Symbol('a') * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.atan((Integer(1) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1))))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.atan((Integer(1) + ((sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(Symbol('e')))**(Integer(-1)))))) * ((sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (Integer(-1) * (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.log((sympy.sqrt(Symbol('e')) + (sympy.sqrt(Symbol('e')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) + (sympy.sqrt(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x))))))))) * ((Integer(2) * sympy.sqrt(Integer(2)) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * Symbol('e') * ((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + (Integer(-1) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2))))))))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Integer(1) + sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Integer(1) + sympy.sin((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-1)) * sympy.sqrt(sympy.sin((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))) + ((Symbol('b') * sympy.elliptic_f((Symbol('c') + (Integer(-1) * (sympy.pi * (Integer(4))**(Integer(-1)))) + (Symbol('d') * x)), Integer(2)) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt(sympy.sin(((Integer(2) * Symbol('c')) + (Integer(2) * Symbol('d') * x))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('e'))**(Integer(2)) * sympy.sqrt((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_308():
    f = sqrt(a + b*sec(c + d*x))*tan(c + d*x)**5
    F = -2*sqrt(a)*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/d - 6*a*(a + b*sec(c + d*x))**(sympy.S(7)/2)/(7*b**4*d) - 2*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - 2*b**2)/(3*b**4*d) + 2*sqrt(a + b*sec(c + d*x))/d + 2*(a + b*sec(c + d*x))**(sympy.S(9)/2)/(9*b**4*d) + (a + b*sec(c + d*x))**(sympy.S(5)/2)*(6*a**2 - 4*b**2)/(5*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_309():
    f = sqrt(a + b*sec(c + d*x))*tan(c + d*x)**3
    F = 2*sqrt(a)*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/d - 2*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)/(3*b**2*d) - 2*sqrt(a + b*sec(c + d*x))/d + 2*(a + b*sec(c + d*x))**(sympy.S(5)/2)/(5*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_310():
    f = sqrt(a + b*sec(c + d*x))*tan(c + d*x)
    F = -2*sqrt(a)*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/d + 2*sqrt(a + b*sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_311():
    f = sqrt(a + b*sec(c + d*x))*cot(c + d*x)
    F = 2*sqrt(a)*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/d - sqrt(a - b)*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a - b))/d - sqrt(a + b)*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a + b))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_312():
    f = sqrt(a + b*sec(c + d*x))*cot(c + d*x)**3
    F = -2*sqrt(a)*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/d + a*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) + a*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a - b))/(d*sqrt(a - b)) + 3*b*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a + b))/(4*d*sqrt(a + b)) - 3*b*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a - b))/(4*d*sqrt(a - b)) - sqrt(a + b*sec(c + d*x))*cot(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_313():
    f = sqrt(a + b*sec(c + d*x))*tan(c + d*x)**2
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_314():
    f = sqrt(a + b*sec(c + d*x))
    F = Integer(-1) * ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Symbol('a') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_315():
    f = sqrt(a + b*sec(c + d*x))*cot(c + d*x)**2
    F = ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Symbol('a') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_316():
    f = tan(c + d*x)**5/sqrt(a + b*sec(c + d*x))
    F = -6*a*(a + b*sec(c + d*x))**(sympy.S(5)/2)/(5*b**4*d) - 2*a*sqrt(a + b*sec(c + d*x))*(a**2 - 2*b**2)/(b**4*d) + 2*(a + b*sec(c + d*x))**(sympy.S(7)/2)/(7*b**4*d) + (a + b*sec(c + d*x))**(sympy.S(3)/2)*(6*a**2 - 4*b**2)/(3*b**4*d) - 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_317():
    f = tan(c + d*x)**3/sqrt(a + b*sec(c + d*x))
    F = -2*a*sqrt(a + b*sec(c + d*x))/(b**2*d) + 2*(a + b*sec(c + d*x))**(sympy.S(3)/2)/(3*b**2*d) + 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_318():
    f = tan(c + d*x)/sqrt(a + b*sec(c + d*x))
    F = -2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_319():
    f = cot(c + d*x)/sqrt(a + b*sec(c + d*x))
    F = -atanh(sqrt(a + b*sec(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) - atanh(sqrt(a + b*sec(c + d*x))/sqrt(a - b))/(d*sqrt(a - b)) + 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_320():
    f = cot(c + d*x)**3/sqrt(a + b*sec(c + d*x))
    F = b*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a + b))/(4*d*(a + b)**(sympy.S(3)/2)) - b*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a - b))/(4*d*(a - b)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))/(d*(4*a - 4*b)*(sec(c + d*x) + 1)) + atanh(sqrt(a + b*sec(c + d*x))/sqrt(a + b))/(d*sqrt(a + b)) + atanh(sqrt(a + b*sec(c + d*x))/sqrt(a - b))/(d*sqrt(a - b)) + sqrt(a + b*sec(c + d*x))/(d*(1 - sec(c + d*x))*(4*a + 4*b)) - 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_321():
    f = tan(c + d*x)**2/sqrt(a + b*sec(c + d*x))
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_322():
    f = 1/sqrt(a + b*sec(c + d*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_323():
    f = cot(c + d*x)**2/sqrt(a + b*sec(c + d*x))
    F = ((sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.cot((Symbol('c') + (Symbol('d') * x))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_324():
    f = tan(c + d*x)**5/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)/(b**4*d) + 2*(a + b*sec(c + d*x))**(sympy.S(5)/2)/(5*b**4*d) + sqrt(a + b*sec(c + d*x))*(6*a**2 - 4*b**2)/(b**4*d) + 2*(a**2 - b**2)**2/(a*b**4*d*sqrt(a + b*sec(c + d*x))) - 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_325():
    f = tan(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*sqrt(a + b*sec(c + d*x))/(b**2*d) + (2*a**2 - 2*b**2)/(a*b**2*d*sqrt(a + b*sec(c + d*x))) + 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_326():
    f = tan(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2/(a*d*sqrt(a + b*sec(c + d*x))) - 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_327():
    f = cot(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -atanh(sqrt(a + b*sec(c + d*x))/sqrt(a + b))/(d*(a + b)**(sympy.S(3)/2)) - atanh(sqrt(a + b*sec(c + d*x))/sqrt(a - b))/(d*(a - b)**(sympy.S(3)/2)) + 2*b**2/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + 2*atanh(sqrt(a + b*sec(c + d*x))/sqrt(a))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_328():
    f = tan(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_329():
    f = (a + b*sec(c + d*x))**(sympy.S(-3)/2)
    F = ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_330():
    f = (d*tan(e + f*x))**n*(a + b*sec(e + f*x))**3
    F = a**3*(d*tan(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(n + 1)) + 3*a**2*b*(d*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + 1)*hyper((n/2 + sympy.S.Half, n/2 + 1), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)*sec(e + f*x)/(d*f*(n + 1)) + 3*a*b**2*(d*tan(e + f*x))**(n + 1)/(d*f*(n + 1)) + b**3*(d*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + 2)*hyper((n/2 + sympy.S.Half, n/2 + 2), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)*sec(e + f*x)**3/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_331():
    f = (d*tan(e + f*x))**n*(a + b*sec(e + f*x))**2
    F = a**2*(d*tan(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(n + 1)) + 2*a*b*(d*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + 1)*hyper((n/2 + sympy.S.Half, n/2 + 1), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)*sec(e + f*x)/(d*f*(n + 1)) + b**2*(d*tan(e + f*x))**(n + 1)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_332():
    f = (d*tan(e + f*x))**n*(a + b*sec(e + f*x))
    F = a*(d*tan(e + f*x))**(n + 1)*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(d*f*(n + 1)) + b*(d*tan(e + f*x))**(n + 1)*(cos(e + f*x)**2)**(n/2 + 1)*hyper((n/2 + sympy.S.Half, n/2 + 1), (n/2 + sympy.S(3)/2,), sin(e + f*x)**2)*sec(e + f*x)/(d*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_333():
    f = (d*tan(e + f*x))**n/(a + b*sec(e + f*x))
    F = d*(d*tan(e + f*x))**(n - 1)*(b*(sec(e + f*x) + 1)/(a + b*sec(e + f*x)))**(sympy.S.Half - n/2)*(-b*(1 - sec(e + f*x))/(a + b*sec(e + f*x)))**(sympy.S.Half - n/2)*appellf1(1 - n, sympy.S.Half - n/2, sympy.S.Half - n/2, 2 - n, (a - b)/(a + b*sec(e + f*x)), (a + b)/(a + b*sec(e + f*x)))/(a*f*(1 - n)) + d*(d*tan(e + f*x))**(n - 1)*tan(e + f*x)**2*hyper((1, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), -tan(e + f*x)**2)/(a*f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_334():
    f = (e*tan(c + d*x))**m*(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_335():
    f = (e*tan(c + d*x))**m*sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_336():
    f = (e*tan(c + d*x))**m/sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m')) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_337():
    f = (e*tan(c + d*x))**m/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m')) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_338():
    f = (e*tan(c + d*x))**m*(a + b*sec(c + d*x))**n
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * ((Symbol('e') * sympy.tan((Symbol('c') + (Symbol('d') * x)))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_339():
    f = (a + b*sec(c + d*x))**n*tan(c + d*x)**5
    F = -a*(a + b*sec(c + d*x))**(n + 1)*(a**2 - 2*b**2)/(b**4*d*(n + 1)) - 3*a*(a + b*sec(c + d*x))**(n + 3)/(b**4*d*(n + 3)) + (a + b*sec(c + d*x))**(n + 2)*(3*a**2 - 2*b**2)/(b**4*d*(n + 2)) + (a + b*sec(c + d*x))**(n + 4)/(b**4*d*(n + 4)) - (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_340():
    f = (a + b*sec(c + d*x))**n*tan(c + d*x)**3
    F = -a*(a + b*sec(c + d*x))**(n + 1)/(b**2*d*(n + 1)) + (a + b*sec(c + d*x))**(n + 2)/(b**2*d*(n + 2)) + (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_341():
    f = (a + b*sec(c + d*x))**n*tan(c + d*x)
    F = -(a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_342():
    f = (a + b*sec(c + d*x))**n*cot(c + d*x)
    F = -(a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a + b))/(d*(2*a + 2*b)*(n + 1)) - (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a - b))/(d*(2*a - 2*b)*(n + 1)) + (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_343():
    f = (a + b*sec(c + d*x))**n*cot(c + d*x)**3
    F = b*(a + b*sec(c + d*x))**(n + 1)*hyper((2, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a + b))/(4*d*(a + b)**2*(n + 1)) - b*(a + b*sec(c + d*x))**(n + 1)*hyper((2, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a - b))/(4*d*(a - b)**2*(n + 1)) + (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a + b))/(d*(2*a + 2*b)*(n + 1)) + (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), (a + b*sec(c + d*x))/(a - b))/(d*(2*a - 2*b)*(n + 1)) - (a + b*sec(c + d*x))**(n + 1)*hyper((1, n + 1), (n + 2,), 1 + b*sec(c + d*x)/a)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_344():
    f = (a + b*sec(c + d*x))**n*tan(c + d*x)**4
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**(Integer(4))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_345():
    f = (a + b*sec(c + d*x))**n*tan(c + d*x)**2
    F = ((sympy.sqrt(Integer(2)) * (Symbol('a') + Symbol('b')) * sympy.appellf1((Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1)), (Integer(-1) + (Integer(-1) * Symbol('n'))), (Integer(3) * (Integer(2))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))), ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))**(Symbol('n')) * (Symbol('b') * Symbol('d') * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(2)) * Symbol('a') * sympy.appellf1((Integer(2))**(Integer(-1)), (Integer(2))**(Integer(-1)), (Integer(-1) * Symbol('n')), (Integer(3) * (Integer(2))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))), ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))**(Symbol('n')) * (Symbol('b') * Symbol('d') * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')), x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_346():
    f = (a + b*sec(c + d*x))**n*cot(c + d*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_347():
    f = (a + b*sec(c + d*x))**n*cot(c + d*x)**4
    F = sympy.Function('Unintegrable')(((sympy.cot((Symbol('c') + (Symbol('d') * x))))**(Integer(4)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_348():
    f = (a + b*sec(c + d*x))**n*tan(c + d*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * (sympy.tan((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_349():
    f = (a + b*sec(c + d*x))**n*sqrt(tan(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * sympy.sqrt(sympy.tan((Symbol('c') + (Symbol('d') * x))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_350():
    f = (a + b*sec(c + d*x))**n/sqrt(tan(c + d*x))
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * (sympy.sqrt(sympy.tan((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_4_d_tan_pow_n_a_plus_b_sec_pow_m_351():
    f = (a + b*sec(c + d*x))**n/tan(c + d*x)**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Symbol('n')) * ((sympy.tan((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F

