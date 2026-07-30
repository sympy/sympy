"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.1.2 (d sec)^n (a+b sec)^m.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_1():
    f = (a*sec(c + d*x) + a)*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + a*tan(c + d*x)/d + 3*a*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_2():
    f = (a*sec(c + d*x) + a)*sec(c + d*x)**3
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*tan(c + d*x)/d + a*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_3():
    f = (a*sec(c + d*x) + a)*sec(c + d*x)**2
    F = a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*tan(c + d*x)/d + a*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_4():
    f = (a*sec(c + d*x) + a)*sec(c + d*x)
    F = a*tan(c + d*x)/d + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_5():
    f = a*sec(c + d*x) + a
    F = a*x + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_6():
    f = (a*sec(c + d*x) + a)*cos(c + d*x)
    F = a*x + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_7():
    f = (a*sec(c + d*x) + a)*cos(c + d*x)**2
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_8():
    f = (a*sec(c + d*x) + a)*cos(c + d*x)**3
    F = a*x/2 - a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)*cos(c + d*x)/(2*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_9():
    f = (a*sec(c + d*x) + a)*cos(c + d*x)**4
    F = 3*a*x/8 - a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_10():
    f = (a*sec(c + d*x) + a)**2*sec(c + d*x)**4
    F = 3*a**2*tan(c + d*x)**3/(5*d) + a**2*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a**2*tan(c + d*x)*sec(c + d*x)**3/(2*d) + 3*a**2*tan(c + d*x)*sec(c + d*x)/(4*d) + 9*a**2*tan(c + d*x)/(5*d) + 3*a**2*atanh(sin(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_11():
    f = (a*sec(c + d*x) + a)**2*sec(c + d*x)**3
    F = 2*a**2*tan(c + d*x)**3/(3*d) + a**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 7*a**2*tan(c + d*x)*sec(c + d*x)/(8*d) + 2*a**2*tan(c + d*x)/d + 7*a**2*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_12():
    f = (a*sec(c + d*x) + a)**2*sec(c + d*x)**2
    F = a**2*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a**2*tan(c + d*x)*sec(c + d*x)/d + 5*a**2*tan(c + d*x)/(3*d) + a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_13():
    f = (a*sec(c + d*x) + a)**2*sec(c + d*x)
    F = a**2*tan(c + d*x)*sec(c + d*x)/(2*d) + 2*a**2*tan(c + d*x)/d + 3*a**2*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_14():
    f = (a*sec(c + d*x) + a)**2
    F = a**2*x + a**2*tan(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_15():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)
    F = 2*a**2*x + a**2*sin(c + d*x)/d + a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_16():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**2
    F = 3*a**2*x/2 + a**2*sin(c + d*x)*cos(c + d*x)/(2*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_17():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**3
    F = a**2*x - a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)*cos(c + d*x)/d + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_18():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**4
    F = 7*a**2*x/8 - 2*a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 7*a**2*sin(c + d*x)*cos(c + d*x)/(8*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_19():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**5
    F = 3*a**2*x/4 + a**2*sin(c + d*x)**5/(5*d) - a**2*sin(c + d*x)**3/d + a**2*sin(c + d*x)*cos(c + d*x)**3/(2*d) + 3*a**2*sin(c + d*x)*cos(c + d*x)/(4*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_20():
    f = (a*sec(c + d*x) + a)**3*sec(c + d*x)**3
    F = a**3*tan(c + d*x)**5/(5*d) + 5*a**3*tan(c + d*x)**3/(3*d) + 3*a**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 13*a**3*tan(c + d*x)*sec(c + d*x)/(8*d) + 4*a**3*tan(c + d*x)/d + 13*a**3*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_21():
    f = (a*sec(c + d*x) + a)**3*sec(c + d*x)**2
    F = a**3*tan(c + d*x)**3/d + a**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 15*a**3*tan(c + d*x)*sec(c + d*x)/(8*d) + 4*a**3*tan(c + d*x)/d + 15*a**3*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_22():
    f = (a*sec(c + d*x) + a)**3*sec(c + d*x)
    F = a**3*tan(c + d*x)**3/(3*d) + 3*a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 4*a**3*tan(c + d*x)/d + 5*a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_23():
    f = (a*sec(c + d*x) + a)**3
    F = a**3*x + 5*a**3*tan(c + d*x)/(2*d) + 7*a**3*atanh(sin(c + d*x))/(2*d) + (a**3*sec(c + d*x) + a**3)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_24():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)
    F = 3*a**3*x + a**3*sin(c + d*x)/d + a**3*tan(c + d*x)/d + 3*a**3*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_25():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**2
    F = 7*a**3*x/2 + a**3*sin(c + d*x)*cos(c + d*x)/(2*d) + 3*a**3*sin(c + d*x)/d + a**3*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_26():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**3
    F = 5*a**3*x/2 - a**3*sin(c + d*x)**3/(3*d) + 3*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) + 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_27():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**4
    F = 15*a**3*x/8 - a**3*sin(c + d*x)**3/d + a**3*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 15*a**3*sin(c + d*x)*cos(c + d*x)/(8*d) + 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_28():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**5
    F = 13*a**3*x/8 + a**3*sin(c + d*x)**5/(5*d) - 5*a**3*sin(c + d*x)**3/(3*d) + 3*a**3*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 13*a**3*sin(c + d*x)*cos(c + d*x)/(8*d) + 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_29():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**6
    F = 23*a**3*x/16 + 3*a**3*sin(c + d*x)**5/(5*d) - 7*a**3*sin(c + d*x)**3/(3*d) + a**3*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 23*a**3*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 23*a**3*sin(c + d*x)*cos(c + d*x)/(16*d) + 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_30():
    f = (a*sec(c + d*x) + a)**4*sec(c + d*x)**3
    F = 4*a**4*tan(c + d*x)**5/(5*d) + 4*a**4*tan(c + d*x)**3/d + a**4*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 41*a**4*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 49*a**4*tan(c + d*x)*sec(c + d*x)/(16*d) + 8*a**4*tan(c + d*x)/d + 49*a**4*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_31():
    f = (a*sec(c + d*x) + a)**4*sec(c + d*x)**2
    F = a**4*tan(c + d*x)**5/(5*d) + 8*a**4*tan(c + d*x)**3/(3*d) + a**4*tan(c + d*x)*sec(c + d*x)**3/d + 7*a**4*tan(c + d*x)*sec(c + d*x)/(2*d) + 8*a**4*tan(c + d*x)/d + 7*a**4*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_32():
    f = (a*sec(c + d*x) + a)**4*sec(c + d*x)
    F = 4*a**4*tan(c + d*x)**3/(3*d) + a**4*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 27*a**4*tan(c + d*x)*sec(c + d*x)/(8*d) + 8*a**4*tan(c + d*x)/d + 35*a**4*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_33():
    f = (a*sec(c + d*x) + a)**4
    F = a**4*x + 5*a**4*tan(c + d*x)/d + 6*a**4*atanh(sin(c + d*x))/d + (a**2*sec(c + d*x) + a**2)**2*tan(c + d*x)/(3*d) + (4*a**4*sec(c + d*x) + 4*a**4)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_34():
    f = (a*sec(c + d*x) + a)**4*cos(c + d*x)
    F = 4*a**4*x + a**4*sin(c + d*x)/d + a**4*tan(c + d*x)*sec(c + d*x)/(2*d) + 4*a**4*tan(c + d*x)/d + 13*a**4*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_35():
    f = (a*sec(c + d*x) + a)**4*cos(c + d*x)**2
    F = 13*a**4*x/2 + a**4*sin(c + d*x)*cos(c + d*x)/(2*d) + 4*a**4*sin(c + d*x)/d + a**4*tan(c + d*x)/d + 4*a**4*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_36():
    f = (a*sec(c + d*x) + a)**4*cos(c + d*x)**3
    F = 6*a**4*x - a**4*sin(c + d*x)**3/(3*d) + 2*a**4*sin(c + d*x)*cos(c + d*x)/d + 7*a**4*sin(c + d*x)/d + a**4*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_37():
    f = (a*sec(c + d*x) + a)**4*cos(c + d*x)**4
    F = 35*a**4*x/8 - 4*a**4*sin(c + d*x)**3/(3*d) + a**4*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 27*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) + 8*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_38():
    f = (a*sec(c + d*x) + a)**4*cos(c + d*x)**5
    F = 7*a**4*x/2 + a**4*sin(c + d*x)**5/(5*d) - 8*a**4*sin(c + d*x)**3/(3*d) + a**4*sin(c + d*x)*cos(c + d*x)**3/d + 7*a**4*sin(c + d*x)*cos(c + d*x)/(2*d) + 8*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_39():
    f = (a*sec(c + d*x) + a)**4*cos(c + d*x)**6
    F = 49*a**4*x/16 + 4*a**4*sin(c + d*x)**5/(5*d) - 4*a**4*sin(c + d*x)**3/d + a**4*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 41*a**4*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 49*a**4*sin(c + d*x)*cos(c + d*x)/(16*d) + 8*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_40():
    f = (a*sec(c + d*x) + a)**4*cos(c + d*x)**7
    F = 11*a**4*x/4 - a**4*sin(c + d*x)**7/(7*d) + 9*a**4*sin(c + d*x)**5/(5*d) - 16*a**4*sin(c + d*x)**3/(3*d) + 2*a**4*sin(c + d*x)*cos(c + d*x)**5/(3*d) + 11*a**4*sin(c + d*x)*cos(c + d*x)**3/(6*d) + 11*a**4*sin(c + d*x)*cos(c + d*x)/(4*d) + 8*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_41():
    f = (a*sec(c + d*x) + a)**5*sec(c + d*x)**3
    F = a**5*tan(c + d*x)**7/(7*d) + 13*a**5*tan(c + d*x)**5/(5*d) + 28*a**5*tan(c + d*x)**3/(3*d) + 5*a**5*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 85*a**5*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 93*a**5*tan(c + d*x)*sec(c + d*x)/(16*d) + 16*a**5*tan(c + d*x)/d + 93*a**5*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_42():
    f = sec(c + d*x)**5/(a*sec(c + d*x) + a)
    F = -tan(c + d*x)*sec(c + d*x)**3/(d*(a*sec(c + d*x) + a)) + 4*tan(c + d*x)**3/(3*a*d) - 3*tan(c + d*x)*sec(c + d*x)/(2*a*d) + 4*tan(c + d*x)/(a*d) - 3*atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_43():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)
    F = -tan(c + d*x)*sec(c + d*x)**2/(d*(a*sec(c + d*x) + a)) + 3*tan(c + d*x)*sec(c + d*x)/(2*a*d) - 2*tan(c + d*x)/(a*d) + 3*atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_44():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)
    F = tan(c + d*x)/(d*(a*sec(c + d*x) + a)) + tan(c + d*x)/(a*d) - atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_45():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)
    F = -tan(c + d*x)/(d*(a*sec(c + d*x) + a)) + atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_46():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)
    F = tan(c + d*x)/(d*(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_47():
    f = 1/(a*sec(c + d*x) + a)
    F = -tan(c + d*x)/(d*(a*sec(c + d*x) + a)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_48():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)) - x/a + 2*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_49():
    f = cos(c + d*x)**2/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)/(d*(a*sec(c + d*x) + a)) + 3*x/(2*a) + 3*sin(c + d*x)*cos(c + d*x)/(2*a*d) - 2*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_50():
    f = cos(c + d*x)**3/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**2/(d*(a*sec(c + d*x) + a)) - 3*x/(2*a) - 4*sin(c + d*x)**3/(3*a*d) - 3*sin(c + d*x)*cos(c + d*x)/(2*a*d) + 4*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_51():
    f = cos(c + d*x)**4/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**3/(d*(a*sec(c + d*x) + a)) + 15*x/(8*a) + 4*sin(c + d*x)**3/(3*a*d) + 5*sin(c + d*x)*cos(c + d*x)**3/(4*a*d) + 15*sin(c + d*x)*cos(c + d*x)/(8*a*d) - 4*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_52():
    f = sec(c + d*x)**5/(a*sec(c + d*x) + a)**2
    F = -tan(c + d*x)*sec(c + d*x)**3/(3*d*(a*sec(c + d*x) + a)**2) + 7*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) - 16*tan(c + d*x)/(3*a**2*d) + 7*atanh(sin(c + d*x))/(2*a**2*d) - 8*tan(c + d*x)*sec(c + d*x)**2/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_53():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**2
    F = -tan(c + d*x)*sec(c + d*x)**2/(3*d*(a*sec(c + d*x) + a)**2) + 4*tan(c + d*x)/(3*a**2*d) - 2*atanh(sin(c + d*x))/(a**2*d) + 2*tan(c + d*x)/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_54():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = tan(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) + atanh(sin(c + d*x))/(a**2*d) - 5*tan(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_55():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**2
    F = 2*tan(c + d*x)/(3*d*(a**2*sec(c + d*x) + a**2)) - tan(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_56():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**2
    F = tan(c + d*x)/(3*d*(a**2*sec(c + d*x) + a**2)) + tan(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_57():
    f = (a*sec(c + d*x) + a)**(-2)
    F = -tan(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) + x/a**2 - 4*tan(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_58():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) - 2*x/a**2 + 10*sin(c + d*x)/(3*a**2*d) - 2*sin(c + d*x)/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_59():
    f = cos(c + d*x)**2/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2) + 7*x/(2*a**2) + 7*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) - 16*sin(c + d*x)/(3*a**2*d) - 8*sin(c + d*x)*cos(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_60():
    f = cos(c + d*x)**3/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**2/(3*d*(a*sec(c + d*x) + a)**2) - 5*x/a**2 - 4*sin(c + d*x)**3/(a**2*d) - 5*sin(c + d*x)*cos(c + d*x)/(a**2*d) + 12*sin(c + d*x)/(a**2*d) - 10*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_61():
    f = sec(c + d*x)**6/(a*sec(c + d*x) + a)**3
    F = -76*tan(c + d*x)*sec(c + d*x)**2/(15*d*(a**3*sec(c + d*x) + a**3)) - tan(c + d*x)*sec(c + d*x)**4/(5*d*(a*sec(c + d*x) + a)**3) - 11*tan(c + d*x)*sec(c + d*x)**3/(15*a*d*(a*sec(c + d*x) + a)**2) + 13*tan(c + d*x)*sec(c + d*x)/(2*a**3*d) - 152*tan(c + d*x)/(15*a**3*d) + 13*atanh(sin(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_62():
    f = sec(c + d*x)**5/(a*sec(c + d*x) + a)**3
    F = 3*tan(c + d*x)/(d*(a**3*sec(c + d*x) + a**3)) - tan(c + d*x)*sec(c + d*x)**3/(5*d*(a*sec(c + d*x) + a)**3) - 3*tan(c + d*x)*sec(c + d*x)**2/(5*a*d*(a*sec(c + d*x) + a)**2) + 9*tan(c + d*x)/(5*a**3*d) - 3*atanh(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_63():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**3
    F = -29*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) - tan(c + d*x)*sec(c + d*x)**2/(5*d*(a*sec(c + d*x) + a)**3) + 7*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2) + atanh(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_64():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**3
    F = 7*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) + tan(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) - 8*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_65():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**3
    F = tan(c + d*x)/(5*d*(a**3*sec(c + d*x) + a**3)) - tan(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) + tan(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_66():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**3
    F = 2*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) + tan(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) + 2*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_67():
    f = (a*sec(c + d*x) + a)**(-3)
    F = -22*tan(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) - tan(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) - 7*tan(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_68():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**3
    F = -3*sin(c + d*x)/(d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) - 3*sin(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**2) - 3*x/a**3 + 24*sin(c + d*x)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_69():
    f = cos(c + d*x)**2/(a*sec(c + d*x) + a)**3
    F = -76*sin(c + d*x)*cos(c + d*x)/(15*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3) - 11*sin(c + d*x)*cos(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2) + 13*x/(2*a**3) + 13*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) - 152*sin(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_70():
    f = sec(c + d*x)**7/(a*sec(c + d*x) + a)**4
    F = -tan(c + d*x)*sec(c + d*x)**5/(7*d*(a*sec(c + d*x) + a)**4) - 2*tan(c + d*x)*sec(c + d*x)**4/(5*a*d*(a*sec(c + d*x) + a)**3) + 21*tan(c + d*x)*sec(c + d*x)/(2*a**4*d) - 576*tan(c + d*x)/(35*a**4*d) + 21*atanh(sin(c + d*x))/(2*a**4*d) - 288*tan(c + d*x)*sec(c + d*x)**2/(35*a**4*d*(sec(c + d*x) + 1)) - 43*tan(c + d*x)*sec(c + d*x)**3/(35*a**4*d*(sec(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_71():
    f = sec(c + d*x)**6/(a*sec(c + d*x) + a)**4
    F = -tan(c + d*x)*sec(c + d*x)**4/(7*d*(a*sec(c + d*x) + a)**4) - 12*tan(c + d*x)*sec(c + d*x)**3/(35*a*d*(a*sec(c + d*x) + a)**3) + 244*tan(c + d*x)/(105*a**4*d) - 4*atanh(sin(c + d*x))/(a**4*d) + 4*tan(c + d*x)/(a**4*d*(sec(c + d*x) + 1)) - 88*tan(c + d*x)*sec(c + d*x)**2/(105*a**4*d*(sec(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_72():
    f = sec(c + d*x)**5/(a*sec(c + d*x) + a)**4
    F = -tan(c + d*x)*sec(c + d*x)**3/(7*d*(a*sec(c + d*x) + a)**4) - 2*tan(c + d*x)*sec(c + d*x)**2/(7*a*d*(a*sec(c + d*x) + a)**3) + atanh(sin(c + d*x))/(a**4*d) - 43*tan(c + d*x)/(21*a**4*d*(sec(c + d*x) + 1)) + 11*tan(c + d*x)/(21*a**4*d*(sec(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_73():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**4
    F = tan(c + d*x)/(5*d*(a**4*sec(c + d*x) + a**4)) - 8*tan(c + d*x)/(35*d*(a**2*sec(c + d*x) + a**2)**2) + tan(c + d*x)*sec(c + d*x)**3/(7*d*(a*sec(c + d*x) + a)**4) + 3*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_74():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**4
    F = 13*tan(c + d*x)/(105*d*(a**4*sec(c + d*x) + a**4)) + 13*tan(c + d*x)/(105*d*(a**2*sec(c + d*x) + a**2)**2) + tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) - 11*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_75():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**4
    F = 8*tan(c + d*x)/(105*d*(a**4*sec(c + d*x) + a**4)) + 8*tan(c + d*x)/(105*d*(a**2*sec(c + d*x) + a**2)**2) - tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) + 4*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_76():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**4
    F = 2*tan(c + d*x)/(35*d*(a**4*sec(c + d*x) + a**4)) + 2*tan(c + d*x)/(35*d*(a**2*sec(c + d*x) + a**2)**2) + tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) + 3*tan(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_77():
    f = (a*sec(c + d*x) + a)**(-4)
    F = -tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) - 2*tan(c + d*x)/(7*a*d*(a*sec(c + d*x) + a)**3) + x/a**4 - 32*tan(c + d*x)/(21*a**4*d*(sec(c + d*x) + 1)) - 11*tan(c + d*x)/(21*a**4*d*(sec(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_78():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**4
    F = -sin(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) - 12*sin(c + d*x)/(35*a*d*(a*sec(c + d*x) + a)**3) - 4*x/a**4 + 664*sin(c + d*x)/(105*a**4*d) - 4*sin(c + d*x)/(a**4*d*(sec(c + d*x) + 1)) - 88*sin(c + d*x)/(105*a**4*d*(sec(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_79():
    f = cos(c + d*x)**2/(a*sec(c + d*x) + a)**4
    F = -sin(c + d*x)*cos(c + d*x)/(7*d*(a*sec(c + d*x) + a)**4) - 2*sin(c + d*x)*cos(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**3) + 21*x/(2*a**4) + 21*sin(c + d*x)*cos(c + d*x)/(2*a**4*d) - 576*sin(c + d*x)/(35*a**4*d) - 288*sin(c + d*x)*cos(c + d*x)/(35*a**4*d*(sec(c + d*x) + 1)) - 43*sin(c + d*x)*cos(c + d*x)/(35*a**4*d*(sec(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_80():
    f = sec(c + d*x)**7/(a*sec(c + d*x) + a)**5
    F = 5*tan(c + d*x)/(d*(a**5*sec(c + d*x) + a**5)) - tan(c + d*x)*sec(c + d*x)**5/(9*d*(a*sec(c + d*x) + a)**5) - 5*tan(c + d*x)*sec(c + d*x)**4/(21*a*d*(a*sec(c + d*x) + a)**4) - 29*tan(c + d*x)*sec(c + d*x)**3/(63*a**2*d*(a*sec(c + d*x) + a)**3) - 67*tan(c + d*x)*sec(c + d*x)**2/(63*a**3*d*(a*sec(c + d*x) + a)**2) + 181*tan(c + d*x)/(63*a**5*d) - 5*atanh(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_81():
    f = sec(c + d*x)**6/(a*sec(c + d*x) + a)**5
    F = -661*tan(c + d*x)/(315*d*(a**5*sec(c + d*x) + a**5)) - tan(c + d*x)*sec(c + d*x)**4/(9*d*(a*sec(c + d*x) + a)**5) - 13*tan(c + d*x)*sec(c + d*x)**3/(63*a*d*(a*sec(c + d*x) + a)**4) - 34*tan(c + d*x)*sec(c + d*x)**2/(105*a**2*d*(a*sec(c + d*x) + a)**3) + 173*tan(c + d*x)/(315*a**3*d*(a*sec(c + d*x) + a)**2) + atanh(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_82():
    f = sec(c + d*x)**5/(a*sec(c + d*x) + a)**5
    F = 4*tan(c + d*x)/(45*d*(a**5*sec(c + d*x) + a**5)) + tan(c + d*x)*sec(c + d*x)**4/(9*d*(a*sec(c + d*x) + a)**5) - 32*tan(c + d*x)/(315*a*d*(a**2*sec(c + d*x) + a**2)**2) + 4*tan(c + d*x)*sec(c + d*x)**3/(63*a*d*(a*sec(c + d*x) + a)**4) + 4*tan(c + d*x)/(105*a**2*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_83():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**5
    F = tan(c + d*x)/(9*d*(a**5*sec(c + d*x) + a**5)) - tan(c + d*x)*sec(c + d*x)**4/(9*d*(a*sec(c + d*x) + a)**5) - 8*tan(c + d*x)/(63*a*d*(a**2*sec(c + d*x) + a**2)**2) + 5*tan(c + d*x)*sec(c + d*x)**3/(63*a*d*(a*sec(c + d*x) + a)**4) + tan(c + d*x)/(21*a**2*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_84():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**5
    F = 2*tan(c + d*x)/(45*d*(a**5*sec(c + d*x) + a**5)) + tan(c + d*x)/(9*d*(a*sec(c + d*x) + a)**5) - 2*tan(c + d*x)/(9*a*d*(a*sec(c + d*x) + a)**4) + tan(c + d*x)/(15*a**2*d*(a*sec(c + d*x) + a)**3) + 2*tan(c + d*x)/(45*a**3*d*(a*sec(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_85():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**5
    F = 2*tan(c + d*x)/(63*d*(a**5*sec(c + d*x) + a**5)) - tan(c + d*x)/(9*d*(a*sec(c + d*x) + a)**5) + 2*tan(c + d*x)/(63*a*d*(a**2*sec(c + d*x) + a**2)**2) + 5*tan(c + d*x)/(63*a*d*(a*sec(c + d*x) + a)**4) + tan(c + d*x)/(21*a**2*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_86():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**5
    F = 8*tan(c + d*x)/(315*d*(a**5*sec(c + d*x) + a**5)) + tan(c + d*x)/(9*d*(a*sec(c + d*x) + a)**5) + 8*tan(c + d*x)/(315*a*d*(a**2*sec(c + d*x) + a**2)**2) + 4*tan(c + d*x)/(63*a*d*(a*sec(c + d*x) + a)**4) + 4*tan(c + d*x)/(105*a**2*d*(a*sec(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_87():
    f = (a*sec(c + d*x) + a)**(-5)
    F = -488*tan(c + d*x)/(315*d*(a**5*sec(c + d*x) + a**5)) - tan(c + d*x)/(9*d*(a*sec(c + d*x) + a)**5) - 13*tan(c + d*x)/(63*a*d*(a*sec(c + d*x) + a)**4) - 34*tan(c + d*x)/(105*a**2*d*(a*sec(c + d*x) + a)**3) - 173*tan(c + d*x)/(315*a**3*d*(a*sec(c + d*x) + a)**2) + x/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_88():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**5
    F = -5*sin(c + d*x)/(d*(a**5*sec(c + d*x) + a**5)) - sin(c + d*x)/(9*d*(a*sec(c + d*x) + a)**5) - 5*sin(c + d*x)/(21*a*d*(a*sec(c + d*x) + a)**4) - 29*sin(c + d*x)/(63*a**2*d*(a*sec(c + d*x) + a)**3) - 67*sin(c + d*x)/(63*a**3*d*(a*sec(c + d*x) + a)**2) - 5*x/a**5 + 496*sin(c + d*x)/(63*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_89():
    f = cos(c + d*x)**2/(a*sec(c + d*x) + a)**5
    F = -3832*sin(c + d*x)*cos(c + d*x)/(315*d*(a**5*sec(c + d*x) + a**5)) - sin(c + d*x)*cos(c + d*x)/(9*d*(a*sec(c + d*x) + a)**5) - 17*sin(c + d*x)*cos(c + d*x)/(63*a*d*(a*sec(c + d*x) + a)**4) - 28*sin(c + d*x)*cos(c + d*x)/(45*a**2*d*(a*sec(c + d*x) + a)**3) - 577*sin(c + d*x)*cos(c + d*x)/(315*a**3*d*(a*sec(c + d*x) + a)**2) + 31*x/(2*a**5) + 31*sin(c + d*x)*cos(c + d*x)/(2*a**5*d) - 7664*sin(c + d*x)/(315*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_90():
    f = sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**4
    F = 2*a*tan(c + d*x)*sec(c + d*x)**3/(7*d*sqrt(a*sec(c + d*x) + a)) + 4*a*tan(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a)) - 8*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(35*d) + 12*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_91():
    f = sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**3
    F = 14*a*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) - 4*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(15*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_92():
    f = sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**2
    F = 2*a*tan(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_93():
    f = sqrt(a*sec(c + d*x) + a)*sec(c + d*x)
    F = 2*a*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_94():
    f = sqrt(a*sec(c + d*x) + a)
    F = 2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_95():
    f = sqrt(a*sec(c + d*x) + a)*cos(c + d*x)
    F = sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_96():
    f = sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**2
    F = 3*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a*sin(c + d*x)*cos(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)) + 3*a*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_97():
    f = sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**3
    F = 5*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a*sin(c + d*x)*cos(c + d*x)**2/(3*d*sqrt(a*sec(c + d*x) + a)) + 5*a*sin(c + d*x)*cos(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)) + 5*a*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_98():
    f = sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**4
    F = 35*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + a*sin(c + d*x)*cos(c + d*x)**3/(4*d*sqrt(a*sec(c + d*x) + a)) + 7*a*sin(c + d*x)*cos(c + d*x)**2/(24*d*sqrt(a*sec(c + d*x) + a)) + 35*a*sin(c + d*x)*cos(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)) + 35*a*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_99():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4
    F = 2*a**2*tan(c + d*x)*sec(c + d*x)**4/(9*d*sqrt(a*sec(c + d*x) + a)) + 34*a**2*tan(c + d*x)*sec(c + d*x)**3/(63*d*sqrt(a*sec(c + d*x) + a)) + 68*a**2*tan(c + d*x)/(45*d*sqrt(a*sec(c + d*x) + a)) - 136*a*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(315*d) + 68*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_100():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3
    F = 152*a**2*tan(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)) + 38*a*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(105*d) - 4*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(35*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_101():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2
    F = 8*a**2*tan(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a)) + 2*a*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(5*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_102():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)
    F = 8*a**2*tan(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*a*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_103():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_104():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)
    F = 3*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a**2*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_105():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2
    F = 7*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a**2*sin(c + d*x)*cos(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)) + 7*a**2*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_106():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3
    F = 11*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a**2*sin(c + d*x)*cos(c + d*x)**2/(3*d*sqrt(a*sec(c + d*x) + a)) + 11*a**2*sin(c + d*x)*cos(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)) + 11*a**2*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_107():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4
    F = 46*a**3*tan(c + d*x)*sec(c + d*x)**4/(99*d*sqrt(a*sec(c + d*x) + a)) + 710*a**3*tan(c + d*x)*sec(c + d*x)**3/(693*d*sqrt(a*sec(c + d*x) + a)) + 284*a**3*tan(c + d*x)/(99*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**4/(11*d) - 568*a**2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(693*d) + 284*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(231*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_108():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**3
    F = 832*a**3*tan(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)) + 208*a**2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(315*d) + 26*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(105*d) - 4*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)/(63*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(7)/2)*tan(c + d*x)/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_109():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 64*a**3*tan(c + d*x)/(21*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(21*d) + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(7*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*tan(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_110():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)
    F = 64*a**3*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(15*d) + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_111():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 14*a**3*tan(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_112():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)
    F = 5*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d - a**3*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_113():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2
    F = 19*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + 9*a**3*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_114():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3
    F = 25*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + 13*a**3*sin(c + d*x)*cos(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)) + 25*a**3*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**2/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_115():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**4
    F = 163*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + 17*a**3*sin(c + d*x)*cos(c + d*x)**2/(24*d*sqrt(a*sec(c + d*x) + a)) + 163*a**3*sin(c + d*x)*cos(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)) + 163*a**3*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_116():
    f = sqrt(-a*sec(c + d*x) + a)*sec(c + d*x)
    F = -2*a*tan(c + d*x)/(d*sqrt(-a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_117():
    f = sqrt(-a*sec(c + d*x) + a)
    F = 2*sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_118():
    f = sqrt(-a*sec(c + d*x) + a)*cos(c + d*x)
    F = -sqrt(a)*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/d + a*sin(c + d*x)/(d*sqrt(-a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_119():
    f = sec(c + d*x)**4/sqrt(a*sec(c + d*x) + a)
    F = 2*tan(c + d*x)*sec(c + d*x)**2/(5*d*sqrt(a*sec(c + d*x) + a)) + 28*tan(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)) - 2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(15*a*d) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_120():
    f = sec(c + d*x)**3/sqrt(a*sec(c + d*x) + a)
    F = -4*tan(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(3*a*d) + sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_121():
    f = sec(c + d*x)**2/sqrt(a*sec(c + d*x) + a)
    F = 2*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_122():
    f = sec(c + d*x)/sqrt(a*sec(c + d*x) + a)
    F = sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_123():
    f = 1/sqrt(a*sec(c + d*x) + a)
    F = 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_124():
    f = cos(c + d*x)/sqrt(a*sec(c + d*x) + a)
    F = sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)) - atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_125():
    f = cos(c + d*x)**2/sqrt(a*sec(c + d*x) + a)
    F = sin(c + d*x)*cos(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)) - sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)) + 7*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*sqrt(a)*d) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_126():
    f = sec(c + d*x)**5/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)*sec(c + d*x)**3/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 9*tan(c + d*x)*sec(c + d*x)**2/(10*a*d*sqrt(a*sec(c + d*x) + a)) + 31*tan(c + d*x)/(5*a*d*sqrt(a*sec(c + d*x) + a)) - 13*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(10*a**2*d) - 15*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_127():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)*sec(c + d*x)**2/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 13*tan(c + d*x)/(3*a*d*sqrt(a*sec(c + d*x) + a)) + 7*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(6*a**2*d) + 11*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_128():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*tan(c + d*x)/(a*d*sqrt(a*sec(c + d*x) + a)) - 7*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_129():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_130():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_131():
    f = (a*sec(c + d*x) + a)**(sympy.S(-3)/2)
    F = -tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_132():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 3*sin(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)) - 3*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 9*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_133():
    f = cos(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*cos(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sin(c + d*x)*cos(c + d*x)/(a*d*sqrt(a*sec(c + d*x) + a)) - 7*sin(c + d*x)/(4*a*d*sqrt(a*sec(c + d*x) + a)) + 19*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*a**(sympy.S(3)/2)*d) - 13*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_134():
    f = sec(c + d*x)**5/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)*sec(c + d*x)**3/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 17*tan(c + d*x)*sec(c + d*x)**2/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) - 197*tan(c + d*x)/(24*a**2*d*sqrt(a*sec(c + d*x) + a)) + 95*sqrt(a*sec(c + d*x) + a)*tan(c + d*x)/(48*a**3*d) + 163*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_135():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)*sec(c + d*x)**2/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 13*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 9*tan(c + d*x)/(4*a**2*d*sqrt(a*sec(c + d*x) + a)) - 75*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_136():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = tan(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 13*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 19*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_137():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 5*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 5*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_138():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = tan(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 3*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_139():
    f = (a*sec(c + d*x) + a)**(sympy.S(-5)/2)
    F = -tan(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 11*tan(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 43*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_140():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 15*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 35*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) - 5*atan(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 115*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_141():
    f = sec(c + d*x)/sqrt(-a*sec(c + d*x) + a)
    F = -sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_142():
    f = 1/sqrt(-a*sec(c + d*x) + a)
    F = 2*atan(sqrt(a)*tan(c + d*x)/sqrt(-a*sec(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(c + d*x)/(2*sqrt(-a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_143():
    f = (a*sec(c + d*x) + a)**(sympy.S(2)/3)*sec(c + d*x)**3
    F = -19*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(160*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)) - 9*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(40*d) + 57*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(80*d*(sec(c + d*x) + 1)) + 3*(a*sec(c + d*x) + a)**(sympy.S(5)/3)*tan(c + d*x)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_144():
    f = (a*sec(c + d*x) + a)**(sympy.S(2)/3)*sec(c + d*x)**2
    F = -2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(10*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)) + 3*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(5*d) + 3*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(5*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_145():
    f = (a*sec(c + d*x) + a)**(sympy.S(2)/3)*sec(c + d*x)
    F = -2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(4*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)) + 3*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_146():
    f = (a*sec(c + d*x) + a)**(sympy.S(2)/3)
    F = 3*sqrt(2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S(7)/6, sympy.S.Half, 1, sympy.S(13)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(7*d*sqrt(1 - sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_147():
    f = (a*sec(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)
    F = -3*sqrt(2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S(7)/6, sympy.S.Half, 2, sympy.S(13)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(7*d*sqrt(1 - sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_148():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/3)*sec(c + d*x)**3
    F = -343*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*a*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(1760*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)) + 147*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(440*d) + 1029*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(880*d*(sec(c + d*x) + 1)) - 9*(a*sec(c + d*x) + a)**(sympy.S(5)/3)*tan(c + d*x)/(88*d) + 3*(a*sec(c + d*x) + a)**(sympy.S(8)/3)*tan(c + d*x)/(11*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_149():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/3)*sec(c + d*x)**2
    F = -7*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*a*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(32*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)) + 3*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(8*d) + 21*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(16*d*(sec(c + d*x) + 1)) + 3*(a*sec(c + d*x) + a)**(sympy.S(5)/3)*tan(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_150():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/3)*sec(c + d*x)
    F = -7*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*a*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(20*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)) + 3*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(5*d) + 21*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(10*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_151():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/3)
    F = 3*sqrt(2)*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(sec(c + d*x) + 1)*tan(c + d*x)*appellf1(sympy.S(13)/6, sympy.S.Half, 1, sympy.S(19)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(13*d*sqrt(1 - sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_152():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/3)*cos(c + d*x)
    F = -3*sqrt(2)*a*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(sec(c + d*x) + 1)*tan(c + d*x)*appellf1(sympy.S(13)/6, sympy.S.Half, 2, sympy.S(19)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(13*d*sqrt(1 - sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_153():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(1)/3)
    F = 37*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(160*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) + 3*tan(c + d*x)*sec(c + d*x)**2/(8*d*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) + 99*tan(c + d*x)/(80*d*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) - 3*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(40*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_154():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(1)/3)
    F = -7*2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(20*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) - 9*tan(c + d*x)/(10*d*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) + 3*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_155():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(1)/3)
    F = 2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(4*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)) + 3*tan(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_156():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(1)/3)
    F = -2**(sympy.S(2)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_157():
    f = (a*sec(c + d*x) + a)**(sympy.S(-1)/3)
    F = 3*sqrt(2)*tan(c + d*x)*appellf1(sympy.S(1)/6, sympy.S.Half, 1, sympy.S(7)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(d*sqrt(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_158():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(1)/3)
    F = -3*sqrt(2)*tan(c + d*x)*appellf1(sympy.S(1)/6, sympy.S.Half, 2, sympy.S(7)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(d*sqrt(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_159():
    f = sec(c + d*x)**4/(a*sec(c + d*x) + a)**(sympy.S(5)/3)
    F = 3*tan(c + d*x)*sec(c + d*x)**2/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/3)) - 33*tan(c + d*x)/(28*d*(a*sec(c + d*x) + a)**(sympy.S(5)/3)) + 135*tan(c + d*x)/(14*a*d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) - 375*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(28*a**2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) - 125*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(56*a**2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) + (375 + 375*sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)/(28*a**2*d*(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_160():
    f = sec(c + d*x)**3/(a*sec(c + d*x) + a)**(sympy.S(5)/3)
    F = 3*tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**(sympy.S(5)/3)) - 36*tan(c + d*x)/(7*a*d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) + 57*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(7*a**2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) + 19*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(14*a**2*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(sec(c + d*x) + 1)**(sympy.S(2)/3)) - (57 + 57*sqrt(3))*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)/(7*a**2*d*(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_161():
    f = sec(c + d*x)**2/(a*sec(c + d*x) + a)**(sympy.S(5)/3)
    F = -3*tan(c + d*x)/(7*d*(a*sec(c + d*x) + a)**(sympy.S(5)/3)) - 15*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)*elliptic_e(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(7*a*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) - 5*2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(14*a*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) + 15*tan(c + d*x)/(7*a*d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) + (15 + 15*sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)/(7*a*d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_162():
    f = sec(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/3)
    F = -6*2**(sympy.S(1)/3)*3**(sympy.S(1)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)*elliptic_e(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(7*a*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(1)/3)*3**(sympy.S(3)/4)*sqrt(((sec(c + d*x) + 1)**(sympy.S(2)/3) + 2**(sympy.S(1)/3)*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(2)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sqrt(3))*(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)*elliptic_f(acos((-(1 - sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))), sqrt(3)/4 + sympy.S.Half)/(7*a*d*sqrt(-(-(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)/(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3))**2)*(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) + 6*tan(c + d*x)/(7*a*d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)) + 3*tan(c + d*x)/(7*a*d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(sec(c + d*x) + 1)) + (6 + 6*sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3)*tan(c + d*x)/(7*a*d*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(-(1 + sqrt(3))*(sec(c + d*x) + 1)**(sympy.S(1)/3) + 2**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_163():
    f = (a*sec(c + d*x) + a)**(sympy.S(-5)/3)
    F = -3*sqrt(2)*tan(c + d*x)*appellf1(sympy.S(-7)/6, sympy.S.Half, 1, sympy.S(-1)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(7*a*d*sqrt(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_164():
    f = cos(c + d*x)/(a*sec(c + d*x) + a)**(sympy.S(5)/3)
    F = 3*sqrt(2)*tan(c + d*x)*appellf1(sympy.S(-7)/6, sympy.S.Half, 2, sympy.S(-1)/6, sec(c + d*x)/2 + sympy.S.Half, sec(c + d*x) + 1)/(7*a*d*sqrt(1 - sec(c + d*x))*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_165():
    f = (a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 6*a*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 6*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_166():
    f = (a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_167():
    f = (a*sec(c + d*x) + a)*sqrt(sec(c + d*x))
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_168():
    f = (a*sec(c + d*x) + a)/sqrt(sec(c + d*x))
    F = 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_169():
    f = (a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_170():
    f = (a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_171():
    f = (a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(7)/2)
    F = 10*a*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 6*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 10*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_172():
    f = (a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 4*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 8*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(7*d) + 12*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 12*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 8*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_173():
    f = (a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 4*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 16*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 16*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_174():
    f = (a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 4*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/d - 4*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 8*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_175():
    f = (a*sec(c + d*x) + a)**2/sqrt(sec(c + d*x))
    F = 2*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/d + 4*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_176():
    f = (a*sec(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 8*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_177():
    f = (a*sec(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(5)/2)
    F = 4*a**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 16*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_178():
    f = (a*sec(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(7)/2)
    F = 8*a**2*sin(c + d*x)/(7*d*sqrt(sec(c + d*x))) + 4*a**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 12*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 8*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_179():
    f = (a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 6*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 52*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 28*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 28*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 52*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_180():
    f = (a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x))
    F = 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/d + 36*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 36*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_181():
    f = (a*sec(c + d*x) + a)**3/sqrt(sec(c + d*x))
    F = 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 6*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/d - 4*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 20*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_182():
    f = (a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/d + 2*a**3*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 20*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_183():
    f = (a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**3*sin(c + d*x)/(d*sqrt(sec(c + d*x))) + 2*a**3*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 36*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_184():
    f = (a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(7)/2)
    F = 52*a**3*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 6*a**3*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**3*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 28*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 52*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_185():
    f = (a*sec(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(9)/2)
    F = 44*a**3*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 68*a**3*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a**3*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**3*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 68*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 44*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_186():
    f = (a*sec(c + d*x) + a)**4*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d) + 8*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 122*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(45*d) + 32*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(7*d) + 152*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) - 152*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 32*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_187():
    f = (a*sec(c + d*x) + a)**4*sqrt(sec(c + d*x))
    F = 2*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 8*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 94*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 64*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 64*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 136*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_188():
    f = (a*sec(c + d*x) + a)**4/sqrt(sec(c + d*x))
    F = 2*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 8*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 66*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 56*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 32*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_189():
    f = (a*sec(c + d*x) + a)**4/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 8*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/d + 2*a**4*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 40*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_190():
    f = (a*sec(c + d*x) + a)**4/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/d + 8*a**4*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a**4*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 56*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 32*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_191():
    f = (a*sec(c + d*x) + a)**4/sec(c + d*x)**(sympy.S(7)/2)
    F = 94*a**4*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 8*a**4*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**4*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 64*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 136*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_192():
    f = (a*sec(c + d*x) + a)**4/sec(c + d*x)**(sympy.S(9)/2)
    F = 32*a**4*sin(c + d*x)/(7*d*sqrt(sec(c + d*x))) + 122*a**4*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 8*a**4*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**4*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 152*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 32*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_193():
    f = (a*sec(c + d*x) + a)**4/sec(c + d*x)**(sympy.S(11)/2)
    F = 904*a**4*sin(c + d*x)/(231*d*sqrt(sec(c + d*x))) + 128*a**4*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 150*a**4*sin(c + d*x)/(77*d*sec(c + d*x)**(sympy.S(5)/2)) + 8*a**4*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a**4*sin(c + d*x)/(11*d*sec(c + d*x)**(sympy.S(9)/2)) + 128*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 904*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(231*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_194():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(d*(a*sec(c + d*x) + a)) + 5*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d) - 3*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d) + 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_195():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*(a*sec(c + d*x) + a)) + 3*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d) - 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) - sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_196():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) + sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_197():
    f = sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) - sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_198():
    f = 1/((a*sec(c + d*x) + a)*sqrt(sec(c + d*x)))
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) + 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) - sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_199():
    f = 1/((a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 5*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) - 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_200():
    f = 1/((a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - 5*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) + 7*sin(c + d*x)/(5*a*d*sec(c + d*x)**(sympy.S(3)/2)) + 21*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_201():
    f = sec(c + d*x)**(sympy.S(9)/2)/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(3*d*(a*sec(c + d*x) + a)**2) + 10*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d) - 7*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d) + 7*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 7*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_202():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*d*(a*sec(c + d*x) + a)**2) + 4*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d) - 4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 5*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_203():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) + sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_204():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**2
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_205():
    f = sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) - sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) + sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_206():
    f = 1/((a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x)))
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*(a*sec(c + d*x) + a)**2) + 4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 5*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_207():
    f = 1/((a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) + 10*sin(c + d*x)/(3*a**2*d*sqrt(sec(c + d*x))) - 7*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 7*sin(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_208():
    f = 1/((a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)) - 5*sin(c + d*x)/(a**2*d*sqrt(sec(c + d*x))) + 56*sin(c + d*x)/(15*a**2*d*sec(c + d*x)**(sympy.S(3)/2)) + 56*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a**2*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - 3*sin(c + d*x)/(a**2*d*(sec(c + d*x) + 1)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_209():
    f = sec(c + d*x)**(sympy.S(11)/2)/(a*sec(c + d*x) + a)**3
    F = -119*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(30*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(3*a*d*(a*sec(c + d*x) + a)**2) + 11*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*a**3*d) - 119*sin(c + d*x)*sqrt(sec(c + d*x))/(10*a**3*d) + 119*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + 11*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_210():
    f = sec(c + d*x)**(sympy.S(9)/2)/(a*sec(c + d*x) + a)**3
    F = -13*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 8*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(15*a*d*(a*sec(c + d*x) + a)**2) + 49*sin(c + d*x)*sqrt(sec(c + d*x))/(10*a**3*d) - 49*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) - 13*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_211():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**3
    F = -9*sin(c + d*x)*sqrt(sec(c + d*x))/(10*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*a*d*(a*sec(c + d*x) + a)**2) + 9*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_212():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**3
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 4*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_213():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**3
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) + sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) - sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) - sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_214():
    f = sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)**3
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(2*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) + 2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*a*d*(a*sec(c + d*x) + a)**2) - 9*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_215():
    f = 1/((a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x)))
    F = -13*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*(a*sec(c + d*x) + a)**3) - 8*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + 49*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) - 13*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_216():
    f = 1/((a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2))
    F = -119*sin(c + d*x)/(30*d*(a**3*sec(c + d*x) + a**3)*sqrt(sec(c + d*x))) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x))) - 2*sin(c + d*x)/(3*a*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) + 11*sin(c + d*x)/(2*a**3*d*sqrt(sec(c + d*x))) - 119*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + 11*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_217():
    f = 1/((a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(5)/2))
    F = -63*sin(c + d*x)/(10*d*(a**3*sec(c + d*x) + a**3)*sec(c + d*x)**(sympy.S(3)/2)) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2)) - 4*sin(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)) - 21*sin(c + d*x)/(2*a**3*d*sqrt(sec(c + d*x))) + 77*sin(c + d*x)/(10*a**3*d*sec(c + d*x)**(sympy.S(3)/2)) + 231*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) - 21*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_218():
    f = sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = 3*sqrt(a)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*sqrt(a*sec(c + d*x) + a)) + 3*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_219():
    f = sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = sqrt(a)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_220():
    f = sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))
    F = 2*sqrt(a)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_221():
    f = sqrt(a*sec(c + d*x) + a)/sqrt(sec(c + d*x))
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_222():
    f = sqrt(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = 4*a*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*a*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_223():
    f = sqrt(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(5)/2)
    F = 16*a*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 8*a*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*a*sin(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_224():
    f = sqrt(a*sec(c + d*x) + a)/sec(c + d*x)**(sympy.S(7)/2)
    F = 32*a*sin(c + d*x)*sqrt(sec(c + d*x))/(35*d*sqrt(a*sec(c + d*x) + a)) + 16*a*sin(c + d*x)/(35*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 12*a*sin(c + d*x)/(35*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(7*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_225():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = 11*a**(sympy.S(3)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(3*d*sqrt(a*sec(c + d*x) + a)) + 11*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(12*d*sqrt(a*sec(c + d*x) + a)) + 11*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(8*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_226():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = 7*a**(sympy.S(3)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*sqrt(a*sec(c + d*x) + a)) + 7*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_227():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = 3*a**(sympy.S(3)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_228():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = 2*a**(sympy.S(3)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 2*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_229():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = 8*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_230():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = 8*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*sqrt(a*sec(c + d*x) + a)) + 2*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(5*d*sqrt(sec(c + d*x))) + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_231():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 208*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + 104*a**2*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 26*a**2*sin(c + d*x)/(35*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*sin(c + d*x)/(7*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_232():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = 544*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 272*a**2*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 68*a**2*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 34*a**2*sin(c + d*x)/(63*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**2*sin(c + d*x)/(9*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_233():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = 163*a**(sympy.S(5)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(64*d) + 17*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(24*d*sqrt(a*sec(c + d*x) + a)) + 163*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(96*d*sqrt(a*sec(c + d*x) + a)) + 163*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(64*d*sqrt(a*sec(c + d*x) + a)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_234():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = 25*a**(sympy.S(5)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(8*d) + 13*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(12*d*sqrt(a*sec(c + d*x) + a)) + 25*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(8*d*sqrt(a*sec(c + d*x) + a)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_235():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = 19*a**(sympy.S(5)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(4*d) + 9*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(a*sec(c + d*x) + a)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_236():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = 5*a**(sympy.S(5)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*sec(c + d*x) + a)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_237():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(5)/2)*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/d + 14*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_238():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = 64*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x))) + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_239():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 64*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(21*d*sqrt(a*sec(c + d*x) + a)) + 16*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_240():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = 1168*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 584*a**3*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 146*a**3*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 38*a**3*sin(c + d*x)/(63*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_241():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(11)/2)
    F = 2272*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(693*d*sqrt(a*sec(c + d*x) + a)) + 1136*a**3*sin(c + d*x)/(693*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 284*a**3*sin(c + d*x)/(231*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 710*a**3*sin(c + d*x)/(693*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + 46*a**3*sin(c + d*x)/(99*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2)) + 2*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(11*d*sec(c + d*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_242():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(1)/4)
    F = 4*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/4)/(d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_243():
    f = sqrt(a*sec(e + f*x) + a)*sqrt(sec(e + f*x))
    F = 2*sqrt(a)*asinh(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_244():
    f = sqrt(-sec(e + f*x))*sqrt(-a*sec(e + f*x) + a)
    F = 2*sqrt(a)*asinh(sqrt(a)*tan(e + f*x)/sqrt(-a*sec(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_245():
    f = sec(c + d*x)**(sympy.S(5)/2)/sqrt(a*sec(c + d*x) + a)
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*sqrt(a*sec(c + d*x) + a)) - asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_246():
    f = sec(c + d*x)**(sympy.S(3)/2)/sqrt(a*sec(c + d*x) + a)
    F = 2*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_247():
    f = sqrt(sec(c + d*x))/sqrt(a*sec(c + d*x) + a)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_248():
    f = 1/(sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x)))
    F = 2*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*sec(c + d*x) + a)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_249():
    f = 1/(sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = -2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) + 2*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_250():
    f = 1/(sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    F = 26*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) - 2*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*sin(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_251():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*a*d*sqrt(a*sec(c + d*x) + a)) - 3*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 9*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_252():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_253():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_254():
    f = sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_255():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 5*sin(c + d*x)*sqrt(sec(c + d*x))/(2*a*d*sqrt(a*sec(c + d*x) + a)) - 7*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_256():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) - 19*sin(c + d*x)*sqrt(sec(c + d*x))/(6*a*d*sqrt(a*sec(c + d*x) + a)) + 7*sin(c + d*x)/(6*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 11*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_257():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) + 49*sin(c + d*x)*sqrt(sec(c + d*x))/(10*a*d*sqrt(a*sec(c + d*x) + a)) - 13*sin(c + d*x)/(10*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 9*sin(c + d*x)/(10*a*d*sqrt(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - 15*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_258():
    f = sec(c + d*x)**(sympy.S(9)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 15*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 35*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) - 5*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 115*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_259():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 11*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 2*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 43*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_260():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_261():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) + 5*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_262():
    f = sqrt(sec(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 9*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 19*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_263():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 13*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 49*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a**2*d*sqrt(a*sec(c + d*x) + a)) - 75*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_264():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - 17*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) - 299*sin(c + d*x)*sqrt(sec(c + d*x))/(48*a**2*d*sqrt(a*sec(c + d*x) + a)) + 95*sin(c + d*x)/(48*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 163*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_265():
    f = sec(c + d*x)**(sympy.S(7)/2)/sqrt(sec(c + d*x) + 1)
    F = -sqrt(2)*asinh(tan(c + d*x)/(sec(c + d*x) + 1))/d + 7*asinh(tan(c + d*x)/sqrt(sec(c + d*x) + 1))/(4*d) + sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*sqrt(sec(c + d*x) + 1)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_266():
    f = sec(c + d*x)**(sympy.S(5)/2)/sqrt(sec(c + d*x) + 1)
    F = sqrt(2)*asinh(tan(c + d*x)/(sec(c + d*x) + 1))/d - asinh(tan(c + d*x)/sqrt(sec(c + d*x) + 1))/d + sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_267():
    f = sec(c + d*x)**(sympy.S(3)/2)/sqrt(sec(c + d*x) + 1)
    F = -sqrt(2)*asinh(tan(c + d*x)/(sec(c + d*x) + 1))/d + 2*asinh(tan(c + d*x)/sqrt(sec(c + d*x) + 1))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_268():
    f = sqrt(sec(c + d*x))/sqrt(sec(c + d*x) + 1)
    F = sqrt(2)*asinh(tan(c + d*x)/(sec(c + d*x) + 1))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_269():
    f = 1/(sqrt(sec(c + d*x) + 1)*sqrt(sec(c + d*x)))
    F = -sqrt(2)*asinh(tan(c + d*x)/(sec(c + d*x) + 1))/d + 2*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_270():
    f = 1/(sqrt(sec(c + d*x) + 1)*sec(c + d*x)**(sympy.S(3)/2))
    F = sqrt(2)*asinh(tan(c + d*x)/(sec(c + d*x) + 1))/d - 2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(sec(c + d*x) + 1)) + 2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x) + 1)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_271():
    f = 1/(sqrt(sec(c + d*x) + 1)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(2)*asinh(tan(c + d*x)/(sec(c + d*x) + 1))/d + 26*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(sec(c + d*x) + 1)) - 2*sin(c + d*x)/(15*d*sqrt(sec(c + d*x) + 1)*sqrt(sec(c + d*x))) + 2*sin(c + d*x)/(5*d*sqrt(sec(c + d*x) + 1)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_272():
    f = (e*sec(c + d*x))**(sympy.S(4)/3)*sqrt(a*sec(c + d*x) + a)
    F = 4*3**(sympy.S(3)/4)*a**2*e*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(5*d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) + 6*a*e*(e*sec(c + d*x))**(sympy.S(1)/3)*tan(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_273():
    f = (e*sec(c + d*x))**(sympy.S(1)/3)*sqrt(a*sec(c + d*x) + a)
    F = 2*3**(sympy.S(3)/4)*a**2*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_274():
    f = sqrt(a*sec(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(2)/3)
    F = 3**(sympy.S(3)/4)*a**2*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(sqrt(3) + 2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(2*d*e*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) + 3*a*tan(c + d*x)/(2*d*(e*sec(c + d*x))**(sympy.S(2)/3)*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_275():
    f = (e*sec(c + d*x))**(sympy.S(8)/3)*sqrt(a*sec(c + d*x) + a)
    F = 120*3**(sympy.S(1)/4)*a**2*e**(sympy.S(7)/3)*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(91*d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) - 80*sqrt(2)*3**(sympy.S(3)/4)*a**2*e**(sympy.S(7)/3)*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(91*d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) - 240*a*e**3*tan(c + d*x)/(91*d*sqrt(a*sec(c + d*x) + a)*(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))) + 60*a*e**2*(e*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)/(91*d*sqrt(a*sec(c + d*x) + a)) + 6*a*e*(e*sec(c + d*x))**(sympy.S(5)/3)*tan(c + d*x)/(13*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_276():
    f = (e*sec(c + d*x))**(sympy.S(5)/3)*sqrt(a*sec(c + d*x) + a)
    F = 12*3**(sympy.S(1)/4)*a**2*e**(sympy.S(4)/3)*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(7*d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) - 8*sqrt(2)*3**(sympy.S(3)/4)*a**2*e**(sympy.S(4)/3)*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(7*d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) - 24*a*e**2*tan(c + d*x)/(7*d*sqrt(a*sec(c + d*x) + a)*(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))) + 6*a*e*(e*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)/(7*d*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_277():
    f = (e*sec(c + d*x))**(sympy.S(2)/3)*sqrt(a*sec(c + d*x) + a)
    F = 3*3**(sympy.S(1)/4)*a**2*e**(sympy.S(1)/3)*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) - 2*sqrt(2)*3**(sympy.S(3)/4)*a**2*e**(sympy.S(1)/3)*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(d*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) - 6*a*e*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_278():
    f = sqrt(a*sec(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*3**(sympy.S(1)/4)*a**2*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(2*d*e**(sympy.S(2)/3)*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) + sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(d*e**(sympy.S(2)/3)*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) + 3*a*tan(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))) + 3*a*tan(c + d*x)/(d*(e*sec(c + d*x))**(sympy.S(1)/3)*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_279():
    f = sqrt(a*sec(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(4)/3)
    F = -15*3**(sympy.S(1)/4)*a**2*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*sqrt(2 - sqrt(3))*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_e(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(16*d*e**(sympy.S(5)/3)*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) + 5*sqrt(2)*3**(sympy.S(3)/4)*a**2*sqrt((e**(sympy.S(2)/3) + e**(sympy.S(1)/3)*(e*sec(c + d*x))**(sympy.S(1)/3) + (e*sec(c + d*x))**(sympy.S(2)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))*tan(c + d*x)*elliptic_f(asin((e**(sympy.S(1)/3)*(1 - sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))), -7 - 4*sqrt(3))/(8*d*e**(sympy.S(5)/3)*sqrt(e**(sympy.S(1)/3)*(e**(sympy.S(1)/3) - (e*sec(c + d*x))**(sympy.S(1)/3))/(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))**2)*(-a*sec(c + d*x) + a)*sqrt(a*sec(c + d*x) + a)) + 3*a*tan(c + d*x)/(4*d*(e*sec(c + d*x))**(sympy.S(4)/3)*sqrt(a*sec(c + d*x) + a)) + 15*a*tan(c + d*x)/(8*d*e*sqrt(a*sec(c + d*x) + a)*(e**(sympy.S(1)/3)*(1 + sqrt(3)) - (e*sec(c + d*x))**(sympy.S(1)/3))) + 15*a*tan(c + d*x)/(8*d*e*(e*sec(c + d*x))**(sympy.S(1)/3)*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_280():
    f = (e*sec(c + d*x))**(sympy.S(2)/3)/sqrt(a*sec(c + d*x) + a)
    F = -3*(e*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S(2)/3, sympy.S.Half, 1, sympy.S(5)/3, sec(c + d*x), -sec(c + d*x))/(2*d*sqrt(1 - sec(c + d*x))*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_281():
    f = (e*sec(c + d*x))**(sympy.S(1)/3)/sqrt(a*sec(c + d*x) + a)
    F = -3*(e*sec(c + d*x))**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S(1)/3, sympy.S.Half, 1, sympy.S(4)/3, sec(c + d*x), -sec(c + d*x))/(d*sqrt(1 - sec(c + d*x))*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_282():
    f = 1/((e*sec(c + d*x))**(sympy.S(1)/3)*sqrt(a*sec(c + d*x) + a))
    F = 3*tan(c + d*x)*appellf1(sympy.S(-1)/3, sympy.S.Half, 1, sympy.S(2)/3, sec(c + d*x), -sec(c + d*x))/(d*(e*sec(c + d*x))**(sympy.S(1)/3)*sqrt(1 - sec(c + d*x))*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_283():
    f = 1/((e*sec(c + d*x))**(sympy.S(2)/3)*sqrt(a*sec(c + d*x) + a))
    F = 3*tan(c + d*x)*appellf1(sympy.S(-2)/3, sympy.S.Half, 1, sympy.S(1)/3, sec(c + d*x), -sec(c + d*x))/(2*d*(e*sec(c + d*x))**(sympy.S(2)/3)*sqrt(1 - sec(c + d*x))*sqrt(a*sec(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_284():
    f = (a*sec(c + d*x) + a)**(sympy.S(1)/3)*sec(c + d*x)**(sympy.S(4)/3)
    F = 2**(sympy.S(5)/6)*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S(1)/6, sympy.S(3)/2, 1 - sec(c + d*x), sympy.S.Half - sec(c + d*x)/2)/(d*(sec(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_285():
    f = (a*sec(c + d*x) + a)**(sympy.S(2)/3)*sec(c + d*x)**(sympy.S(4)/3)
    F = 2*2**(sympy.S(1)/6)*(a*sec(c + d*x) + a)**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S(-1)/6, sympy.S(3)/2, 1 - sec(c + d*x), sympy.S.Half - sec(c + d*x)/2)/(d*(sec(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_286():
    f = (a*sec(c + d*x) + a)**(sympy.S(2)/3)*sec(c + d*x)**(sympy.S(5)/3)
    F = -3*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/3)/(2*d*(a*(sec(c + d*x) + 1))**(sympy.S(1)/3)) + (a*(sec(c + d*x) + 1))**(sympy.S(2)/3)*(cos(c + d*x)*sec(c/2 + d*x/2)**4)**(sympy.S(1)/3)*tan(c + d*x)*hyper((sympy.S(1)/4, sympy.S(1)/3), (sympy.S(5)/4,), tan(c/2 + d*x/2)**4)/(8*d*(sec(c + d*x) + 1)**(sympy.S(4)/3)*(1/(cos(c + d*x) + 1))**(sympy.S(1)/3)) - 5*(a*(sec(c + d*x) + 1))**(sympy.S(2)/3)*(cos(c + d*x)*sec(c/2 + d*x/2)**4)**(sympy.S(1)/3)*tan(c + d*x)**3*hyper((sympy.S(1)/3, sympy.S(3)/4), (sympy.S(7)/4,), tan(c/2 + d*x/2)**4)/(8*d*(sec(c + d*x) + 1)**(sympy.S(10)/3)*(1/(cos(c + d*x) + 1))**(sympy.S(1)/3)) + 9*(a*(sec(c + d*x) + 1))**(sympy.S(2)/3)*sin(c + d*x)*sec(c + d*x)**(sympy.S(2)/3)/(4*d) - 9*(a*(sec(c + d*x) + 1))**(sympy.S(2)/3)*tan(c + d*x)/(4*d*(sec(c + d*x) + 1)**(sympy.S(7)/3)*(1/(cos(c + d*x) + 1))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_287():
    f = (a*sec(c + d*x) + a)**(sympy.S(4)/3)/sec(c + d*x)**(sympy.S(1)/3)
    F = 2*2**(sympy.S(5)/6)*a*(a*sec(c + d*x) + a)**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/6, sympy.S(4)/3, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, 1 - sec(c + d*x))/(d*(sec(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_288():
    f = (a*sec(e + f*x) + a)**4*sec(e + f*x)**n
    F = a**4*(4*n**2 + 21*n + 30)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(n + 1)*(n + 2)*(n + 3)) - a**4*(8*n**2 + 24*n + 3)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(f*(1 - n)*(n + 1)*(n + 3)*sqrt(sin(e + f*x)**2)) + 4*a**4*(2*n + 3)*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**n/(f*n*(n + 2)*sqrt(sin(e + f*x)**2)) + (a**2*sec(e + f*x) + a**2)**2*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(n + 3)) + (2*n + 8)*(a**4*sec(e + f*x) + a**4)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(n + 2)*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_289():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)**n
    F = a**3*(2*n + 5)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(n + 1)*(n + 2)) - a**3*(4*n + 1)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(f*(1 - n**2)*sqrt(sin(e + f*x)**2)) + a**3*(4*n + 7)*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**n/(f*n*(n + 2)*sqrt(sin(e + f*x)**2)) + (a**3*sec(e + f*x) + a**3)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_290():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)**n
    F = a**2*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(n + 1)) - a**2*(2*n + 1)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(f*(1 - n**2)*sqrt(sin(e + f*x)**2)) + 2*a**2*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**n/(f*n*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_291():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)**n
    F = -a*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(f*(1 - n)*sqrt(sin(e + f*x)**2)) + a*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**n/(f*n*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_292():
    f = sec(e + f*x)**n/(a*sec(e + f*x) + a)
    F = sin(e + f*x)*sec(e + f*x)**n/(f*(a*sec(e + f*x) + a)) + (1 - n)*sin(e + f*x)*hyper((sympy.S.Half, 1 - n/2), (2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 2)/(a*f*(2 - n)*sqrt(sin(e + f*x)**2)) - sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(a*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_293():
    f = sec(e + f*x)**n/(a*sec(e + f*x) + a)**2
    F = -sin(e + f*x)*sec(e + f*x)**(n + 1)/(3*f*(a*sec(e + f*x) + a)**2) - (3 - 2*n)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**(n - 1)/(3*a**2*f*sqrt(sin(e + f*x)**2)) + (4 - 2*n)*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)*sec(e + f*x)**n/(3*a**2*f*sqrt(sin(e + f*x)**2)) - (4 - 2*n)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(3*a**2*f*(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_294():
    f = (sec(e + f*x) + 1)**(sympy.S(5)/2)*sec(e + f*x)**n
    F = 2*sqrt(sec(e + f*x) + 1)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(2*n + 3)) + (8*n + 14)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(2*n + 1)*(2*n + 3)*sqrt(sec(e + f*x) + 1)) + (32*n**2 + 48*n + 6)*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*(2*n + 3)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_295():
    f = (sec(e + f*x) + 1)**(sympy.S(3)/2)*sec(e + f*x)**n
    F = (8*n + 2)*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*sqrt(sec(e + f*x) + 1)) + 2*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(2*n + 1)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_296():
    f = sqrt(sec(e + f*x) + 1)*sec(e + f*x)**n
    F = 2*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))/(f*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_297():
    f = sec(e + f*x)**n/sqrt(sec(e + f*x) + 1)
    F = tan(e + f*x)*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/(f*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_298():
    f = sec(e + f*x)**n/(sec(e + f*x) + 1)**(sympy.S(3)/2)
    F = tan(e + f*x)*appellf1(sympy.S.Half, 2, 1 - n, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/(2*f*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_299():
    f = (-sec(e + f*x))**n*(sec(e + f*x) + 1)**(sympy.S(3)/2)
    F = 2*(-sec(e + f*x))**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(sec(e + f*x) + 1)) - (-sec(e + f*x))**n*(4*n + 1)*tan(e + f*x)*hyper((sympy.S.Half, n), (n + 1,), sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*(2*n + 1)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_300():
    f = (-sec(e + f*x))**n*sqrt(sec(e + f*x) + 1)
    F = -(-sec(e + f*x))**n*tan(e + f*x)*hyper((sympy.S.Half, n), (n + 1,), sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_301():
    f = (-sec(e + f*x))**n/sqrt(sec(e + f*x) + 1)
    F = -(-sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 1, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_302():
    f = (-sec(e + f*x))**n/(sec(e + f*x) + 1)**(sympy.S(3)/2)
    F = -(-sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 2, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_303():
    f = (d*sec(e + f*x))**n*(sec(e + f*x) + 1)**(sympy.S(3)/2)
    F = 2*(d*sec(e + f*x))**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(sec(e + f*x) + 1)) - (d*sec(e + f*x))**n*(4*n + 1)*tan(e + f*x)*hyper((sympy.S.Half, n), (n + 1,), sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*(2*n + 1)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_304():
    f = (d*sec(e + f*x))**n*sqrt(sec(e + f*x) + 1)
    F = -(d*sec(e + f*x))**n*tan(e + f*x)*hyper((sympy.S.Half, n), (n + 1,), sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_305():
    f = (d*sec(e + f*x))**n/sqrt(sec(e + f*x) + 1)
    F = -(d*sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 1, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_306():
    f = (d*sec(e + f*x))**n/(sec(e + f*x) + 1)**(sympy.S(3)/2)
    F = -(d*sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 2, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_307():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)**n
    F = 2*a**3*(4*n + 7)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(2*n + 1)*(2*n + 3)*sqrt(a*sec(e + f*x) + a)) + 2*a**3*(16*n**2 + 24*n + 3)*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*(2*n + 3)*sqrt(a*sec(e + f*x) + a)) + 2*a**2*sqrt(a*sec(e + f*x) + a)*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(2*n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_308():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)**n
    F = 2*a**2*(4*n + 1)*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) + 2*a**2*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_309():
    f = sqrt(a*sec(e + f*x) + a)*sec(e + f*x)**n
    F = 2*a*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_310():
    f = sec(e + f*x)**n/sqrt(a*sec(e + f*x) + a)
    F = tan(e + f*x)*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_311():
    f = sec(e + f*x)**n/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = tan(e + f*x)*appellf1(sympy.S.Half, 2, 1 - n, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/(2*a*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_312():
    f = (-sec(e + f*x))**n*(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*a**2*(-sec(e + f*x))**n*(4*n + 1)*sin(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))*sec(e + f*x)**(1 - n)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(-sec(e + f*x))**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_313():
    f = (-sec(e + f*x))**n*sqrt(a*sec(e + f*x) + a)
    F = 2*a*(-sec(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))*sec(e + f*x)**(1 - n)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_314():
    f = (-sec(e + f*x))**n/sqrt(a*sec(e + f*x) + a)
    F = -(-sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 1, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_315():
    f = (-sec(e + f*x))**n/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -(-sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 2, n + 1, sec(e + f*x), -sec(e + f*x))/(a*f*n*sqrt(1 - sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_316():
    f = (d*sec(e + f*x))**n*(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*a**2*(d*sec(e + f*x))**n*(4*n + 1)*sin(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))*sec(e + f*x)**(1 - n)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(d*sec(e + f*x))**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_317():
    f = (d*sec(e + f*x))**n*sqrt(a*sec(e + f*x) + a)
    F = 2*a*(d*sec(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), 1 - sec(e + f*x))*sec(e + f*x)**(1 - n)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_318():
    f = (d*sec(e + f*x))**n/sqrt(a*sec(e + f*x) + a)
    F = -(d*sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 1, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_319():
    f = (d*sec(e + f*x))**n/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -(d*sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, 2, n + 1, sec(e + f*x), -sec(e + f*x))/(a*f*n*sqrt(1 - sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_320():
    f = (-sec(e + f*x))**n*(-a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = 2*a**3*(-sec(e + f*x))**n*(4*n + 7)*tan(e + f*x)/(f*(2*n + 1)*(2*n + 3)*sqrt(-a*sec(e + f*x) + a)) + 2*a**3*(16*n**2 + 24*n + 3)*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), sec(e + f*x) + 1)/(f*(2*n + 1)*(2*n + 3)*sqrt(-a*sec(e + f*x) + a)) + 2*a**2*(-sec(e + f*x))**n*sqrt(-a*sec(e + f*x) + a)*tan(e + f*x)/(f*(2*n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_321():
    f = (-sec(e + f*x))**n*(-a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*a**2*(-sec(e + f*x))**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(-a*sec(e + f*x) + a)) + 2*a**2*(4*n + 1)*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), sec(e + f*x) + 1)/(f*(2*n + 1)*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_322():
    f = (-sec(e + f*x))**n*sqrt(-a*sec(e + f*x) + a)
    F = 2*a*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), sec(e + f*x) + 1)/(f*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_323():
    f = (-sec(e + f*x))**n/sqrt(-a*sec(e + f*x) + a)
    F = tan(e + f*x)*appellf1(sympy.S.Half, 1, 1 - n, sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/(f*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_324():
    f = (-sec(e + f*x))**n/(-a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = tan(e + f*x)*appellf1(sympy.S.Half, 2, 1 - n, sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/(2*a*f*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_325():
    f = (-a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)**n
    F = 2*a**2*sin(e + f*x)*sec(e + f*x)**(n + 1)/(f*(2*n + 1)*sqrt(-a*sec(e + f*x) + a)) + 2*a**2*(4*n + 1)*sin(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), sec(e + f*x) + 1)*sec(e + f*x)**(n + 1)/(f*(-sec(e + f*x))**n*(2*n + 1)*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_326():
    f = sqrt(-a*sec(e + f*x) + a)*sec(e + f*x)**n
    F = 2*a*sin(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), sec(e + f*x) + 1)*sec(e + f*x)**(n + 1)/(f*(-sec(e + f*x))**n*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_327():
    f = (d*sec(e + f*x))**n*(-a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*a**2*(d*sec(e + f*x))**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(-a*sec(e + f*x) + a)) + 2*a**2*(d*sec(e + f*x))**n*(4*n + 1)*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), sec(e + f*x) + 1)/(f*(-sec(e + f*x))**n*(2*n + 1)*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_328():
    f = (d*sec(e + f*x))**n*sqrt(-a*sec(e + f*x) + a)
    F = 2*a*(d*sec(e + f*x))**n*tan(e + f*x)*hyper((sympy.S.Half, 1 - n), (sympy.S(3)/2,), sec(e + f*x) + 1)/(f*(-sec(e + f*x))**n*sqrt(-a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_329():
    f = (sec(e + f*x) + 1)**m*sec(e + f*x)**n
    F = 2**(m + sympy.S.Half)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half - m, 1 - n, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/(f*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_330():
    f = (1 - sec(e + f*x))**m*sec(e + f*x)**n
    F = sqrt(2)*(1 - sec(e + f*x))**m*tan(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1 - n, m + sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/(f*(2*m + 1)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_331():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)**n
    F = 2**(m + sympy.S.Half)*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half - m, 1 - n, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_332():
    f = (-a*sec(e + f*x) + a)**m*sec(e + f*x)**n
    F = sqrt(2)*(-a*sec(e + f*x) + a)**m*tan(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1 - n, m + sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, 1 - sec(e + f*x))/(f*(2*m + 1)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_333():
    f = (-sec(e + f*x))**n*(sec(e + f*x) + 1)**m
    F = sqrt(2)*(sec(e + f*x) + 1)**m*tan(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1 - n, m + sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/(f*sqrt(1 - sec(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_334():
    f = (-sec(e + f*x))**n*(1 - sec(e + f*x))**m
    F = 2**(m + sympy.S.Half)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half - m, 1 - n, sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/(f*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_335():
    f = (-sec(e + f*x))**n*(a*sec(e + f*x) + a)**m
    F = sqrt(2)*(a*sec(e + f*x) + a)**m*tan(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1 - n, m + sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/(f*sqrt(1 - sec(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_336():
    f = (-sec(e + f*x))**n*(-a*sec(e + f*x) + a)**m
    F = 2**(m + sympy.S.Half)*(1 - sec(e + f*x))**(-m + sympy.S(-1)/2)*(-a*sec(e + f*x) + a)**m*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half - m, 1 - n, sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_337():
    f = (d*sec(e + f*x))**n*(sec(e + f*x) + 1)**m
    F = -(d*sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, sympy.S.Half - m, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_338():
    f = (d*sec(e + f*x))**n*(1 - sec(e + f*x))**m
    F = -(d*sec(e + f*x))**n*tan(e + f*x)*appellf1(n, sympy.S.Half, sympy.S.Half - m, n + 1, -sec(e + f*x), sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x))*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_339():
    f = (d*sec(e + f*x))**n*(a*sec(e + f*x) + a)**m
    F = -(d*sec(e + f*x))**n*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(n, sympy.S.Half, sympy.S.Half - m, n + 1, sec(e + f*x), -sec(e + f*x))/(f*n*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_340():
    f = (d*sec(e + f*x))**n*(-a*sec(e + f*x) + a)**m
    F = -(d*sec(e + f*x))**n*(1 - sec(e + f*x))**(-m + sympy.S(-1)/2)*(-a*sec(e + f*x) + a)**m*tan(e + f*x)*appellf1(n, sympy.S.Half, sympy.S.Half - m, n + 1, -sec(e + f*x), sec(e + f*x))/(f*n*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_341():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)**4
    F = 2**(m + sympy.S.Half)*m*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(m**2 + 3*m + 5)*tan(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sec(e + f*x)/2)/(f*(m + 1)*(m + 2)*(m + 3)) + (a*sec(e + f*x) + a)**m*tan(e + f*x)*sec(e + f*x)**2/(f*(m + 3)) + (m + 4)*(a*sec(e + f*x) + a)**m*tan(e + f*x)/(f*(m + 1)*(m + 2)*(m + 3)) + m*(a*sec(e + f*x) + a)**(m + 1)*tan(e + f*x)/(a*f*(m**2 + 5*m + 6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_342():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)**3
    F = 2**(m + sympy.S.Half)*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*(m**2 + m + 1)*tan(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sec(e + f*x)/2)/(f*(m + 1)*(m + 2)) - (a*sec(e + f*x) + a)**m*tan(e + f*x)/(f*(m**2 + 3*m + 2)) + (a*sec(e + f*x) + a)**(m + 1)*tan(e + f*x)/(a*f*(m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_343():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)**2
    F = 2**(m + sympy.S.Half)*m*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sec(e + f*x)/2)/(f*(m + 1)) + (a*sec(e + f*x) + a)**m*tan(e + f*x)/(f*(m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_344():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)
    F = 2**(m + sympy.S.Half)*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m), (sympy.S(3)/2,), sympy.S.Half - sec(e + f*x)/2)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_345():
    f = (a*sec(e + f*x) + a)**m
    F = sqrt(2)*(a*sec(e + f*x) + a)**m*tan(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 1, m + sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/(f*sqrt(1 - sec(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_346():
    f = (a*sec(e + f*x) + a)**m*cos(e + f*x)
    F = -sqrt(2)*(a*sec(e + f*x) + a)**m*tan(e + f*x)*appellf1(m + sympy.S.Half, sympy.S.Half, 2, m + sympy.S(3)/2, sec(e + f*x)/2 + sympy.S.Half, sec(e + f*x) + 1)/(f*sqrt(1 - sec(e + f*x))*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_347():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)*(a*sec(e + f*x) + a)**m
    F = -2*(d*sec(e + f*x))**(sympy.S(3)/2)*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(sympy.S(3)/2, sympy.S.Half, sympy.S.Half - m, sympy.S(5)/2, sec(e + f*x), -sec(e + f*x))/(3*f*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_348():
    f = sqrt(d*sec(e + f*x))*(a*sec(e + f*x) + a)**m
    F = -2*sqrt(d*sec(e + f*x))*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S.Half - m, sympy.S(3)/2, sec(e + f*x), -sec(e + f*x))/(f*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_349():
    f = (a*sec(e + f*x) + a)**m/sqrt(d*sec(e + f*x))
    F = 2*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(sympy.S(-1)/2, sympy.S.Half, sympy.S.Half - m, sympy.S.Half, sec(e + f*x), -sec(e + f*x))/(f*sqrt(d*sec(e + f*x))*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_350():
    f = (a*sec(e + f*x) + a)**m/(d*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(sympy.S(-3)/2, sympy.S.Half, sympy.S.Half - m, sympy.S(-1)/2, sec(e + f*x), -sec(e + f*x))/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2)*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_351():
    f = (a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 10*a*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 10*a*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_352():
    f = (a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_353():
    f = (a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_354():
    f = (a*sec(c + d*x) + a)*sqrt(cos(c + d*x))
    F = 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_355():
    f = (a*sec(c + d*x) + a)/sqrt(cos(c + d*x))
    F = 2*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_356():
    f = (a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) - 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_357():
    f = (a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = 6*a*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 2*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_358():
    f = (a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(7)/2)
    F = 6*a*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 10*a*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) - 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 10*a*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_359():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 4*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 32*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 20*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 32*a**2*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 20*a**2*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_360():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 4*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 8*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(7*d) + 12*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 8*a**2*elliptic_f(c/2 + d*x/2, 2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_361():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 4*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 16*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_362():
    f = (a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 4*a**2*elliptic_e(c/2 + d*x/2, 2)/d + 8*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_363():
    f = (a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x))
    F = 2*a**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 4*a**2*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_364():
    f = (a*sec(c + d*x) + a)**2/sqrt(cos(c + d*x))
    F = 4*a**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) - 4*a**2*elliptic_e(c/2 + d*x/2, 2)/d + 8*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_365():
    f = (a*sec(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 16*a**2*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 4*a**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 16*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_366():
    f = (a*sec(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(5)/2)
    F = 12*a**2*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 8*a**2*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**2*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) - 12*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 8*a**2*elliptic_f(c/2 + d*x/2, 2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_367():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 6*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 68*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 44*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 68*a**3*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 44*a**3*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_368():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 6*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 52*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 28*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 52*a**3*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_369():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/d + 36*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**3*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_370():
    f = (a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a**3*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 4*a**3*elliptic_e(c/2 + d*x/2, 2)/d + 20*a**3*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_371():
    f = (a*sec(c + d*x) + a)**3*sqrt(cos(c + d*x))
    F = 6*a**3*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a**3*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) - 4*a**3*elliptic_e(c/2 + d*x/2, 2)/d + 20*a**3*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_372():
    f = (a*sec(c + d*x) + a)**3/sqrt(cos(c + d*x))
    F = 36*a**3*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 2*a**3*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**3*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 36*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**3*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_373():
    f = (a*sec(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(3)/2)
    F = 28*a**3*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 52*a**3*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 6*a**3*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**3*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) - 28*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 52*a**3*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_374():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(d*(a*sec(c + d*x) + a)) + 7*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a*d) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) + 21*elliptic_e(c/2 + d*x/2, 2)/(5*a*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_375():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*sec(c + d*x) + a)) + 5*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) - 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) + 5*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_376():
    f = sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) - elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_377():
    f = 1/((a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - elliptic_e(c/2 + d*x/2, 2)/(a*d) + elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_378():
    f = 1/((a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + elliptic_e(c/2 + d*x/2, 2)/(a*d) + elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_379():
    f = 1/((a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 3*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) - 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) - elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_380():
    f = 1/((a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) - 3*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) + 5*sin(c + d*x)/(3*a*d*cos(c + d*x)**(sympy.S(3)/2)) + 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) + 5*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_381():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) + 56*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a**2*d) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d) + 56*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(a**2*d) - 3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_382():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*sec(c + d*x) + a)**2) + 10*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d) - 7*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 10*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 7*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_383():
    f = sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**2
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x))) + 4*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 5*sin(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_384():
    f = 1/((a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) - elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 2*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) + sin(c + d*x)/(a**2*d*(sec(c + d*x) + 1)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_385():
    f = 1/((a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) + elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_386():
    f = 1/((a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) + elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 2*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - sin(c + d*x)/(a**2*d*(sec(c + d*x) + 1)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_387():
    f = 1/((a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2)) + 4*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) - 4*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 5*sin(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_388():
    f = 1/((a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(9)/2))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(7)/2)) - 7*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) + 10*sin(c + d*x)/(3*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + 7*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 10*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 7*sin(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_389():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**3
    F = -63*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(10*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 4*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a*d*(a*sec(c + d*x) + a)**2) + 77*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(10*a**3*d) - 21*sin(c + d*x)*sqrt(cos(c + d*x))/(2*a**3*d) + 231*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) - 21*elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_390():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**3
    F = -119*sin(c + d*x)*sqrt(cos(c + d*x))/(30*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*sec(c + d*x) + a)**3) - 2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d*(a*sec(c + d*x) + a)**2) + 11*sin(c + d*x)*sqrt(cos(c + d*x))/(2*a**3*d) - 119*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + 11*elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_391():
    f = sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**3
    F = -13*sin(c + d*x)/(6*d*(a**3*sec(c + d*x) + a**3)*sqrt(cos(c + d*x))) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*sqrt(cos(c + d*x))) - 8*sin(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x))) + 49*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) - 13*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_392():
    f = 1/((a*sec(c + d*x) + a)**3*sqrt(cos(c + d*x)))
    F = sin(c + d*x)/(2*d*(a**3*sec(c + d*x) + a**3)*sqrt(cos(c + d*x))) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)) + 2*sin(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x))) - 9*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_393():
    f = 1/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)/(6*d*(a**3*sec(c + d*x) + a**3)*sqrt(cos(c + d*x))) + sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)) - sin(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x))) - elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_394():
    f = 1/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(5)/2))
    F = sin(c + d*x)/(6*d*(a**3*sec(c + d*x) + a**3)*sqrt(cos(c + d*x))) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)) - 4*sin(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2*sqrt(cos(c + d*x))) + elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_395():
    f = 1/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(7)/2))
    F = -9*sin(c + d*x)/(10*d*(a**3*sec(c + d*x) + a**3)*sqrt(cos(c + d*x))) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(5)/2)) - 2*sin(c + d*x)/(5*a*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) + 9*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_396():
    f = 1/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(9)/2))
    F = -13*sin(c + d*x)/(6*d*(a**3*sec(c + d*x) + a**3)*cos(c + d*x)**(sympy.S(3)/2)) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(7)/2)) - 8*sin(c + d*x)/(15*a*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2)) + 49*sin(c + d*x)/(10*a**3*d*sqrt(cos(c + d*x))) - 49*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) - 13*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_397():
    f = 1/((a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(11)/2))
    F = -119*sin(c + d*x)/(30*d*(a**3*sec(c + d*x) + a**3)*cos(c + d*x)**(sympy.S(5)/2)) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(9)/2)) - 2*sin(c + d*x)/(3*a*d*(a*sec(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(7)/2)) - 119*sin(c + d*x)/(10*a**3*d*sqrt(cos(c + d*x))) + 11*sin(c + d*x)/(2*a**3*d*cos(c + d*x)**(sympy.S(3)/2)) + 119*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + 11*elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_398():
    f = sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d*sqrt(a*sec(c + d*x) + a)) + 12*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d*sqrt(a*sec(c + d*x) + a)) + 16*a*sin(c + d*x)*sqrt(cos(c + d*x))/(35*d*sqrt(a*sec(c + d*x) + a)) + 32*a*sin(c + d*x)/(35*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_399():
    f = sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*sqrt(a*sec(c + d*x) + a)) + 8*a*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 16*a*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_400():
    f = sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) + 4*a*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_401():
    f = sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))
    F = 2*a*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_402():
    f = sqrt(a*sec(c + d*x) + a)/sqrt(cos(c + d*x))
    F = 2*sqrt(a)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_403():
    f = sqrt(a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = sqrt(a)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + a*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_404():
    f = sqrt(a*sec(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = 3*sqrt(a)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + 3*a*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a*sin(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_405():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d*sqrt(a*sec(c + d*x) + a)) + 26*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d*sqrt(a*sec(c + d*x) + a)) + 104*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(105*d*sqrt(a*sec(c + d*x) + a)) + 208*a**2*sin(c + d*x)/(105*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_406():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 8*a**2*sin(c + d*x)/(5*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_407():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 8*a**2*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_408():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = 2*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + 2*a**2*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_409():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = 3*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + a**2*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_410():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = 7*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + 7*a**2*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a**2*sin(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_411():
    f = (a*sec(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 11*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + 11*a**2*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 11*a**2*sin(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + a**2*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_412():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(9)/2)
    F = 38*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d*sqrt(a*sec(c + d*x) + a)) + 146*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(105*d*sqrt(a*sec(c + d*x) + a)) + 584*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(315*d*sqrt(a*sec(c + d*x) + a)) + 1168*a**3*sin(c + d*x)/(315*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_413():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 64*a**3*sin(c + d*x)/(21*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 16*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_414():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 64*a**3*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 16*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d) + 2*a*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_415():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + 14*a**3*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_416():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = 5*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/d + a**3*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_417():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = 19*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + 9*a**3*sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_418():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = 25*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + 25*a**3*sin(c + d*x)/(8*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 13*a**3*sin(c + d*x)/(12*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_419():
    f = (a*sec(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 163*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(64*d) + 163*a**3*sin(c + d*x)/(64*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 163*a**3*sin(c + d*x)/(96*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 17*a**3*sin(c + d*x)/(24*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + a**2*sqrt(a*sec(c + d*x) + a)*sin(c + d*x)/(4*d*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_420():
    f = cos(c + d*x)**(sympy.S(5)/2)/sqrt(a*sec(c + d*x) + a)
    F = 2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*sqrt(a*sec(c + d*x) + a)) - 2*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d*sqrt(a*sec(c + d*x) + a)) + 26*sin(c + d*x)/(15*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_421():
    f = cos(c + d*x)**(sympy.S(3)/2)/sqrt(a*sec(c + d*x) + a)
    F = 2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*sqrt(a*sec(c + d*x) + a)) - 2*sin(c + d*x)/(3*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_422():
    f = sqrt(cos(c + d*x))/sqrt(a*sec(c + d*x) + a)
    F = 2*sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_423():
    f = 1/(sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d) - sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_424():
    f = 1/(sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = sin(c + d*x)/(d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d) + sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_425():
    f = 1/(sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(4*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + sin(c + d*x)/(2*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 7*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(4*sqrt(a)*d) - sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_426():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 9*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(10*a*d*sqrt(a*sec(c + d*x) + a)) - 13*sin(c + d*x)*sqrt(cos(c + d*x))/(10*a*d*sqrt(a*sec(c + d*x) + a)) + 49*sin(c + d*x)/(10*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - 15*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_427():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 7*sin(c + d*x)*sqrt(cos(c + d*x))/(6*a*d*sqrt(a*sec(c + d*x) + a)) - 19*sin(c + d*x)/(6*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 11*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_428():
    f = sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + 5*sin(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - 7*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_429():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + 3*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_430():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_431():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + 2*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_432():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(2*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)) + 3*sin(c + d*x)/(2*a*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - 3*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) + 9*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_433():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)) - 17*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)) + 95*sin(c + d*x)*sqrt(cos(c + d*x))/(48*a**2*d*sqrt(a*sec(c + d*x) + a)) - 299*sin(c + d*x)/(48*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) + 163*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_434():
    f = sqrt(cos(c + d*x))/(a*sec(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))) - 13*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + 49*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)*sqrt(cos(c + d*x))) - 75*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_435():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)) - 9*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + 19*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_436():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)) + 5*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + 5*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_437():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)) + 3*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + 3*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_438():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)) - 11*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + 2*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) - 43*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_439():
    f = 1/((a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(9)/2))
    F = -sin(c + d*x)/(4*d*(a*sec(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2)) - 15*sin(c + d*x)/(16*a*d*(a*sec(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)) + 35*sin(c + d*x)/(16*a**2*d*sqrt(a*sec(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - 5*sqrt(cos(c + d*x))*asinh(sqrt(a)*tan(c + d*x)/sqrt(a*sec(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) + 115*sqrt(2)*sqrt(cos(c + d*x))*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*sqrt(a*sec(c + d*x) + a)))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_440():
    f = (d*cos(e + f*x))**n*(a*sec(e + f*x) + a)**3
    F = -a**3*(d*cos(e + f*x))**n*(1 - 4*n)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(1 - n)*(n + 1)*sqrt(sin(e + f*x)**2)) + a**3*(d*cos(e + f*x))**n*(5 - 2*n)*tan(e + f*x)/(f*(1 - n)*(2 - n)) - a**3*(d*cos(e + f*x))**n*(7 - 4*n)*sin(e + f*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(e + f*x)**2)/(f*n*(2 - n)*sqrt(sin(e + f*x)**2)) + (d*cos(e + f*x))**n*(a**3*sec(e + f*x) + a**3)*tan(e + f*x)/(f*(2 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_441():
    f = (d*cos(e + f*x))**n*(a*sec(e + f*x) + a)**2
    F = -a**2*(d*cos(e + f*x))**n*(1 - 2*n)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(1 - n)*(n + 1)*sqrt(sin(e + f*x)**2)) + a**2*(d*cos(e + f*x))**n*tan(e + f*x)/(f*(1 - n)) - 2*a**2*(d*cos(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(e + f*x)**2)/(f*n*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_442():
    f = (d*cos(e + f*x))**n*(a*sec(e + f*x) + a)
    F = -a*(d*cos(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(e + f*x)**2)/(f*n*sqrt(sin(e + f*x)**2)) - a*(d*cos(e + f*x))**(n + 1)*sin(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_443():
    f = (d*cos(e + f*x))**n/(a*sec(e + f*x) + a)
    F = (d*cos(e + f*x))**n*sin(e + f*x)/(f*(a*sec(e + f*x) + a)) + (d*cos(e + f*x))**n*(n + 1)*sin(e + f*x)*cos(e + f*x)**2*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(e + f*x)**2)/(a*f*(n + 2)*sqrt(sin(e + f*x)**2)) - (d*cos(e + f*x))**n*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_444():
    f = (d*cos(e + f*x))**n/(a*sec(e + f*x) + a)**2
    F = -(d*cos(e + f*x))**n*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - (d*cos(e + f*x))**n*(2*n + 3)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(3*a**2*f*sqrt(sin(e + f*x)**2)) + (d*cos(e + f*x))**n*(2*n + 4)*sin(e + f*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(e + f*x)**2)/(3*a**2*f*sqrt(sin(e + f*x)**2)) - (d*cos(e + f*x))**n*(2*n + 4)*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_445():
    f = (a + b*sec(c + d*x))*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*b*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*b*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_446():
    f = (a + b*sec(c + d*x))*sec(c + d*x)**3
    F = a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*atanh(sin(c + d*x))/(2*d) + b*tan(c + d*x)**3/(3*d) + b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_447():
    f = (a + b*sec(c + d*x))*sec(c + d*x)**2
    F = a*tan(c + d*x)/d + b*tan(c + d*x)*sec(c + d*x)/(2*d) + b*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_448():
    f = (a + b*sec(c + d*x))*sec(c + d*x)
    F = a*atanh(sin(c + d*x))/d + b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_449():
    f = a + b*sec(c + d*x)
    F = a*x + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_450():
    f = (a + b*sec(c + d*x))*cos(c + d*x)
    F = a*sin(c + d*x)/d + b*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_451():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**2
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) + b*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_452():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**3
    F = -a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d + b*x/2 + b*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_453():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**4
    F = 3*a*x/8 + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - b*sin(c + d*x)**3/(3*d) + b*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_454():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**5
    F = a*sin(c + d*x)**5/(5*d) - 2*a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d + 3*b*x/8 + b*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*b*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_455():
    f = (a + b*sec(c + d*x))**2*sec(c + d*x)**4
    F = a*b*tan(c + d*x)*sec(c + d*x)**3/(2*d) + 3*a*b*tan(c + d*x)*sec(c + d*x)/(4*d) + 3*a*b*atanh(sin(c + d*x))/(4*d) + b**2*tan(c + d*x)*sec(c + d*x)**4/(5*d) + (5*a**2 + 4*b**2)*tan(c + d*x)**3/(15*d) + (5*a**2 + 4*b**2)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_456():
    f = (a + b*sec(c + d*x))**2*sec(c + d*x)**3
    F = 2*a*b*tan(c + d*x)**3/(3*d) + 2*a*b*tan(c + d*x)/d + b**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + (4*a**2 + 3*b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + (4*a**2 + 3*b**2)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_457():
    f = (a + b*sec(c + d*x))**2*sec(c + d*x)**2
    F = a*b*tan(c + d*x)*sec(c + d*x)/d + a*b*atanh(sin(c + d*x))/d + b**2*tan(c + d*x)*sec(c + d*x)**2/(3*d) + (3*a**2 + 2*b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_458():
    f = (a + b*sec(c + d*x))**2*sec(c + d*x)
    F = 2*a*b*tan(c + d*x)/d + b**2*tan(c + d*x)*sec(c + d*x)/(2*d) + (2*a**2 + b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_459():
    f = (a + b*sec(c + d*x))**2
    F = a**2*x + 2*a*b*atanh(sin(c + d*x))/d + b**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_460():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)
    F = a**2*sin(c + d*x)/d + 2*a*b*x + b**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_461():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**2
    F = a**2*sin(c + d*x)*cos(c + d*x)/(2*d) + 2*a*b*sin(c + d*x)/d + x*(a**2/2 + b**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_462():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**3
    F = -a**2*sin(c + d*x)**3/(3*d) + a*b*x + a*b*sin(c + d*x)*cos(c + d*x)/d + (a**2 + b**2)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_463():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**4
    F = a**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) - 2*a*b*sin(c + d*x)**3/(3*d) + 2*a*b*sin(c + d*x)/d + x*(3*a**2/8 + b**2/2) + (3*a**2 + 4*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_464():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**5
    F = a**2*sin(c + d*x)**5/(5*d) + 3*a*b*x/4 + a*b*sin(c + d*x)*cos(c + d*x)**3/(2*d) + 3*a*b*sin(c + d*x)*cos(c + d*x)/(4*d) + (a**2 + b**2)*sin(c + d*x)/d - (2*a**2 + b**2)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_465():
    f = (a + b*sec(c + d*x))**3*sec(c + d*x)**3
    F = a*(4*a**2 + 9*b**2)*atanh(sin(c + d*x))/(8*d) - a*(6*a**2 - 71*b**2)*tan(c + d*x)*sec(c + d*x)/(120*d) - a*(a + b*sec(c + d*x))**3*tan(c + d*x)/(20*b*d) + (a + b*sec(c + d*x))**4*tan(c + d*x)/(5*b*d) - (a + b*sec(c + d*x))**2*(3*a**2 - 16*b**2)*tan(c + d*x)/(60*b*d) - (3*a**4 - 52*a**2*b**2 - 16*b**4)*tan(c + d*x)/(30*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_466():
    f = (a + b*sec(c + d*x))**3*sec(c + d*x)**2
    F = a*(a + b*sec(c + d*x))**2*tan(c + d*x)/(4*d) + a*(a**2 + 4*b**2)*tan(c + d*x)/(2*d) + b*(2*a**2 + 3*b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*b*(4*a**2 + b**2)*atanh(sin(c + d*x))/(8*d) + (a + b*sec(c + d*x))**3*tan(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_467():
    f = (a + b*sec(c + d*x))**3*sec(c + d*x)
    F = 5*a*b**2*tan(c + d*x)*sec(c + d*x)/(6*d) + a*(2*a**2 + 3*b**2)*atanh(sin(c + d*x))/(2*d) + b*(a + b*sec(c + d*x))**2*tan(c + d*x)/(3*d) + 2*b*(4*a**2 + b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_468():
    f = (a + b*sec(c + d*x))**3
    F = a**3*x + 5*a*b**2*tan(c + d*x)/(2*d) + b**2*(a + b*sec(c + d*x))*tan(c + d*x)/(2*d) + b*(6*a**2 + b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_469():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)
    F = 3*a**2*b*x + 3*a*b**2*atanh(sin(c + d*x))/d + a*(a**2 - b**2)*sin(c + d*x)/d + b**2*(a + b*sec(c + d*x))*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_470():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**2
    F = 5*a**2*b*sin(c + d*x)/(2*d) + a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)/(2*d) + a*x*(a**2 + 6*b**2)/2 + b**3*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_471():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**3
    F = 7*a**2*b*sin(c + d*x)*cos(c + d*x)/(6*d) + a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a*(2*a**2 + 9*b**2)*sin(c + d*x)/(3*d) + b*x*(3*a**2 + 2*b**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_472():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**4
    F = -3*a**2*b*sin(c + d*x)**3/(4*d) + a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*x*(a**2 + 4*b**2)/8 + 3*a*(a**2 + 4*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d) + b*(11*a**2 + 4*b**2)*sin(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_473():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**5
    F = 11*a**2*b*sin(c + d*x)*cos(c + d*x)**3/(20*d) + a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**4/(5*d) - a*(4*a**2 + 15*b**2)*sin(c + d*x)**3/(15*d) + a*(4*a**2 + 15*b**2)*sin(c + d*x)/(5*d) + b*x*(9*a**2 + 4*b**2)/8 + b*(9*a**2 + 4*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_474():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**6
    F = 13*a**2*b*sin(c + d*x)**5/(30*d) + a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**5/(6*d) + a*x*(5*a**2 + 18*b**2)/16 + a*(5*a**2 + 18*b**2)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + a*(5*a**2 + 18*b**2)*sin(c + d*x)*cos(c + d*x)/(16*d) - b*(5*a**2 + b**2)*sin(c + d*x)**3/(3*d) + b*(17*a**2 + 6*b**2)*sin(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_475():
    f = (a + b*sec(c + d*x))**4*sec(c + d*x)**3
    F = -a*(a + b*sec(c + d*x))**4*tan(c + d*x)/(30*b*d) - a*(a + b*sec(c + d*x))**2*(4*a**2 - 53*b**2)*tan(c + d*x)/(120*b*d) - a*(4*a**4 - 121*a**2*b**2 - 128*b**4)*tan(c + d*x)/(60*b*d) - (8*a**4 - 178*a**2*b**2 - 75*b**4)*tan(c + d*x)*sec(c + d*x)/(240*d) + (8*a**4 + 36*a**2*b**2 + 5*b**4)*atanh(sin(c + d*x))/(16*d) + (a + b*sec(c + d*x))**5*tan(c + d*x)/(6*b*d) - (a + b*sec(c + d*x))**3*(4*a**2 - 25*b**2)*tan(c + d*x)/(120*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_476():
    f = (a + b*sec(c + d*x))**4*sec(c + d*x)**2
    F = a*b*(4*a**2 + 3*b**2)*atanh(sin(c + d*x))/(2*d) + a*b*(6*a**2 + 29*b**2)*tan(c + d*x)*sec(c + d*x)/(30*d) + a*(a + b*sec(c + d*x))**3*tan(c + d*x)/(5*d) + (a + b*sec(c + d*x))**4*tan(c + d*x)/(5*d) + (a + b*sec(c + d*x))**2*(3*a**2 + 4*b**2)*tan(c + d*x)/(15*d) + (6*a**4 + 56*a**2*b**2 + 8*b**4)*tan(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_477():
    f = (a + b*sec(c + d*x))**4*sec(c + d*x)
    F = 7*a*b*(a + b*sec(c + d*x))**2*tan(c + d*x)/(12*d) + a*b*(19*a**2 + 16*b**2)*tan(c + d*x)/(6*d) + b**2*(26*a**2 + 9*b**2)*tan(c + d*x)*sec(c + d*x)/(24*d) + b*(a + b*sec(c + d*x))**3*tan(c + d*x)/(4*d) + (8*a**4 + 24*a**2*b**2 + 3*b**4)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_478():
    f = (a + b*sec(c + d*x))**4
    F = a**4*x + 4*a*b**3*tan(c + d*x)*sec(c + d*x)/(3*d) + 2*a*b*(2*a**2 + b**2)*atanh(sin(c + d*x))/d + b**2*(a + b*sec(c + d*x))**2*tan(c + d*x)/(3*d) + b**2*(17*a**2 + 2*b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_479():
    f = (a + b*sec(c + d*x))**4*cos(c + d*x)
    F = 4*a**3*b*x + a**2*(2*a**2 - b**2)*sin(c + d*x)/(2*d) + 3*a*b**3*tan(c + d*x)/d + b**2*(a + b*sec(c + d*x))**2*sin(c + d*x)/(2*d) + b**2*(12*a**2 + b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_480():
    f = (a + b*sec(c + d*x))**4*cos(c + d*x)**2
    F = 3*a**3*b*sin(c + d*x)/d + a**2*x*(a**2 + 12*b**2)/2 + a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)/(2*d) + 4*a*b**3*atanh(sin(c + d*x))/d - b**2*(a**2 - 2*b**2)*tan(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_481():
    f = (a + b*sec(c + d*x))**4*cos(c + d*x)**3
    F = 4*a**3*b*sin(c + d*x)*cos(c + d*x)/(3*d) + a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a**2*(2*a**2 + 17*b**2)*sin(c + d*x)/(3*d) + 2*a*b*x*(a**2 + 2*b**2) + b**4*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_482():
    f = (a + b*sec(c + d*x))**4*cos(c + d*x)**4
    F = 5*a**3*b*sin(c + d*x)*cos(c + d*x)**2/(6*d) + a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + a**2*(3*a**2 + 22*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d) + 4*a*b*(2*a**2 + 3*b**2)*sin(c + d*x)/(3*d) + x*(3*a**4/8 + 3*a**2*b**2 + b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_483():
    f = (a + b*sec(c + d*x))**4*cos(c + d*x)**5
    F = 3*a**3*b*sin(c + d*x)*cos(c + d*x)**3/(5*d) + a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**4/(5*d) - a**2*(4*a**2 + 27*b**2)*sin(c + d*x)**3/(15*d) + a*b*x*(3*a**2 + 4*b**2)/2 + a*b*(3*a**2 + 4*b**2)*sin(c + d*x)*cos(c + d*x)/(2*d) + (4*a**4 + 29*a**2*b**2 + 5*b**4)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_484():
    f = (a + b*sec(c + d*x))**4*cos(c + d*x)**6
    F = 7*a**3*b*sin(c + d*x)*cos(c + d*x)**4/(15*d) + a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) + a**2*(5*a**2 + 32*b**2)*sin(c + d*x)*cos(c + d*x)**3/(24*d) - 4*a*b*(4*a**2 + 5*b**2)*sin(c + d*x)**3/(15*d) + 4*a*b*(4*a**2 + 5*b**2)*sin(c + d*x)/(5*d) + x*(5*a**4/16 + 9*a**2*b**2/4 + b**4/2) + (5*a**4 + 36*a**2*b**2 + 8*b**4)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_485():
    f = (a + b*sec(c + d*x))**5
    F = a**5*x + 11*a*b**2*(a + b*sec(c + d*x))**2*tan(c + d*x)/(12*d) + a*b**2*(53*a**2 + 20*b**2)*tan(c + d*x)/(6*d) + b**3*(58*a**2 + 9*b**2)*tan(c + d*x)*sec(c + d*x)/(24*d) + b**2*(a + b*sec(c + d*x))**3*tan(c + d*x)/(4*d) + b*(40*a**4 + 40*a**2*b**2 + 3*b**4)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_486():
    f = sec(c + d*x)**5/(a + b*sec(c + d*x))
    F = 2*a**4*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*sqrt(a - b)*sqrt(a + b)) - a*tan(c + d*x)*sec(c + d*x)/(2*b**2*d) - a*(2*a**2 + b**2)*atanh(sin(c + d*x))/(2*b**4*d) + tan(c + d*x)*sec(c + d*x)**2/(3*b*d) + (3*a**2 + 2*b**2)*tan(c + d*x)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_487():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))
    F = -2*a**3*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*sqrt(a - b)*sqrt(a + b)) - a*tan(c + d*x)/(b**2*d) + tan(c + d*x)*sec(c + d*x)/(2*b*d) + (2*a**2 + b**2)*atanh(sin(c + d*x))/(2*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_488():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))
    F = 2*a**2*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*sqrt(a - b)*sqrt(a + b)) - a*atanh(sin(c + d*x))/(b**2*d) + tan(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_489():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))
    F = -2*a*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b*d*sqrt(a - b)*sqrt(a + b)) + atanh(sin(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_490():
    f = sec(c + d*x)/(a + b*sec(c + d*x))
    F = 2*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_491():
    f = 1/(a + b*sec(c + d*x))
    F = -2*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*sqrt(a - b)*sqrt(a + b)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_492():
    f = cos(c + d*x)/(a + b*sec(c + d*x))
    F = sin(c + d*x)/(a*d) + 2*b**2*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*sqrt(a - b)*sqrt(a + b)) - b*x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_493():
    f = cos(c + d*x)**2/(a + b*sec(c + d*x))
    F = sin(c + d*x)*cos(c + d*x)/(2*a*d) - b*sin(c + d*x)/(a**2*d) - 2*b**3*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*sqrt(a - b)*sqrt(a + b)) + x*(a**2 + 2*b**2)/(2*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_494():
    f = cos(c + d*x)**3/(a + b*sec(c + d*x))
    F = sin(c + d*x)*cos(c + d*x)**2/(3*a*d) - b*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) + (2*a**2 + 3*b**2)*sin(c + d*x)/(3*a**3*d) + 2*b**4*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*sqrt(a - b)*sqrt(a + b)) - b*x*(a**2 + 2*b**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_495():
    f = cos(c + d*x)**4/(a + b*sec(c + d*x))
    F = sin(c + d*x)*cos(c + d*x)**3/(4*a*d) - b*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d) + (3*a**2 + 4*b**2)*sin(c + d*x)*cos(c + d*x)/(8*a**3*d) - b*(2*a**2 + 3*b**2)*sin(c + d*x)/(3*a**4*d) - 2*b**5*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*sqrt(a - b)*sqrt(a + b)) + x*(3*a**4 + 4*a**2*b**2 + 8*b**4)/(8*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_496():
    f = sec(c + d*x)**5/(a + b*sec(c + d*x))**2
    F = -2*a**3*(3*a**2 - 4*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - a**2*tan(c + d*x)*sec(c + d*x)**2/(b*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - a*(3*a**2 - 2*b**2)*tan(c + d*x)/(b**3*d*(a**2 - b**2)) + (3*a**2 - b**2)*tan(c + d*x)*sec(c + d*x)/(2*b**2*d*(a**2 - b**2)) + (6*a**2 + b**2)*atanh(sin(c + d*x))/(2*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_497():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))**2
    F = -a**2*tan(c + d*x)*sec(c + d*x)/(b*d*(a + b*sec(c + d*x))*(a**2 - b**2)) + 2*a**2*(2*a**2 - 3*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - 2*a*atanh(sin(c + d*x))/(b**3*d) + (2*a**2 - b**2)*tan(c + d*x)/(b**2*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_498():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = -a**2*tan(c + d*x)/(b*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - 2*a*(a**2 - 2*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + atanh(sin(c + d*x))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_499():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = a*tan(c + d*x)/(d*(a + b*sec(c + d*x))*(a**2 - b**2)) - 2*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_500():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**2
    F = 2*a*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - b*tan(c + d*x)/(d*(a + b*sec(c + d*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_501():
    f = (a + b*sec(c + d*x))**(-2)
    F = b**2*tan(c + d*x)/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) - 2*b*(2*a**2 - b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x/a**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_502():
    f = cos(c + d*x)/(a + b*sec(c + d*x))**2
    F = b**2*sin(c + d*x)/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) + (a**2 - 2*b**2)*sin(c + d*x)/(a**2*d*(a**2 - b**2)) + 2*b**2*(3*a**2 - 2*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - 2*b*x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_503():
    f = cos(c + d*x)**2/(a + b*sec(c + d*x))**2
    F = b**2*sin(c + d*x)*cos(c + d*x)/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) + (a**2 - 3*b**2)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d*(a**2 - b**2)) - b*(2*a**2 - 3*b**2)*sin(c + d*x)/(a**3*d*(a**2 - b**2)) - 2*b**3*(4*a**2 - 3*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x*(a**2 + 6*b**2)/(2*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_504():
    f = cos(c + d*x)**3/(a + b*sec(c + d*x))**2
    F = b**2*sin(c + d*x)*cos(c + d*x)**2/(a*d*(a + b*sec(c + d*x))*(a**2 - b**2)) + (a**2 - 4*b**2)*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d*(a**2 - b**2)) - b*(a**2 - 2*b**2)*sin(c + d*x)*cos(c + d*x)/(a**3*d*(a**2 - b**2)) + (2*a**4 + 7*a**2*b**2 - 12*b**4)*sin(c + d*x)/(3*a**4*d*(a**2 - b**2)) + 2*b**4*(5*a**2 - 4*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - b*x*(a**2 + 4*b**2)/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_505():
    f = sec(c + d*x)**5/(a + b*sec(c + d*x))**3
    F = 3*a**3*(a**2 - 2*b**2)*tan(c + d*x)/(2*b**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - a**2*tan(c + d*x)*sec(c + d*x)**2/(2*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + 3*a**2*(2*a**4 - 5*a**2*b**2 + 4*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - 3*a*atanh(sin(c + d*x))/(b**4*d) + (3*a**2 - 2*b**2)*tan(c + d*x)/(2*b**3*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_506():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))**3
    F = -a**2*tan(c + d*x)*sec(c + d*x)/(2*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) - a**2*(2*a**2 - 5*b**2)*tan(c + d*x)/(2*b**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - a*(2*a**4 - 5*a**2*b**2 + 6*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + atanh(sin(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_507():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))**3
    F = -a**2*tan(c + d*x)/(2*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + a*(a**2 - 4*b**2)*tan(c + d*x)/(2*b*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) + (a**2 + 2*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_508():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))**3
    F = -3*a*b*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a*tan(c + d*x)/(d*(a + b*sec(c + d*x))**2*(2*a**2 - 2*b**2)) + (a**2 + 2*b**2)*tan(c + d*x)/(2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_509():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**3
    F = -3*a*b*tan(c + d*x)/(2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - b*tan(c + d*x)/(d*(a + b*sec(c + d*x))**2*(2*a**2 - 2*b**2)) + (2*a**2 + b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_510():
    f = (a + b*sec(c + d*x))**(-3)
    F = b**2*tan(c + d*x)/(2*a*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + b**2*(5*a**2 - 2*b**2)*tan(c + d*x)/(2*a**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) - b*(6*a**4 - 5*a**2*b**2 + 2*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_511():
    f = cos(c + d*x)/(a + b*sec(c + d*x))**3
    F = b**2*sin(c + d*x)/(2*a*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + 3*b**2*(2*a**2 - b**2)*sin(c + d*x)/(2*a**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) + (2*a**4 - 11*a**2*b**2 + 6*b**4)*sin(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) + 3*b**2*(4*a**4 - 5*a**2*b**2 + 2*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - 3*b*x/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_512():
    f = cos(c + d*x)**2/(a + b*sec(c + d*x))**3
    F = b**2*sin(c + d*x)*cos(c + d*x)/(2*a*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)) + b**2*(7*a**2 - 4*b**2)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**2) + (a**4 - 10*a**2*b**2 + 6*b**4)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) - 3*b*(2*a**4 - 7*a**2*b**2 + 4*b**4)*sin(c + d*x)/(2*a**4*d*(a**2 - b**2)**2) - b**3*(20*a**4 - 29*a**2*b**2 + 12*b**4)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + x*(a**2 + 12*b**2)/(2*a**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_513():
    f = sec(c + d*x)**6/(a + b*sec(c + d*x))**4
    F = a**3*(4*a**4 - 11*a**2*b**2 + 12*b**4)*tan(c + d*x)/(2*b**4*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - a**2*tan(c + d*x)*sec(c + d*x)**3/(3*b*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) - a**2*(4*a**2 - 9*b**2)*tan(c + d*x)*sec(c + d*x)**2/(6*b**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + a**2*(8*a**6 - 28*a**4*b**2 + 35*a**2*b**4 - 20*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - 4*a*atanh(sin(c + d*x))/(b**5*d) + (12*a**4 - 23*a**2*b**2 + 6*b**4)*tan(c + d*x)/(6*b**4*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_514():
    f = sec(c + d*x)**5/(a + b*sec(c + d*x))**4
    F = a**3*(3*a**2 - 8*b**2)*tan(c + d*x)/(6*b**3*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) - a**2*tan(c + d*x)*sec(c + d*x)**2/(3*b*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) - a**2*(9*a**4 - 28*a**2*b**2 + 34*b**4)*tan(c + d*x)/(6*b**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - a*(2*a**6 - 7*a**4*b**2 + 8*a**2*b**4 - 8*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + atanh(sin(c + d*x))/(b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_515():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))**4
    F = -a**2*tan(c + d*x)*sec(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) - a**2*(2*a**2 - 7*b**2)*tan(c + d*x)/(6*b**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + a*(2*a**4 - 5*a**2*b**2 + 18*b**4)*tan(c + d*x)/(6*b**2*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - b*(3*a**2 + 2*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_516():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))**4
    F = -a**2*tan(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + a*(a**2 + 4*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + a*(a**2 - 6*b**2)*tan(c + d*x)/(6*b*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + (a**4 - 10*a**2*b**2 - 6*b**4)*tan(c + d*x)/(6*b*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_517():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))**4
    F = a*(2*a**2 + 13*b**2)*tan(c + d*x)/(6*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) + a*tan(c + d*x)/(d*(a + b*sec(c + d*x))**3*(3*a**2 - 3*b**2)) - b*(4*a**2 + b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (2*a**2 + 3*b**2)*tan(c + d*x)/(6*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_518():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**4
    F = -5*a*b*tan(c + d*x)/(6*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + a*(2*a**2 + 3*b**2)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - b*(11*a**2 + 4*b**2)*tan(c + d*x)/(6*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - b*tan(c + d*x)/(d*(a + b*sec(c + d*x))**3*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_519():
    f = (a + b*sec(c + d*x))**(-4)
    F = b**2*tan(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + b**2*(8*a**2 - 3*b**2)*tan(c + d*x)/(6*a**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + b**2*(26*a**4 - 17*a**2*b**2 + 6*b**4)*tan(c + d*x)/(6*a**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) - b*(8*a**6 - 8*a**4*b**2 + 7*a**2*b**4 - 2*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + x/a**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_520():
    f = cos(c + d*x)/(a + b*sec(c + d*x))**4
    F = b**2*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + b**2*(9*a**2 - 4*b**2)*sin(c + d*x)/(6*a**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + b**2*(12*a**4 - 11*a**2*b**2 + 4*b**4)*sin(c + d*x)/(2*a**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) + (6*a**6 - 65*a**4*b**2 + 68*a**2*b**4 - 24*b**6)*sin(c + d*x)/(6*a**4*d*(a**2 - b**2)**3) + b**2*(20*a**6 - 35*a**4*b**2 + 28*a**2*b**4 - 8*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - 4*b*x/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_521():
    f = cos(c + d*x)**2/(a + b*sec(c + d*x))**4
    F = b**2*sin(c + d*x)*cos(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**3*(a**2 - b**2)) + 5*b**2*(2*a**2 - b**2)*sin(c + d*x)*cos(c + d*x)/(6*a**2*d*(a + b*sec(c + d*x))**2*(a**2 - b**2)**2) + b**2*(48*a**4 - 53*a**2*b**2 + 20*b**4)*sin(c + d*x)*cos(c + d*x)/(6*a**3*d*(a + b*sec(c + d*x))*(a**2 - b**2)**3) + (a**6 - 23*a**4*b**2 + 27*a**2*b**4 - 10*b**6)*sin(c + d*x)*cos(c + d*x)/(2*a**4*d*(a**2 - b**2)**3) - b*(24*a**6 - 146*a**4*b**2 + 167*a**2*b**4 - 60*b**6)*sin(c + d*x)/(6*a**5*d*(a**2 - b**2)**3) - b**3*(40*a**6 - 84*a**4*b**2 + 69*a**2*b**4 - 20*b**6)*atanh(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**6*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + x*(a**2 + 20*b**2)/(2*a**6)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_522():
    f = 1/(5*sec(c + d*x) + 3)
    F = -x/12 + 5*atan(sin(c + d*x)/(cos(c + d*x) + 3))/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_523():
    f = (5*sec(c + d*x) + 3)**(-2)
    F = 29*x/576 + 35*atan(sin(c + d*x)/(cos(c + d*x) + 3))/(288*d) - 25*tan(c + d*x)/(48*d*(5*sec(c + d*x) + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_524():
    f = (5*sec(c + d*x) + 3)**(-3)
    F = -1007*x/55296 + 3055*atan(sin(c + d*x)/(cos(c + d*x) + 3))/(27648*d) - 125*tan(c + d*x)/(4608*d*(5*sec(c + d*x) + 3)) - 25*tan(c + d*x)/(96*d*(5*sec(c + d*x) + 3)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_525():
    f = (5*sec(c + d*x) + 3)**(-4)
    F = 21553*x/2654208 + 11215*atan(sin(c + d*x)/(cos(c + d*x) + 3))/(1327104*d) - 16925*tan(c + d*x)/(221184*d*(5*sec(c + d*x) + 3)) - 25*tan(c + d*x)/(4608*d*(5*sec(c + d*x) + 3)**2) - 25*tan(c + d*x)/(144*d*(5*sec(c + d*x) + 3)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_526():
    f = 1/(3*sec(c + d*x) + 5)
    F = x/5 + 3*log(-sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(20*d) - 3*log(sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_527():
    f = (3*sec(c + d*x) + 5)**(-2)
    F = x/25 + 123*log(-sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(1600*d) - 123*log(sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(1600*d) + 9*tan(c + d*x)/(80*d*(3*sec(c + d*x) + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_528():
    f = (3*sec(c + d*x) + 5)**(-3)
    F = x/125 + 8361*log(-sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(256000*d) - 8361*log(sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(256000*d) + 963*tan(c + d*x)/(12800*d*(3*sec(c + d*x) + 5)) + 9*tan(c + d*x)/(160*d*(3*sec(c + d*x) + 5)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_529():
    f = (3*sec(c + d*x) + 5)**(-4)
    F = x/625 + 278151*log(-sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(20480000*d) - 278151*log(sin(c/2 + d*x/2) + 2*cos(c/2 + d*x/2))/(20480000*d) + 38733*tan(c + d*x)/(1024000*d*(3*sec(c + d*x) + 5)) + 519*tan(c + d*x)/(12800*d*(3*sec(c + d*x) + 5)**2) + 3*tan(c + d*x)/(80*d*(3*sec(c + d*x) + 5)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_530():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**3
    F = -4*a*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(15*b*d) + 2*(a + b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(2*a + 9*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**2*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(2*a**2 - 9*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_531():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**2
    F = -2*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**2*d) + 2*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(3*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_532():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)
    F = sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_533():
    f = sqrt(a + b*sec(c + d*x))
    F = (Integer(-2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Symbol('a') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + Symbol('b'))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_534():
    f = sqrt(a + b*sec(c + d*x))*cos(c + d*x)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_535():
    f = sqrt(a + b*sec(c + d*x))*cos(c + d*x)**2
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a')) + Symbol('b')) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_536():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**4
    F = -8*a*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(63*b**2*d) + 2*a*sqrt(a + b*sec(c + d*x))*(8*a**2 + 39*b**2)*tan(c + d*x)/(315*b**2*d) + 2*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)*sec(c + d*x)/(9*b*d) + (a + b*sec(c + d*x))**(sympy.S(3)/2)*(16*a**2 + 98*b**2)*tan(c + d*x)/(315*b**2*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(8*a**3 + 6*a**2*b + 39*a*b**2 - 147*b**3)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(8*a**4 + 33*a**2*b**2 + 147*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_537():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**3
    F = -4*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(35*b*d) + 4*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*(3*a**2 - 41*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**3*d) + 2*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*b*d) - sqrt(a + b*sec(c + d*x))*(12*a**2 - 50*b**2)*tan(c + d*x)/(105*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(6*a**2 + 57*a*b - 25*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_538():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2
    F = 2*a*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(5*d) + 2*(a + b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*(2*a - 6*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(5*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(a**2 + 3*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(5*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_539():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)
    F = -8*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b*d) + 2*b*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3*a - b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_540():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(-2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + ((Integer(2) * ((Integer(2) * Symbol('a')) + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_541():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)
    F = ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_542():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2
    F = ((Integer(5) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a')) + (Integer(5) * Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((Integer(5) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_543():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**4
    F = -8*a*(a + b*sec(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)/(99*b**2*d) + 2*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(8*a**2 + 67*b**2)*tan(c + d*x)/(693*b**2*d) - 2*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*(8*a**4 + 51*a**2*b**2 + 741*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(693*b**4*d) + 2*(a + b*sec(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)*sec(c + d*x)/(11*b*d) + (a + b*sec(c + d*x))**(sympy.S(5)/2)*(16*a**2 + 162*b**2)*tan(c + d*x)/(693*b**2*d) + sqrt(a + b*sec(c + d*x))*(16*a**4 + 114*a**2*b**2 + 270*b**4)*tan(c + d*x)/(693*b**2*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(8*a**4 + 6*a**3*b + 57*a**2*b**2 - 606*a*b**3 + 135*b**4)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(693*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_544():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**3
    F = -4*a*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(63*b*d) - 4*a*sqrt(a + b*sec(c + d*x))*(5*a**2 - 57*b**2)*tan(c + d*x)/(315*b*d) + 2*(a + b*sec(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)/(9*b*d) - (a + b*sec(c + d*x))**(sympy.S(3)/2)*(20*a**2 - 98*b**2)*tan(c + d*x)/(315*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(10*a**3 + 165*a**2*b - 114*a*b**2 + 147*b**3)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**2*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(10*a**4 - 279*a**2*b**2 - 147*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(315*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_545():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 2*a*(a + b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(7*d) - 2*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*(3*a**2 + 29*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(21*b**2*d) + 2*(a + b*sec(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(7*d) + sqrt(a + b*sec(c + d*x))*(6*a**2 + 10*b**2)*tan(c + d*x)/(21*d) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3*a**2 - 24*a*b + 5*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(21*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_546():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)
    F = 16*a*b*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(15*d) + 2*b*(a + b*sec(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(5*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(23*a**2 + 9*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(15*a**2 - 8*a*b + 9*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_547():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(-14) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('b'))) + (Symbol('b'))**(Integer(2))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_548():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(6) * Symbol('a') * Symbol('b')) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_549():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2
    F = ((Integer(9) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(9) * Symbol('a') * Symbol('b')) + (Integer(8) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(15) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + ((Integer(9) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_550():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**3
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2))) + (Integer(33) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(16) * (Symbol('a'))**(Integer(2))) + (Integer(26) * Symbol('a') * Symbol('b')) + (Integer(33) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(8) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((((Integer(16) * (Symbol('a'))**(Integer(2))) + (Integer(33) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_551():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**4
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(284) * (Symbol('a'))**(Integer(2))) + (Integer(15) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(192) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(72) * (Symbol('a'))**(Integer(3))) + (Integer(284) * (Symbol('a'))**(Integer(2)) * Symbol('b')) + (Integer(118) * Symbol('a') * (Symbol('b'))**(Integer(2))) + (Integer(15) * (Symbol('b'))**(Integer(3)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(192) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(48) * (Symbol('a'))**(Integer(4))) + (Integer(120) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(4))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(64) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Integer(284) * (Symbol('a'))**(Integer(2))) + (Integer(15) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((((Integer(36) * (Symbol('a'))**(Integer(2))) + (Integer(59) * (Symbol('b'))**(Integer(2)))) * sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(96) * Symbol('d')))**(Integer(-1))) + ((Integer(17) * Symbol('a') * Symbol('b') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**(Integer(3)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_552():
    f = (a + b*sec(c + d*x))**(sympy.S(7)/2)
    F = ((Integer(-2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(58) * (Symbol('a'))**(Integer(2))) + (Integer(9) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(60) * (Symbol('a'))**(Integer(3))) + (Integer(-1) * (Integer(58) * (Symbol('a'))**(Integer(2)) * Symbol('b'))) + (Integer(22) * Symbol('a') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(9) * (Symbol('b'))**(Integer(3))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(26) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_553():
    f = sec(c + d*x)**5/sqrt(a + b*sec(c + d*x))
    F = -12*a*sqrt(a + b*sec(c + d*x))*tan(c + d*x)*sec(c + d*x)/(35*b**2*d) + 8*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*(12*a**2 + 11*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**5*d) + 2*sqrt(a + b*sec(c + d*x))*tan(c + d*x)*sec(c + d*x)**2/(7*b*d) + sqrt(a + b*sec(c + d*x))*(48*a**2 + 50*b**2)*tan(c + d*x)/(105*b**3*d) + 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(48*a**3 - 12*a**2*b + 44*a*b**2 + 25*b**3)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(105*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_554():
    f = sec(c + d*x)**4/sqrt(a + b*sec(c + d*x))
    F = -8*a*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(15*b**2*d) + 2*sqrt(a + b*sec(c + d*x))*tan(c + d*x)*sec(c + d*x)/(5*b*d) - 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(8*a**2 - 2*a*b + 9*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**3*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*(8*a**2 + 9*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(15*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_555():
    f = sec(c + d*x)**3/sqrt(a + b*sec(c + d*x))
    F = 4*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**3*d) + 2*sqrt(a + b*sec(c + d*x))*tan(c + d*x)/(3*b*d) + 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_556():
    f = sec(c + d*x)**2/sqrt(a + b*sec(c + d*x))
    F = -2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-2*a + 2*b)*sqrt(a + b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_557():
    f = sec(c + d*x)/sqrt(a + b*sec(c + d*x))
    F = 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_558():
    f = 1/sqrt(a + b*sec(c + d*x))
    F = (Integer(-2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_559():
    f = cos(c + d*x)/sqrt(a + b*sec(c + d*x))
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_560():
    f = cos(c + d*x)**2/sqrt(a + b*sec(c + d*x))
    F = ((Integer(-3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * Symbol('a')) + (Integer(-1) * (Integer(3) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_561():
    f = sec(c + d*x)**5/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*tan(c + d*x)*sec(c + d*x)**2/(b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - 2*a*sqrt(a + b*sec(c + d*x))*(8*a**2 - 3*b**2)*tan(c + d*x)/(5*b**3*d*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(12*a**2 - 2*b**2)*tan(c + d*x)*sec(c + d*x)/(5*b**2*d*(a**2 - b**2)) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(8*a + 6*b)*(4*a**2 + b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(5*b**4*d*sqrt(a + b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-32*a**4 + 16*a**2*b**2 + 6*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(5*b**5*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_562():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*tan(c + d*x)*sec(c + d*x)/(b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + 2*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(8*a**2 - 5*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**4*d*sqrt(a + b)) + sqrt(a + b*sec(c + d*x))*(8*a**2 - 2*b**2)*tan(c + d*x)/(3*b**2*d*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(4*a + b)*(4*a + 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(3*b**3*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_563():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*tan(c + d*x)/(b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(4*a + 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d*sqrt(a + b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-4*a**2 + 2*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**3*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_564():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*a*tan(c + d*x)/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + 2*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d*sqrt(a + b)) + 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_565():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*b*tan(c + d*x)/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*sqrt(a + b)) + 2*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_566():
    f = (a + b*sec(c + d*x))**(sympy.S(-3)/2)
    F = ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_567():
    f = cos(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (((Symbol('a') + (Integer(3) * Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (sympy.sin((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_568():
    f = cos(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(-1) * (((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('b'))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(15) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_569():
    f = sec(c + d*x)**5/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 4*a**3*(3*a**2 - 5*b**2)*tan(c + d*x)/(3*b**3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) - 2*a**2*tan(c + d*x)*sec(c + d*x)**2/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 8*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(4*a**4 - 7*a**2*b**2 + 2*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**5*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b)) + sqrt(a + b*sec(c + d*x))*(4*a**2 - 2*b**2)*tan(c + d*x)/(3*b**3*d*(a**2 - b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(32*a**4 + 24*a**3*b - 32*a**2*b**2 - 18*a*b**3 - 2*b**4)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**4*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_570():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*tan(c + d*x)*sec(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - 8*a**2*(a**2 - 2*b**2)*tan(c + d*x)/(3*b**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) - sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(16*a**3 + 12*a**2*b - 18*a*b**2 - 6*b**3)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**3*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(-16*a**4 + 30*a**2*b**2 - 6*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**4*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_571():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*tan(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 4*a*(a**2 - 3*b**2)*tan(c + d*x)/(3*b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + 4*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(a**2 - 3*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**3*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(4*a**2 + 6*a*b - 6*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_572():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*a*tan(c + d*x)/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + (2*a**2 + 6*b**2)*tan(c + d*x)/(3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(2*a - 6*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(2*a**2 + 6*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b**2*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_573():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -8*a*b*tan(c + d*x)/(3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) - 8*a*sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b)) - 2*b*tan(c + d*x)/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + sqrt(b*(1 - sec(c + d*x))/(a + b))*sqrt(-b*(sec(c + d*x) + 1)/(a - b))*(6*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*sec(c + d*x))/sqrt(a + b)), (a + b)/(a - b))/(b*d*(a + b)**(sympy.S(3)/2)*(3*a - 3*b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_574():
    f = (a + b*sec(c + d*x))**(sympy.S(-5)/2)
    F = ((Integer(2) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(6) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * Symbol('b'))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_575():
    f = cos(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(3))) + (Integer(21) * (Symbol('a'))**(Integer(2)) * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(3))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(5) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + (sympy.sin((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_576():
    f = cos(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(-1) * (((Integer(33) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(170) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(105) * (Symbol('b'))**(Integer(4)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (((Symbol('a') + (Integer(3) * Symbol('b'))) * ((Integer(6) * (Symbol('a'))**(Integer(3))) + (Integer(-1) * (Integer(45) * (Symbol('a'))**(Integer(2)) * Symbol('b'))) + (Integer(35) * (Symbol('b'))**(Integer(3)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(35) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * (Symbol('a'))**(Integer(5)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.cos((Symbol('c') + (Symbol('d') * x))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(27) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(35) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(33) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(170) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(105) * (Symbol('b'))**(Integer(4)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_577():
    f = (a + b*sec(c + d*x))**(sympy.S(-7)/2)
    F = ((Integer(2) * ((Integer(58) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(41) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(45) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(13) * (Symbol('a'))**(Integer(3)) * Symbol('b'))) + (Integer(-1) * (Integer(36) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(13) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(58) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(41) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_578():
    f = (a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 6*b*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 6*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_579():
    f = (a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_580():
    f = (a + b*sec(c + d*x))*sqrt(sec(c + d*x))
    F = 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_581():
    f = (a + b*sec(c + d*x))/sqrt(sec(c + d*x))
    F = 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_582():
    f = (a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_583():
    f = (a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_584():
    f = (a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(7)/2)
    F = 10*a*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 10*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*b*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_585():
    f = (a + b*sec(c + d*x))**2*sec(c + d*x)**(sympy.S(5)/2)
    F = 4*a*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 12*a*b*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 12*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + (14*a**2 + 10*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + (14*a**2 + 10*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_586():
    f = (a + b*sec(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 4*a*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 4*a*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*b**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + (10*a**2 + 6*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (10*a**2 + 6*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_587():
    f = (a + b*sec(c + d*x))**2*sqrt(sec(c + d*x))
    F = 4*a*b*sin(c + d*x)*sqrt(sec(c + d*x))/d - 4*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + (6*a**2 + 2*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_588():
    f = (a + b*sec(c + d*x))**2/sqrt(sec(c + d*x))
    F = 4*a*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b**2*sin(c + d*x)*sqrt(sec(c + d*x))/d + (2*a**2 - 2*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_589():
    f = (a + b*sec(c + d*x))**2/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (2*a**2 + 6*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_590():
    f = (a + b*sec(c + d*x))**2/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*a*b*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (6*a**2 + 10*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_591():
    f = (a + b*sec(c + d*x))**2/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 4*a*b*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 12*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (10*a**2 + 14*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (10*a**2 + 14*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_592():
    f = (a + b*sec(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2)
    F = 32*a*b**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + 2*a*(5*a**2 + 9*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 2*a*(5*a**2 + 9*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*(a + b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*b*(21*a**2 + 5*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 2*b*(21*a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_593():
    f = (a + b*sec(c + d*x))**3*sqrt(sec(c + d*x))
    F = 8*a*b**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a*(a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b**2*(a + b*sec(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*b*(5*a**2 + b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 6*b*(5*a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_594():
    f = (a + b*sec(c + d*x))**3/sqrt(sec(c + d*x))
    F = 16*a*b**2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*a*(a**2 - 3*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b**2*(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*b*(9*a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_595():
    f = (a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*(a**2 + 9*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - 2*b*(a**2 - 3*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*b*(3*a**2 - b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_596():
    f = (a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(5)/2)
    F = 8*a**2*b*sin(c + d*x)/(5*d*sqrt(sec(c + d*x))) + 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a*(a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b*(a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_597():
    f = (a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(7)/2)
    F = 32*a**2*b*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a*(5*a**2 + 21*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*(5*a**2 + 21*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*b*(9*a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_598():
    f = (a + b*sec(c + d*x))**3/sec(c + d*x)**(sympy.S(9)/2)
    F = 40*a**2*b*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a*(7*a**2 + 27*b**2)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*(7*a**2 + 27*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 2*b*(15*a**2 + 7*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*b*(15*a**2 + 7*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_599():
    f = (a + b*sec(c + d*x))**4*sec(c + d*x)**(sympy.S(3)/2)
    F = 44*a*b**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d) + 8*a*b*(7*a**2 + 5*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 8*a*b*(7*a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*b**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(9*d) + 14*b**2*(7*a**2 + b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(45*d) + (30*a**4 + 108*a**2*b**2 + 14*b**4)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) - (30*a**4 + 108*a**2*b**2 + 14*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_600():
    f = (a + b*sec(c + d*x))**4*sqrt(sec(c + d*x))
    F = 36*a*b**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + 8*a*b*(5*a**2 + 3*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 8*a*b*(5*a**2 + 3*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*b**2*(39*a**2 + 5*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + (42*a**4 + 84*a**2*b**2 + 10*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_601():
    f = (a + b*sec(c + d*x))**4/sqrt(sec(c + d*x))
    F = 28*a*b**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d) + 8*a*b*(3*a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*b**2*(a + b*sec(c + d*x))**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*(29*a**2 + 3*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + (10*a**4 - 60*a**2*b**2 - 6*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_602():
    f = (a + b*sec(c + d*x))**4/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) - 4*a*b*(a**2 - 6*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 8*a*b*(a**2 - b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d - 2*b**2*(a**2 - b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + (2*a**4 + 36*a**2*b**2 + 2*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_603():
    f = (a + b*sec(c + d*x))**4/sec(c + d*x)**(sympy.S(5)/2)
    F = 28*a**3*b*sin(c + d*x)/(15*d*sqrt(sec(c + d*x))) + 2*a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 8*a*b*(a**2 + 3*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - 2*b**2*(a**2 - 5*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + (6*a**4 + 60*a**2*b**2 - 10*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_604():
    f = (a + b*sec(c + d*x))**4/sec(c + d*x)**(sympy.S(7)/2)
    F = 36*a**3*b*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(5*a**2 + 39*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 8*a*b*(3*a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (10*a**4 + 84*a**2*b**2 + 42*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_605():
    f = (a + b*sec(c + d*x))**4/sec(c + d*x)**(sympy.S(9)/2)
    F = 44*a**3*b*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 14*a**2*(a**2 + 7*b**2)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 8*a*b*(5*a**2 + 7*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 8*a*b*(5*a**2 + 7*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (14*a**4 + 108*a**2*b**2 + 30*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_606():
    f = (a + b*sec(c + d*x))**4/sec(c + d*x)**(sympy.S(11)/2)
    F = 52*a**3*b*sin(c + d*x)/(99*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*a**2*(a + b*sec(c + d*x))**2*sin(c + d*x)/(11*d*sec(c + d*x)**(sympy.S(9)/2)) + 2*a**2*(9*a**2 + 59*b**2)*sin(c + d*x)/(77*d*sec(c + d*x)**(sympy.S(5)/2)) + 8*a*b*(7*a**2 + 9*b**2)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 8*a*b*(7*a**2 + 9*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + (90*a**4 + 660*a**2*b**2 + 154*b**4)*sin(c + d*x)/(231*d*sqrt(sec(c + d*x))) + (90*a**4 + 660*a**2*b**2 + 154*b**4)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(231*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_607():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))
    F = ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_608():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))
    F = ((Integer(-2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_609():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))
    F = (Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_610():
    f = sqrt(sec(c + d*x))/(a + b*sec(c + d*x))
    F = ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_611():
    f = 1/((a + b*sec(c + d*x))*sqrt(sec(c + d*x)))
    F = ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_612():
    f = 1/((a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(-2) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_613():
    f = sec(c + d*x)**(sympy.S(9)/2)/(a + b*sec(c + d*x))**2
    F = ((Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_614():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_615():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**2
    F = ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_616():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_617():
    f = sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**2
    F = ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_618():
    f = 1/((a + b*sec(c + d*x))**2*sqrt(sec(c + d*x)))
    F = ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_619():
    f = 1/((a + b*sec(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(4))) + (Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(4))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_620():
    f = sec(c + d*x)**(sympy.S(9)/2)/(a + b*sec(c + d*x))**3
    F = ((Integer(-1) * (((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(11) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(11) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_621():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))**3
    F = ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_622():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**3
    F = ((((Symbol('a'))**(Integer(2)) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_623():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**3
    F = (Integer(-1) * ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_624():
    f = sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**3
    F = ((Integer(3) * Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(3) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Integer(5) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Symbol('b'))**(Integer(4))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_625():
    f = 1/((a + b*sec(c + d*x))**3*sqrt(sec(c + d*x)))
    F = ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Integer(35) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_626():
    f = 1/((a + b*sec(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(-1) * (Symbol('b') * ((Integer(24) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(6))) + (Integer(128) * (Symbol('a'))**(Integer(4)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(223) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(4)))) + (Integer(105) * (Symbol('b'))**(Integer(6)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(12) * (Symbol('a'))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Integer(63) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(5)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(13) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_627():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = ((Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_628():
    f = sqrt(a + b*sec(c + d*x))*sqrt(sec(c + d*x))
    F = ((Integer(2) * Symbol('a') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_629():
    f = sqrt(a + b*sec(c + d*x))/sqrt(sec(c + d*x))
    F = 2*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_630():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*b*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_631():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(15*a*d*sqrt(sec(c + d*x))) - 4*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a**2*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*a**2 - 4*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_632():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(35*a*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(50*a**2 - 8*b**2)*sin(c + d*x)/(105*a**2*d*sqrt(sec(c + d*x))) + 2*b*sqrt(a + b*sec(c + d*x))*(19*a**2 + 8*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(50*a**4 - 34*a**2*b**2 - 16*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(105*a**3*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_633():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(7) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(5) * Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_634():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_635():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_636():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 8*b*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_637():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*sqrt(sec(c + d*x))) + 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(5*a*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(6*a**2 + 2*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(5*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_638():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 16*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(50*a**2 + 6*b**2)*sin(c + d*x)/(105*a*d*sqrt(sec(c + d*x))) + 4*b*sqrt(a + b*sec(c + d*x))*(41*a**2 - 3*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(50*a**4 - 62*a**2*b**2 + 12*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(105*a**2*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_639():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = ((Symbol('b') * ((Integer(59) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(24) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(8) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_640():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = ((Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(11) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * Symbol('a') * Symbol('b') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(9) * Symbol('a') * Symbol('b') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_641():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_642():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(14) * Symbol('a') * Symbol('b') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_643():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 22*a*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(15*d*sqrt(sec(c + d*x))) + 16*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*a**2 + 46*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_644():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 6*a*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(10*a**2 + 18*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*b*sqrt(a + b*sec(c + d*x))*(29*a**2 + 3*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(21*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(10*a**4 - 4*a**2*b**2 - 6*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(21*a*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_645():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 38*a*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*sec(c + d*x))*(98*a**2 + 150*b**2)*sin(c + d*x)/(315*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*b*sqrt(a + b*sec(c + d*x))*(163*a**2 + 5*b**2)*sin(c + d*x)/(315*a*d*sqrt(sec(c + d*x))) + 4*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(57*a**4 - 62*a**2*b**2 + 5*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(315*a**2*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(294*a**4 + 558*a**2*b**2 - 20*b**4)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(315*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_646():
    f = sec(c + d*x)**(sympy.S(7)/2)/sqrt(a + b*sec(c + d*x))
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_647():
    f = sec(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*sec(c + d*x))
    F = ((sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_648():
    f = sec(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*sec(c + d*x))
    F = (Integer(2) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_649():
    f = sqrt(sec(c + d*x))/sqrt(a + b*sec(c + d*x))
    F = 2*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_650():
    f = 1/(sqrt(a + b*sec(c + d*x))*sqrt(sec(c + d*x)))
    F = -2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x))) + 2*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_651():
    f = 1/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) - 4*b*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 + 4*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_652():
    f = 1/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(5*a*d*sec(c + d*x)**(sympy.S(3)/2)) - 8*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)/(15*a**2*d*sqrt(sec(c + d*x))) - 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(7*a**2 + 8*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a**3*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*a**2 + 16*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_653():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_654():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(2) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_655():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - 2*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_656():
    f = sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*b*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + 2*b*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x))) + 2*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_657():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = 2*b**2*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - 4*b*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(a**2*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*a**2 - 4*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_658():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*b**2*sin(c + d*x)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*a**2 - 8*b**2)*sin(c + d*x)/(3*a**2*d*(a**2 - b**2)*sqrt(sec(c + d*x))) - 2*b*sqrt(a + b*sec(c + d*x))*(5*a**2 - 8*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 + 16*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**3*d*sqrt(a + b*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_659():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*b**2*sin(c + d*x)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(2*a**2 - 12*b**2)*sin(c + d*x)/(5*a**2*d*(a**2 - b**2)*sec(c + d*x)**(sympy.S(3)/2)) - 2*b*sqrt(a + b*sec(c + d*x))*(3*a**2 - 8*b**2)*sin(c + d*x)/(5*a**3*d*(a**2 - b**2)*sqrt(sec(c + d*x))) - 8*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 + 4*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(5*a**4*d*sqrt(a + b*sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(6*a**4 + 16*a**2*b**2 - 32*b**4)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(5*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_660():
    f = sec(c + d*x)**(sympy.S(9)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(3) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(9) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(3) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_661():
    f = sec(c + d*x)**(sympy.S(7)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_662():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*a*(a**2 - 5*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + 8*b*sqrt(a + b*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x))) + 2*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(d*sqrt(a + b*sec(c + d*x))*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_663():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + (4*a**2 + 4*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) - 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) - sqrt(a + b*sec(c + d*x))*(6*a**2 + 2*b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_664():
    f = sqrt(sec(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*b*sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) - 2*b*(5*a**2 - b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + 4*b*sqrt(a + b*sec(c + d*x))*(3*a**2 - b**2)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(6*a**2 - 4*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_665():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = 2*b**2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 8*b**2*(2*a**2 - b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) - 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(9*a**2 - 8*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(6*a**4 - 30*a**2*b**2 + 16*b**4)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_666():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*b**2*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(sec(c + d*x))) + 4*b**2*(5*a**2 - 3*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(sec(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*a**4 - 26*a**2*b**2 + 16*b**4)*sin(c + d*x)/(3*a**3*d*(a**2 - b**2)**2*sqrt(sec(c + d*x))) - 8*b*sqrt(a + b*sec(c + d*x))*(2*a**4 - 7*a**2*b**2 + 4*b**4)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**4 + 32*a**2*b**2 - 32*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(3*a**4*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_667():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*b**2*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sec(c + d*x)**(sympy.S(3)/2)) + 8*b**2*(3*a**2 - 2*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*sec(c + d*x))*(6*a**4 - 142*a**2*b**2 + 96*b**4)*sin(c + d*x)/(15*a**3*d*(a**2 - b**2)**2*sec(c + d*x)**(sympy.S(3)/2)) - 4*b*sqrt(a + b*sec(c + d*x))*(7*a**4 - 49*a**2*b**2 + 32*b**4)*sin(c + d*x)/(15*a**4*d*(a**2 - b**2)**2*sqrt(sec(c + d*x))) - 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(17*a**4 + 116*a**2*b**2 - 128*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))*sqrt(sec(c + d*x))/(15*a**5*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(18*a**6 + 110*a**4*b**2 - 424*a**2*b**4 + 256*b**6)*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**5*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_668():
    f = 1/(sqrt(3*sec(c + d*x) + 2)*sqrt(sec(c + d*x)))
    F = -3*sqrt(5)*sqrt(2*cos(c + d*x) + 3)*elliptic_f(c/2 + d*x/2, sympy.S(4)/5)*sqrt(sec(c + d*x))/(5*d*sqrt(3*sec(c + d*x) + 2)) + sqrt(5)*sqrt(3*sec(c + d*x) + 2)*elliptic_e(c/2 + d*x/2, sympy.S(4)/5)/(d*sqrt(2*cos(c + d*x) + 3)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_669():
    f = 1/(sqrt(3*sec(c + d*x) - 2)*sqrt(sec(c + d*x)))
    F = 3*sqrt(3 - 2*cos(c + d*x))*elliptic_f(c/2 + d*x/2, -4)*sqrt(sec(c + d*x))/(d*sqrt(3*sec(c + d*x) - 2)) - sqrt(3*sec(c + d*x) - 2)*elliptic_e(c/2 + d*x/2, -4)/(d*sqrt(3 - 2*cos(c + d*x))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_670():
    f = 1/(sqrt(2 - 3*sec(c + d*x))*sqrt(sec(c + d*x)))
    F = sqrt(2 - 3*sec(c + d*x))*elliptic_e(c/2 + d*x/2, -4)/(d*sqrt(3 - 2*cos(c + d*x))*sqrt(sec(c + d*x))) + 3*sqrt(3 - 2*cos(c + d*x))*elliptic_f(c/2 + d*x/2, -4)*sqrt(sec(c + d*x))/(d*sqrt(2 - 3*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_671():
    f = 1/(sqrt(-3*sec(c + d*x) - 2)*sqrt(sec(c + d*x)))
    F = -3*sqrt(5)*sqrt(2*cos(c + d*x) + 3)*elliptic_f(c/2 + d*x/2, sympy.S(4)/5)*sqrt(sec(c + d*x))/(5*d*sqrt(-3*sec(c + d*x) - 2)) - sqrt(5)*sqrt(-3*sec(c + d*x) - 2)*elliptic_e(c/2 + d*x/2, sympy.S(4)/5)/(d*sqrt(2*cos(c + d*x) + 3)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_672():
    f = 1/(sqrt(2*sec(c + d*x) + 3)*sqrt(sec(c + d*x)))
    F = -4*sqrt(5)*sqrt(3*cos(c + d*x) + 2)*elliptic_f(c/2 + d*x/2, sympy.S(6)/5)*sqrt(sec(c + d*x))/(15*d*sqrt(2*sec(c + d*x) + 3)) + 2*sqrt(5)*sqrt(2*sec(c + d*x) + 3)*elliptic_e(c/2 + d*x/2, sympy.S(6)/5)/(3*d*sqrt(3*cos(c + d*x) + 2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_673():
    f = 1/(sqrt(3 - 2*sec(c + d*x))*sqrt(sec(c + d*x)))
    F = 2*sqrt(3 - 2*sec(c + d*x))*elliptic_e(c/2 + d*x/2, 6)/(3*d*sqrt(3*cos(c + d*x) - 2)*sqrt(sec(c + d*x))) + 4*sqrt(3*cos(c + d*x) - 2)*elliptic_f(c/2 + d*x/2, 6)*sqrt(sec(c + d*x))/(3*d*sqrt(3 - 2*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_674():
    f = 1/(sqrt(2*sec(c + d*x) - 3)*sqrt(sec(c + d*x)))
    F = 4*sqrt(5)*sqrt(2 - 3*cos(c + d*x))*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(6)/5)*sqrt(sec(c + d*x))/(15*d*sqrt(2*sec(c + d*x) - 3)) - 2*sqrt(5)*sqrt(2*sec(c + d*x) - 3)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(6)/5)/(3*d*sqrt(2 - 3*cos(c + d*x))*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_675():
    f = 1/(sqrt(-2*sec(c + d*x) - 3)*sqrt(sec(c + d*x)))
    F = -4*sqrt(-3*cos(c + d*x) - 2)*elliptic_f(c/2 + d*x/2 + pi/2, 6)*sqrt(sec(c + d*x))/(3*d*sqrt(-2*sec(c + d*x) - 3)) - 2*sqrt(-2*sec(c + d*x) - 3)*elliptic_e(c/2 + d*x/2 + pi/2, 6)/(3*d*sqrt(-3*cos(c + d*x) - 2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_676():
    f = sqrt(sec(c + d*x))/sqrt(3*sec(c + d*x) + 2)
    F = 2*sqrt(5)*sqrt(2*cos(c + d*x) + 3)*elliptic_f(c/2 + d*x/2, sympy.S(4)/5)*sqrt(sec(c + d*x))/(5*d*sqrt(3*sec(c + d*x) + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_677():
    f = sqrt(sec(c + d*x))/sqrt(3*sec(c + d*x) - 2)
    F = 2*sqrt(3 - 2*cos(c + d*x))*elliptic_f(c/2 + d*x/2, -4)*sqrt(sec(c + d*x))/(d*sqrt(3*sec(c + d*x) - 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_678():
    f = sqrt(sec(c + d*x))/sqrt(2 - 3*sec(c + d*x))
    F = 2*sqrt(3 - 2*cos(c + d*x))*elliptic_f(c/2 + d*x/2, -4)*sqrt(sec(c + d*x))/(d*sqrt(2 - 3*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_679():
    f = sqrt(sec(c + d*x))/sqrt(-3*sec(c + d*x) - 2)
    F = 2*sqrt(5)*sqrt(2*cos(c + d*x) + 3)*elliptic_f(c/2 + d*x/2, sympy.S(4)/5)*sqrt(sec(c + d*x))/(5*d*sqrt(-3*sec(c + d*x) - 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_680():
    f = sqrt(sec(c + d*x))/sqrt(2*sec(c + d*x) + 3)
    F = 2*sqrt(5)*sqrt(3*cos(c + d*x) + 2)*elliptic_f(c/2 + d*x/2, sympy.S(6)/5)*sqrt(sec(c + d*x))/(5*d*sqrt(2*sec(c + d*x) + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_681():
    f = sqrt(sec(c + d*x))/sqrt(3 - 2*sec(c + d*x))
    F = 2*sqrt(3*cos(c + d*x) - 2)*elliptic_f(c/2 + d*x/2, 6)*sqrt(sec(c + d*x))/(d*sqrt(3 - 2*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_682():
    f = sqrt(sec(c + d*x))/sqrt(2*sec(c + d*x) - 3)
    F = 2*sqrt(5)*sqrt(2 - 3*cos(c + d*x))*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(6)/5)*sqrt(sec(c + d*x))/(5*d*sqrt(2*sec(c + d*x) - 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_683():
    f = sqrt(sec(c + d*x))/sqrt(-2*sec(c + d*x) - 3)
    F = 2*sqrt(-3*cos(c + d*x) - 2)*elliptic_f(c/2 + d*x/2 + pi/2, 6)*sqrt(sec(c + d*x))/(d*sqrt(-2*sec(c + d*x) - 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_684():
    f = (a + b*sec(c + d*x))**(sympy.S(1)/3)*sec(c + d*x)
    F = sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_685():
    f = (a + b*sec(c + d*x))**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_686():
    f = (a + b*sec(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)**4
    F = -9*a*(a + b*sec(c + d*x))**(sympy.S(5)/3)*tan(c + d*x)/(44*b**2*d) + sqrt(2)*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(18*a**2 + 49*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(220*b**3*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(5)/3)*tan(c + d*x)*sec(c + d*x)/(11*b*d) + (a + b*sec(c + d*x))**(sympy.S(2)/3)*(27*a**2 + 96*b**2)*tan(c + d*x)/(220*b**2*d) - sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(9*a**4 + 23*a**2*b**2 - 32*b**4)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(110*b**3*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_687():
    f = (a + b*sec(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)**3
    F = -9*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)/(40*b*d) + 3*sqrt(2)*a*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(a**2 - b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(20*b**2*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(5)/3)*tan(c + d*x)/(8*b*d) - sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(6*a**2 - 25*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(40*b**2*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_688():
    f = (a + b*sec(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)**2
    F = 2*sqrt(2)*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(5*b*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)/(5*d) - 2*sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(a**2 - b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(5*b*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_689():
    f = (a + b*sec(c + d*x))**(sympy.S(2)/3)*sec(c + d*x)
    F = sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_690():
    f = (a + b*sec(c + d*x))**(sympy.S(2)/3)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_691():
    f = (a + b*sec(c + d*x))**(sympy.S(4)/3)*sec(c + d*x)
    F = sqrt(2)*(a + b)*(a + b*sec(c + d*x))**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_692():
    f = (a + b*sec(c + d*x))**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_693():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/3)*sec(c + d*x)**4
    F = -9*a*(a + b*sec(c + d*x))**(sympy.S(8)/3)*tan(c + d*x)/(77*b**2*d) + 3*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(18*a**2 + 97*b**2)*tan(c + d*x)/(1232*b**2*d) - sqrt(2)*a*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(18*a**4 + 79*a**2*b**2 - 97*b**4)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(616*b**3*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(8)/3)*tan(c + d*x)*sec(c + d*x)/(14*b*d) + (a + b*sec(c + d*x))**(sympy.S(5)/3)*(54*a**2 + 363*b**2)*tan(c + d*x)/(1232*b**2*d) + sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(36*a**4 + 164*a**2*b**2 + 605*b**4)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(1232*b**3*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_694():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/3)*sec(c + d*x)**3
    F = -9*a*(a + b*sec(c + d*x))**(sympy.S(5)/3)*tan(c + d*x)/(88*b*d) - sqrt(2)*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(30*a**2 - 373*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(440*b**2*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(8)/3)*tan(c + d*x)/(11*b*d) - (a + b*sec(c + d*x))**(sympy.S(2)/3)*(45*a**2 - 192*b**2)*tan(c + d*x)/(440*b*d) + sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(15*a**4 - 79*a**2*b**2 + 64*b**4)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(220*b**2*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_695():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/3)*sec(c + d*x)**2
    F = 3*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)/(8*d) - sqrt(2)*a*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(a**2 - b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(4*b*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(5)/3)*tan(c + d*x)/(8*d) + sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(2*a**2 + 5*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(8*b*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_696():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/3)*sec(c + d*x)
    F = sqrt(2)*(a + b)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_697():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(3))**(Integer(-1)))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_698():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))**(sympy.S(1)/3)
    F = -9*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)/(20*b**2*d) - sqrt(2)*a*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(9*a**2 + 11*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(20*b**3*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)*sec(c + d*x)/(8*b*d) + sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(18*a**2 + 25*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(40*b**3*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_699():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(1)/3)
    F = -3*sqrt(2)*a*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(5*b**2*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1)) + 3*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)/(5*b*d) + sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(3*a**2 + 2*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(5*b**2*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_700():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(1)/3)
    F = -sqrt(2)*a*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(b*d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1)) + sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(b*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_701():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(1)/3)
    F = sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(d*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_702():
    f = (a + b*sec(c + d*x))**(sympy.S(-1)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3))**(Integer(-1))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_703():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(2)/3)
    F = sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(d*(a + b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_704():
    f = (a + b*sec(c + d*x))**(sympy.S(-2)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_705():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(4)/3)
    F = sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(4)/3, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(d*(a + b)*(a + b*sec(c + d*x))**(sympy.S(1)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_706():
    f = (a + b*sec(c + d*x))**(sympy.S(-4)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_707():
    f = sec(c + d*x)**4/(a + b*sec(c + d*x))**(sympy.S(5)/3)
    F = -3*a**2*tan(c + d*x)*sec(c + d*x)/(2*b*d*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(a**2 - b**2)) - sqrt(2)*a*(a + b*sec(c + d*x))**(sympy.S(1)/3)*(9*a**2 - 7*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(4*b**3*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(a**2 - b**2)*sqrt(sec(c + d*x) + 1)) + (a + b*sec(c + d*x))**(sympy.S(1)/3)*(9*a**2 - 3*b**2)*tan(c + d*x)/(4*b**2*d*(a**2 - b**2)) + sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*(9*a**4 - 10*a**2*b**2 - b**4)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(4*b**3*d*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(a**2 - b**2)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_708():
    f = sec(c + d*x)**3/(a + b*sec(c + d*x))**(sympy.S(5)/3)
    F = -3*a**2*tan(c + d*x)/(2*b*d*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(a**2 - b**2)) - sqrt(2)*a*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*(3*a**2 - 4*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(2*b**2*d*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(a**2 - b**2)*sqrt(sec(c + d*x) + 1)) + sqrt(2)*(a + b*sec(c + d*x))**(sympy.S(1)/3)*(3*a**2 - 2*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(2*b**2*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(a**2 - b**2)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_709():
    f = sec(c + d*x)**2/(a + b*sec(c + d*x))**(sympy.S(5)/3)
    F = 3*a*tan(c + d*x)/(d*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(2*a**2 - 2*b**2)) - sqrt(2)*a*(a + b*sec(c + d*x))**(sympy.S(1)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - sec(c + d*x))/(a + b), sympy.S.Half - sec(c + d*x)/2)/(2*b*d*((a + b*sec(c + d*x))/(a + b))**(sympy.S(1)/3)*(a**2 - b**2)*sqrt(sec(c + d*x) + 1)) + sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*(a**2 - 2*b**2)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(2*b*d*(a + b*sec(c + d*x))**(sympy.S(2)/3)*(a**2 - b**2)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_710():
    f = sec(c + d*x)/(a + b*sec(c + d*x))**(sympy.S(5)/3)
    F = sqrt(2)*((a + b*sec(c + d*x))/(a + b))**(sympy.S(2)/3)*tan(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(5)/3, sympy.S(3)/2, sympy.S.Half - sec(c + d*x)/2, b*(1 - sec(c + d*x))/(a + b))/(d*(a + b)*(a + b*sec(c + d*x))**(sympy.S(2)/3)*sqrt(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_711():
    f = (a + b*sec(c + d*x))**(sympy.S(-5)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(3))**(Integer(-1)))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_712():
    f = sec(c + d*x)**(sympy.S(2)/3)/(a + b*sec(c + d*x))
    F = a*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*(cos(c + d*x)**2)**(sympy.S(1)/6)*sec(c + d*x)**(sympy.S(1)/3)) - b*(cos(c + d*x)**2)**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))*sec(c + d*x)**(sympy.S(2)/3)/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_713():
    f = sec(c + d*x)**(sympy.S(1)/3)/(a + b*sec(c + d*x))
    F = a*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*(cos(c + d*x)**2)**(sympy.S(1)/3)*sec(c + d*x)**(sympy.S(2)/3)) - b*(cos(c + d*x)**2)**(sympy.S(1)/6)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))*sec(c + d*x)**(sympy.S(1)/3)/(d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_714():
    f = 1/((a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(1)/3))
    F = a*(cos(c + d*x)**2)**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))*sec(c + d*x)**(sympy.S(2)/3)/(d*(a**2 - b**2)) - b*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*(cos(c + d*x)**2)**(sympy.S(1)/6)*sec(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_715():
    f = 1/((a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(2)/3))
    F = a*(cos(c + d*x)**2)**(sympy.S(1)/6)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))*sec(c + d*x)**(sympy.S(1)/3)/(d*(a**2 - b**2)) - b*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, a**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*(cos(c + d*x)**2)**(sympy.S(1)/3)*sec(c + d*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_716():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_717():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_718():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_719():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(2)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_720():
    f = sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_721():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_722():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(2)/3)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_723():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_724():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_725():
    f = sqrt(a + b*sec(c + d*x))/sec(c + d*x)**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_726():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_727():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_728():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_729():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(2)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_730():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_731():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_732():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(2)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_733():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_734():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_735():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_736():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_737():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_738():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_739():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(2)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_740():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_741():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(1)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_742():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(2)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_743():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_744():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_745():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_746():
    f = sec(c + d*x)**(sympy.S(7)/3)/sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_747():
    f = sec(c + d*x)**(sympy.S(5)/3)/sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_748():
    f = sec(c + d*x)**(sympy.S(4)/3)/sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_749():
    f = sec(c + d*x)**(sympy.S(2)/3)/sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_750():
    f = sec(c + d*x)**(sympy.S(1)/3)/sqrt(a + b*sec(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_751():
    f = 1/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(1)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_752():
    f = 1/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(2)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_753():
    f = 1/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(4)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_754():
    f = 1/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(5)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_755():
    f = 1/(sqrt(a + b*sec(c + d*x))*sec(c + d*x)**(sympy.S(7)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_756():
    f = sec(c + d*x)**(sympy.S(7)/3)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_757():
    f = sec(c + d*x)**(sympy.S(5)/3)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_758():
    f = sec(c + d*x)**(sympy.S(4)/3)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_759():
    f = sec(c + d*x)**(sympy.S(2)/3)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_760():
    f = sec(c + d*x)**(sympy.S(1)/3)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_761():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(1)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_762():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(2)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_763():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(4)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_764():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_765():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(7)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_766():
    f = sec(c + d*x)**(sympy.S(7)/3)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_767():
    f = sec(c + d*x)**(sympy.S(5)/3)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_768():
    f = sec(c + d*x)**(sympy.S(4)/3)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_769():
    f = sec(c + d*x)**(sympy.S(2)/3)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_770():
    f = sec(c + d*x)**(sympy.S(1)/3)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_771():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(1)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_772():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(2)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_773():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(4)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_774():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_775():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/3))
    F = sympy.Function('Unintegrable')((((sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_776():
    f = (d*sec(e + f*x))**n*(a + b*sec(e + f*x))**3
    F = a*b**2*(d*sec(e + f*x))**n*(2*n + 5)*tan(e + f*x)/(f*(n + 1)*(n + 2)) - a*d*(d*sec(e + f*x))**(n - 1)*(a**2*(n + 1) + 3*b**2*n)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n**2)*sqrt(sin(e + f*x)**2)) + b**2*(d*sec(e + f*x))**n*(a + b*sec(e + f*x))*tan(e + f*x)/(f*(n + 2)) + b*(d*sec(e + f*x))**n*(3*a**2*(n + 2) + b**2*(n + 1))*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)/(f*n*(n + 2)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_777():
    f = (d*sec(e + f*x))**n*(a + b*sec(e + f*x))**2
    F = 2*a*b*(d*sec(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)/(f*n*sqrt(sin(e + f*x)**2)) + b**2*(d*sec(e + f*x))**n*tan(e + f*x)/(f*(n + 1)) - d*(d*sec(e + f*x))**(n - 1)*(a**2*(n + 1) + b**2*n)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n**2)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_778():
    f = (d*sec(e + f*x))**n*(a + b*sec(e + f*x))
    F = -a*d*(d*sec(e + f*x))**(n - 1)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), cos(e + f*x)**2)/(f*(1 - n)*sqrt(sin(e + f*x)**2)) + b*(d*sec(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, -n/2), (1 - n/2,), cos(e + f*x)**2)/(f*n*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_779():
    f = (d*sec(e + f*x))**n/(a + b*sec(e + f*x))
    F = a*(d*sec(e + f*x))**n*(cos(e + f*x)**2)**(n/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 1, n/2 + sympy.S(-1)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)) - b*(d*sec(e + f*x))**n*(cos(e + f*x)**2)**(n/2)*sin(e + f*x)*appellf1(sympy.S.Half, 1, n/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_780():
    f = (d*sec(e + f*x))**n/(a + b*sec(e + f*x))**2
    F = a**2*(d*sec(e + f*x))**n*(cos(e + f*x)**2)**(n/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, n/2 + sympy.S(-3)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2) - 2*a*b*(d*sec(e + f*x))**n*(cos(e + f*x)**2)**(n/2)*sin(e + f*x)*appellf1(sympy.S.Half, 2, n/2 - 1, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2) + b**2*(d*sec(e + f*x))**n*(cos(e + f*x)**2)**(n/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, n/2 + sympy.S(-1)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_781():
    f = (d*sec(e + f*x))**n*(a + b*sec(e + f*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_782():
    f = (d*sec(e + f*x))**n*sqrt(a + b*sec(e + f*x))
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_783():
    f = (d*sec(e + f*x))**n/sqrt(a + b*sec(e + f*x))
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_784():
    f = (d*sec(e + f*x))**n/(a + b*sec(e + f*x))**(sympy.S(3)/2)
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')) * (((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_785():
    f = (a + b*sec(e + f*x))**m*sec(e + f*x)**n
    F = sympy.Function('Unintegrable')(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Symbol('n')) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_786():
    f = (d*sec(e + f*x))**n*(a + b*sec(e + f*x))**m
    F = sympy.Function('Unintegrable')((((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('n')) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_787():
    f = (a + b*sec(e + f*x))**m*sec(e + f*x)**3
    F = -sqrt(2)*a*(a + b)*(a + b*sec(e + f*x))**m*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, b*(1 - sec(e + f*x))/(a + b))/(b**2*f*((a + b*sec(e + f*x))/(a + b))**m*(m + 2)*sqrt(sec(e + f*x) + 1)) + (a + b*sec(e + f*x))**(m + 1)*tan(e + f*x)/(b*f*(m + 2)) + sqrt(2)*(a + b*sec(e + f*x))**m*(a**2 + b**2*(m + 1))*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, b*(1 - sec(e + f*x))/(a + b))/(b**2*f*((a + b*sec(e + f*x))/(a + b))**m*(m + 2)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_788():
    f = (a + b*sec(e + f*x))**m*sec(e + f*x)**2
    F = -sqrt(2)*a*(a + b*sec(e + f*x))**m*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, b*(1 - sec(e + f*x))/(a + b))/(b*f*((a + b*sec(e + f*x))/(a + b))**m*sqrt(sec(e + f*x) + 1)) + sqrt(2)*(a + b)*(a + b*sec(e + f*x))**m*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m - 1, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, b*(1 - sec(e + f*x))/(a + b))/(b*f*((a + b*sec(e + f*x))/(a + b))**m*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_789():
    f = (a + b*sec(e + f*x))**m*sec(e + f*x)
    F = sqrt(2)*(a + b*sec(e + f*x))**m*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S.Half, -m, sympy.S(3)/2, sympy.S.Half - sec(e + f*x)/2, b*(1 - sec(e + f*x))/(a + b))/(f*((a + b*sec(e + f*x))/(a + b))**m*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_790():
    f = (a + b*sec(e + f*x))**m
    F = sympy.Function('Unintegrable')(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m')), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_791():
    f = (a + b*sec(e + f*x))**m*cos(e + f*x)
    F = sympy.Function('Unintegrable')((sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_792():
    f = (a + b*sec(e + f*x))**m*cos(e + f*x)**2
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m'))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_793():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 14*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 14*a*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 2*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 10*b*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 10*b*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_794():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 10*a*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 10*a*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*b*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_795():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*b*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*b*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_796():
    f = (a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*b*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_797():
    f = (a + b*sec(c + d*x))*sqrt(cos(c + d*x))
    F = 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*b*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_798():
    f = (a + b*sec(c + d*x))/sqrt(cos(c + d*x))
    F = 2*a*elliptic_f(c/2 + d*x/2, 2)/d + 2*b*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*b*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_799():
    f = (a + b*sec(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*b*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*b*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_800():
    f = (a + b*sec(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 6*b*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 2*b*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 6*b*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_801():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 4*a*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 20*a*b*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 20*a*b*elliptic_f(c/2 + d*x/2, 2)/(21*d) + (14*a**2 + 18*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + (14*a**2 + 18*b**2)*elliptic_e(c/2 + d*x/2, 2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_802():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 4*a*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 12*a*b*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (10*a**2 + 14*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (10*a**2 + 14*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_803():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 4*a*b*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 4*a*b*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (6*a**2 + 10*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_804():
    f = (a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 4*a*b*elliptic_e(c/2 + d*x/2, 2)/d + (2*a**2 + 6*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_805():
    f = (a + b*sec(c + d*x))**2*sqrt(cos(c + d*x))
    F = 4*a*b*elliptic_f(c/2 + d*x/2, 2)/d + 2*b**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + (2*a**2 - 2*b**2)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_806():
    f = (a + b*sec(c + d*x))**2/sqrt(cos(c + d*x))
    F = 4*a*b*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 4*a*b*elliptic_e(c/2 + d*x/2, 2)/d + 2*b**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + (6*a**2 + 2*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_807():
    f = (a + b*sec(c + d*x))**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 4*a*b*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a*b*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*b**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + (10*a**2 + 6*b**2)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - (10*a**2 + 6*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_808():
    f = (a + b*sec(c + d*x))**2/cos(c + d*x)**(sympy.S(5)/2)
    F = 12*a*b*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 4*a*b*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 12*a*b*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*b**2*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + (14*a**2 + 10*b**2)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + (14*a**2 + 10*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_809():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(9)/2)
    F = 40*a**2*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d) + 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 2*a*(7*a**2 + 27*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 2*a*(7*a**2 + 27*b**2)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 2*b*(15*a**2 + 7*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*b*(15*a**2 + 7*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_810():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(7)/2)
    F = 32*a**2*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*(5*a**2 + 21*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*a*(5*a**2 + 21*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*b*(9*a**2 + 5*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_811():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(5)/2)
    F = 8*a**2*b*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*a*(a**2 + 5*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*b*(a**2 + b**2)*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_812():
    f = (a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*(a**2 + 9*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d) - 2*b*(a**2 - 3*b**2)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x))) + 2*b*(3*a**2 - b**2)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_813():
    f = (a + b*sec(c + d*x))**3*sqrt(cos(c + d*x))
    F = 16*a*b**2*sin(c + d*x)/(3*d*sqrt(cos(c + d*x))) + 2*a*(a**2 - 3*b**2)*elliptic_e(c/2 + d*x/2, 2)/d + 2*b**2*(a + b*sec(c + d*x))*sin(c + d*x)/(3*d*sqrt(cos(c + d*x))) + 2*b*(9*a**2 + b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_814():
    f = (a + b*sec(c + d*x))**3/sqrt(cos(c + d*x))
    F = 8*a*b**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(a**2 + b**2)*elliptic_f(c/2 + d*x/2, 2)/d + 2*b**2*(a + b*sec(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(3)/2)) + 6*b*(5*a**2 + b**2)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 6*b*(5*a**2 + b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_815():
    f = (a + b*sec(c + d*x))**3/cos(c + d*x)**(sympy.S(3)/2)
    F = 32*a*b**2*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a*(5*a**2 + 9*b**2)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 2*a*(5*a**2 + 9*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*b**2*(a + b*sec(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*b*(21*a**2 + 5*b**2)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*b*(21*a**2 + 5*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_816():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))
    F = ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(4)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_817():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))
    F = ((Integer(-2) * Symbol('b') * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_818():
    f = sqrt(cos(c + d*x))/(a + b*sec(c + d*x))
    F = ((Integer(2) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_819():
    f = 1/((a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    F = ((Integer(2) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_820():
    f = 1/((a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(2) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_821():
    f = 1/((a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(-2) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_822():
    f = 1/((a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = ((Integer(2) * Symbol('a') * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_823():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**2
    F = (Integer(-1) * ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(4))) + (Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(4))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(4)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_824():
    f = sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**2
    F = ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_825():
    f = 1/((a + b*sec(c + d*x))**2*sqrt(cos(c + d*x)))
    F = ((Symbol('b') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_826():
    f = 1/((a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * (sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_827():
    f = 1/((a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_828():
    f = 1/((a + b*sec(c + d*x))**2*cos(c + d*x)**(sympy.S(7)/2))
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_829():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**3
    F = (Integer(-1) * ((Symbol('b') * ((Integer(24) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(6))) + (Integer(128) * (Symbol('a'))**(Integer(4)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(223) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(4)))) + (Integer(105) * (Symbol('b'))**(Integer(6)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(12) * (Symbol('a'))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(3)) * ((Integer(63) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(5)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(13) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_830():
    f = sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**3
    F = ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * ((Integer(35) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_831():
    f = 1/((a + b*sec(c + d*x))**3*sqrt(cos(c + d*x)))
    F = ((Integer(3) * Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(3) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * ((Integer(5) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Symbol('b'))**(Integer(4))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_832():
    f = 1/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_833():
    f = 1/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(5)/2))
    F = ((((Symbol('a'))**(Integer(2)) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_834():
    f = 1/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(7)/2))
    F = ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_835():
    f = 1/((a + b*sec(c + d*x))**3*cos(c + d*x)**(sympy.S(9)/2))
    F = (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(11) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(11) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_836():
    f = sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d) - 4*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*a**2 - 4*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_837():
    f = sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*b*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_838():
    f = sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))
    F = 2*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_839():
    f = sqrt(a + b*sec(c + d*x))/sqrt(cos(c + d*x))
    F = ((Integer(2) * Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_840():
    f = sqrt(a + b*sec(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = ((Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_841():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 16*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + sqrt(a + b*sec(c + d*x))*(50*a**2 + 6*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*a*d) + 4*b*sqrt(a + b*sec(c + d*x))*(41*a**2 - 3*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(105*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(50*a**4 - 62*a**2*b**2 + 12*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(105*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_842():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 4*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(5*a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(6*a**2 + 2*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(5*a*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_843():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 8*b*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*d*sqrt((a*cos(c + d*x) + b)/(a + b))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_844():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_845():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_846():
    f = (a + b*sec(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(7) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_847():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 38*a*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d) + sqrt(a + b*sec(c + d*x))*(98*a**2 + 150*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(315*d) + 2*b*sqrt(a + b*sec(c + d*x))*(163*a**2 + 5*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(315*a*d) + 4*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(57*a**4 - 62*a**2*b**2 + 5*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(315*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(294*a**4 + 558*a**2*b**2 - 20*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(315*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_848():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 6*a*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(7*d) + sqrt(a + b*sec(c + d*x))*(10*a**2 + 18*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*b*sqrt(a + b*sec(c + d*x))*(29*a**2 + 3*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(21*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(10*a**4 - 4*a**2*b**2 - 6*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(21*a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_849():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 22*a*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d) + 16*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*a**2 + 46*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_850():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(14) * Symbol('a') * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_851():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_852():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = ((Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(11) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * Symbol('a') * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(9) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_853():
    f = (a + b*sec(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = ((Symbol('b') * ((Integer(59) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(5) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_854():
    f = cos(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*sec(c + d*x))
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a*d) - 8*b*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a**2*d) - 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(7*a**2 + 8*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(15*a**3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(18*a**2 + 16*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(15*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_855():
    f = cos(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*sec(c + d*x))
    F = 2*sqrt(a + b*sec(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) - 4*b*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 + 4*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_856():
    f = sqrt(cos(c + d*x))/sqrt(a + b*sec(c + d*x))
    F = -2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_857():
    f = 1/(sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    F = 2*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_858():
    f = 1/(sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(2) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_859():
    f = 1/(sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = ((sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_860():
    f = 1/(sqrt(a + b*sec(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = ((Integer(-1) * (Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_861():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(2*a**2 - 12*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a**2*d*(a**2 - b**2)) - 2*b*sqrt(a + b*sec(c + d*x))*(3*a**2 - 8*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*a**3*d*(a**2 - b**2)) - 8*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 + 4*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(5*a**4*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(6*a**4 + 16*a**2*b**2 - 32*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(5*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_862():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sin(c + d*x)*sqrt(cos(c + d*x))/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)) + sqrt(a + b*sec(c + d*x))*(2*a**2 - 8*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*(a**2 - b**2)) - 2*b*sqrt(a + b*sec(c + d*x))*(5*a**2 - 8*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**2 + 16*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_863():
    f = sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sin(c + d*x)/(a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - 4*b*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(a**2*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(2*a**2 - 4*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_864():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = -2*b*sin(c + d*x)/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + 2*b*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)) + 2*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(a*d*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_865():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*a*sin(c + d*x)/(d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - 2*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_866():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(2) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_867():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = ((sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_868():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*b**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 4*b**2*(5*a**2 - 3*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2) + sqrt(a + b*sec(c + d*x))*(2*a**4 - 26*a**2*b**2 + 16*b**4)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**3*d*(a**2 - b**2)**2) - 8*b*sqrt(a + b*sec(c + d*x))*(2*a**4 - 7*a**2*b**2 + 4*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**4*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2) + sqrt((a*cos(c + d*x) + b)/(a + b))*(2*a**4 + 32*a**2*b**2 - 32*b**4)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**4*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_869():
    f = sqrt(cos(c + d*x))/(a + b*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*b**2*sin(c + d*x)/(3*a*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(cos(c + d*x))) + 8*b**2*(2*a**2 - b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) - 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*(9*a**2 - 8*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + sqrt(a + b*sec(c + d*x))*(6*a**4 - 30*a**2*b**2 + 16*b**4)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_870():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = -2*b*sin(c + d*x)/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*sqrt(cos(c + d*x))) - 2*b*(5*a**2 - b**2)*sin(c + d*x)/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + 4*b*sqrt(a + b*sec(c + d*x))*(3*a**2 - b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2) + sqrt((a*cos(c + d*x) + b)/(a + b))*(6*a**2 - 4*b**2)*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a**2*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_871():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*a*sin(c + d*x)/(d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*sqrt(cos(c + d*x))) + (4*a**2 + 4*b**2)*sin(c + d*x)/(3*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) - 2*b*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - sqrt(a + b*sec(c + d*x))*(6*a**2 + 2*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*a*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_872():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -2*a**2*sin(c + d*x)/(3*b*d*(a + b*sec(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(cos(c + d*x))) + 2*a*(a**2 - 5*b**2)*sin(c + d*x)/(3*b*d*sqrt(a + b*sec(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + 8*b*sqrt(a + b*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*a/(a + b))/(3*d*sqrt((a*cos(c + d*x) + b)/(a + b))*(a**2 - b**2)**2) + 2*sqrt((a*cos(c + d*x) + b)/(a + b))*elliptic_f(c/2 + d*x/2, 2*a/(a + b))/(d*sqrt(a + b*sec(c + d*x))*(3*a**2 - 3*b**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_873():
    f = 1/((a + b*sec(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = ((Integer(-2) * Symbol('a') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_874():
    f = (d*cos(e + f*x))**n*(a + b*sec(e + f*x))**3
    F = a*b**2*(d*cos(e + f*x))**n*(5 - 2*n)*tan(e + f*x)/(f*(1 - n)*(2 - n)) - a*(d*cos(e + f*x))**n*(a**2*(1 - n) - 3*b**2*n)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(1 - n)*(n + 1)*sqrt(sin(e + f*x)**2)) + b**2*(d*cos(e + f*x))**n*(a + b*sec(e + f*x))*tan(e + f*x)/(f*(2 - n)) - b*(d*cos(e + f*x))**n*(3*a**2*(2 - n) + b**2*(1 - n))*sin(e + f*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(e + f*x)**2)/(f*n*(2 - n)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_875():
    f = (d*cos(e + f*x))**n*(a + b*sec(e + f*x))**2
    F = -2*a*b*(d*cos(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(e + f*x)**2)/(f*n*sqrt(sin(e + f*x)**2)) + b**2*(d*cos(e + f*x))**n*tan(e + f*x)/(f*(1 - n)) - (d*cos(e + f*x))**n*(a**2*(1 - n) - b**2*n)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(1 - n)*(n + 1)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_876():
    f = (d*cos(e + f*x))**n*(a + b*sec(e + f*x))
    F = -a*(d*cos(e + f*x))**(n + 1)*sin(e + f*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(d*f*(n + 1)*sqrt(sin(e + f*x)**2)) - b*(d*cos(e + f*x))**n*sin(e + f*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(e + f*x)**2)/(f*n*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_877():
    f = (d*cos(e + f*x))**n/(a + b*sec(e + f*x))
    F = a*(d*cos(e + f*x))**n*(cos(e + f*x)**2)**(-n/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 1, -n/2 + sympy.S(-1)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)) - b*(d*cos(e + f*x))**n*sin(e + f*x)*appellf1(sympy.S.Half, 1, -n/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)*(cos(e + f*x)**2)**(n/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_1_2_d_sec_pow_n_a_plus_b_sec_pow_m_878():
    f = (d*cos(e + f*x))**n/(a + b*sec(e + f*x))**2
    F = a**2*(d*cos(e + f*x))**n*(cos(e + f*x)**2)**(-n/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, -n/2 + sympy.S(-3)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2) - 2*a*b*(d*cos(e + f*x))**n*sin(e + f*x)*appellf1(sympy.S.Half, 2, -n/2 - 1, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2*(cos(e + f*x)**2)**(n/2)) + b**2*(d*cos(e + f*x))**n*(cos(e + f*x)**2)**(-n/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, -n/2 + sympy.S(-1)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2)
    assert integrate(f, x) == F

