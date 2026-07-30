"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.2.1 (a+b cos)^m (c+d cos)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, n = symbols('A B a b c d e f m n')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_1():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)**5
    F = 5*a*x/16 + a*sin(c + d*x)**5/(5*d) - 2*a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a*sin(c + d*x)*cos(c + d*x)/(16*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_2():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)**4
    F = 3*a*x/8 + a*sin(c + d*x)**5/(5*d) - 2*a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_3():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)**3
    F = 3*a*x/8 - a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_4():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)**2
    F = a*x/2 - a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)*cos(c + d*x)/(2*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_5():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_6():
    f = a*cos(c + d*x) + a
    F = a*x + a*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_7():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)
    F = a*x + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_8():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**2
    F = a*tan(c + d*x)/d + a*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_9():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**3
    F = a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*tan(c + d*x)/d + a*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_10():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*tan(c + d*x)/d + a*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_11():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**5
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + a*tan(c + d*x)/d + 3*a*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_12():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**6
    F = a*tan(c + d*x)**5/(5*d) + 2*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + a*tan(c + d*x)/d + 3*a*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_13():
    f = (a*cos(c + d*x) + a)**2*cos(c + d*x)**4
    F = 11*a**2*x/16 + 2*a**2*sin(c + d*x)**5/(5*d) - 4*a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 11*a**2*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 11*a**2*sin(c + d*x)*cos(c + d*x)/(16*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_14():
    f = (a*cos(c + d*x) + a)**2*cos(c + d*x)**3
    F = 3*a**2*x/4 + a**2*sin(c + d*x)**5/(5*d) - a**2*sin(c + d*x)**3/d + a**2*sin(c + d*x)*cos(c + d*x)**3/(2*d) + 3*a**2*sin(c + d*x)*cos(c + d*x)/(4*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_15():
    f = (a*cos(c + d*x) + a)**2*cos(c + d*x)**2
    F = 7*a**2*x/8 - 2*a**2*sin(c + d*x)**3/(3*d) + a**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 7*a**2*sin(c + d*x)*cos(c + d*x)/(8*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_16():
    f = (a*cos(c + d*x) + a)**2
    F = 3*a**2*x/2 + a**2*sin(c + d*x)*cos(c + d*x)/(2*d) + 2*a**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_17():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)
    F = 2*a**2*x + a**2*sin(c + d*x)/d + a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_18():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)**2
    F = a**2*x + a**2*tan(c + d*x)/d + 2*a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_19():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)**3
    F = a**2*tan(c + d*x)*sec(c + d*x)/(2*d) + 2*a**2*tan(c + d*x)/d + 3*a**2*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_20():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)**4
    F = a**2*tan(c + d*x)**3/(3*d) + a**2*tan(c + d*x)*sec(c + d*x)/d + 2*a**2*tan(c + d*x)/d + a**2*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_21():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)**5
    F = 2*a**2*tan(c + d*x)**3/(3*d) + a**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 7*a**2*tan(c + d*x)*sec(c + d*x)/(8*d) + 2*a**2*tan(c + d*x)/d + 7*a**2*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_22():
    f = (a*cos(c + d*x) + a)**3*cos(c + d*x)**3
    F = 23*a**3*x/16 + 3*a**3*sin(c + d*x)**5/(5*d) - 7*a**3*sin(c + d*x)**3/(3*d) + a**3*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 23*a**3*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 23*a**3*sin(c + d*x)*cos(c + d*x)/(16*d) + 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_23():
    f = (a*cos(c + d*x) + a)**3*cos(c + d*x)**2
    F = 13*a**3*x/8 + a**3*sin(c + d*x)**5/(5*d) - 5*a**3*sin(c + d*x)**3/(3*d) + 3*a**3*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 13*a**3*sin(c + d*x)*cos(c + d*x)/(8*d) + 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_24():
    f = (a*cos(c + d*x) + a)**3
    F = 5*a**3*x/2 - a**3*sin(c + d*x)**3/(3*d) + 3*a**3*sin(c + d*x)*cos(c + d*x)/(2*d) + 4*a**3*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_25():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)
    F = 7*a**3*x/2 + a**3*sin(c + d*x)*cos(c + d*x)/(2*d) + 3*a**3*sin(c + d*x)/d + a**3*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_26():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**2
    F = 3*a**3*x + a**3*sin(c + d*x)/d + a**3*tan(c + d*x)/d + 3*a**3*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_27():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**3
    F = a**3*x + a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 3*a**3*tan(c + d*x)/d + 7*a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_28():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**4
    F = a**3*tan(c + d*x)**3/(3*d) + 3*a**3*tan(c + d*x)*sec(c + d*x)/(2*d) + 4*a**3*tan(c + d*x)/d + 5*a**3*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_29():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**5
    F = a**3*tan(c + d*x)**3/d + a**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 15*a**3*tan(c + d*x)*sec(c + d*x)/(8*d) + 4*a**3*tan(c + d*x)/d + 15*a**3*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_30():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**6
    F = a**3*tan(c + d*x)**5/(5*d) + 5*a**3*tan(c + d*x)**3/(3*d) + 3*a**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 13*a**3*tan(c + d*x)*sec(c + d*x)/(8*d) + 4*a**3*tan(c + d*x)/d + 13*a**3*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_31():
    f = (a*cos(c + d*x) + a)**4*cos(c + d*x)**2
    F = 49*a**4*x/16 + 4*a**4*sin(c + d*x)**5/(5*d) - 4*a**4*sin(c + d*x)**3/d + a**4*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 41*a**4*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 49*a**4*sin(c + d*x)*cos(c + d*x)/(16*d) + 8*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_32():
    f = (a*cos(c + d*x) + a)**4
    F = 35*a**4*x/8 - 4*a**4*sin(c + d*x)**3/(3*d) + a**4*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 27*a**4*sin(c + d*x)*cos(c + d*x)/(8*d) + 8*a**4*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_33():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)
    F = 6*a**4*x - a**4*sin(c + d*x)**3/(3*d) + 2*a**4*sin(c + d*x)*cos(c + d*x)/d + 7*a**4*sin(c + d*x)/d + a**4*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_34():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**2
    F = 13*a**4*x/2 + a**4*sin(c + d*x)*cos(c + d*x)/(2*d) + 4*a**4*sin(c + d*x)/d + a**4*tan(c + d*x)/d + 4*a**4*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_35():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**3
    F = 4*a**4*x + a**4*sin(c + d*x)/d + a**4*tan(c + d*x)*sec(c + d*x)/(2*d) + 4*a**4*tan(c + d*x)/d + 13*a**4*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_36():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**4
    F = a**4*x + a**4*tan(c + d*x)**3/(3*d) + 2*a**4*tan(c + d*x)*sec(c + d*x)/d + 7*a**4*tan(c + d*x)/d + 6*a**4*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_37():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**5
    F = 4*a**4*tan(c + d*x)**3/(3*d) + a**4*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 27*a**4*tan(c + d*x)*sec(c + d*x)/(8*d) + 8*a**4*tan(c + d*x)/d + 35*a**4*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_38():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**6
    F = a**4*tan(c + d*x)**5/(5*d) + 8*a**4*tan(c + d*x)**3/(3*d) + a**4*tan(c + d*x)*sec(c + d*x)**3/d + 7*a**4*tan(c + d*x)*sec(c + d*x)/(2*d) + 8*a**4*tan(c + d*x)/d + 7*a**4*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_39():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**7
    F = 4*a**4*tan(c + d*x)**5/(5*d) + 4*a**4*tan(c + d*x)**3/d + a**4*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 41*a**4*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 49*a**4*tan(c + d*x)*sec(c + d*x)/(16*d) + 8*a**4*tan(c + d*x)/d + 49*a**4*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_40():
    f = cos(c + d*x)**5/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**4/(d*(a*cos(c + d*x) + a)) + 15*x/(8*a) + 4*sin(c + d*x)**3/(3*a*d) + 5*sin(c + d*x)*cos(c + d*x)**3/(4*a*d) + 15*sin(c + d*x)*cos(c + d*x)/(8*a*d) - 4*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_41():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**3/(d*(a*cos(c + d*x) + a)) - 3*x/(2*a) - 4*sin(c + d*x)**3/(3*a*d) - 3*sin(c + d*x)*cos(c + d*x)/(2*a*d) + 4*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_42():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**2/(d*(a*cos(c + d*x) + a)) + 3*x/(2*a) + 3*sin(c + d*x)*cos(c + d*x)/(2*a*d) - 2*sin(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_43():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + a)
    F = -x/a + sin(c + d*x)/(a*d) + sin(c + d*x)/(a*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_44():
    f = cos(c + d*x)/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)/(d*(a*cos(c + d*x) + a)) + x/a
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_45():
    f = 1/(a*cos(c + d*x) + a)
    F = sin(c + d*x)/(d*(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_46():
    f = sec(c + d*x)/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)/(d*(a*cos(c + d*x) + a)) + atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_47():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + a)
    F = -tan(c + d*x)/(d*(a*cos(c + d*x) + a)) + 2*tan(c + d*x)/(a*d) - atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_48():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + a)
    F = -tan(c + d*x)*sec(c + d*x)/(d*(a*cos(c + d*x) + a)) + 3*tan(c + d*x)*sec(c + d*x)/(2*a*d) - 2*tan(c + d*x)/(a*d) + 3*atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_49():
    f = sec(c + d*x)**4/(a*cos(c + d*x) + a)
    F = -tan(c + d*x)*sec(c + d*x)**2/(d*(a*cos(c + d*x) + a)) + 4*tan(c + d*x)**3/(3*a*d) - 3*tan(c + d*x)*sec(c + d*x)/(2*a*d) + 4*tan(c + d*x)/(a*d) - 3*atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_50():
    f = cos(c + d*x)**5/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**4/(3*d*(a*cos(c + d*x) + a)**2) - 5*x/a**2 - 4*sin(c + d*x)**3/(a**2*d) - 5*sin(c + d*x)*cos(c + d*x)/(a**2*d) + 12*sin(c + d*x)/(a**2*d) - 10*sin(c + d*x)*cos(c + d*x)**3/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_51():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**3/(3*d*(a*cos(c + d*x) + a)**2) + 7*x/(2*a**2) + 7*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) - 16*sin(c + d*x)/(3*a**2*d) - 8*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_52():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**2/(3*d*(a*cos(c + d*x) + a)**2) - 2*x/a**2 + 4*sin(c + d*x)/(3*a**2*d) + 2*sin(c + d*x)/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_53():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + a)**2
    F = sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) + x/a**2 - 5*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_54():
    f = cos(c + d*x)/(a*cos(c + d*x) + a)**2
    F = 2*sin(c + d*x)/(3*d*(a**2*cos(c + d*x) + a**2)) - sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_55():
    f = (a*cos(c + d*x) + a)**(-2)
    F = sin(c + d*x)/(3*d*(a**2*cos(c + d*x) + a**2)) + sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_56():
    f = sec(c + d*x)/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) + atanh(sin(c + d*x))/(a**2*d) - 4*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_57():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + a)**2
    F = -tan(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) + 10*tan(c + d*x)/(3*a**2*d) - 2*atanh(sin(c + d*x))/(a**2*d) - 2*tan(c + d*x)/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_58():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + a)**2
    F = -tan(c + d*x)*sec(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) + 7*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) - 16*tan(c + d*x)/(3*a**2*d) + 7*atanh(sin(c + d*x))/(2*a**2*d) - 8*tan(c + d*x)*sec(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_59():
    f = sec(c + d*x)**4/(a*cos(c + d*x) + a)**2
    F = -tan(c + d*x)*sec(c + d*x)**2/(3*d*(a*cos(c + d*x) + a)**2) + 4*tan(c + d*x)**3/(a**2*d) - 5*tan(c + d*x)*sec(c + d*x)/(a**2*d) + 12*tan(c + d*x)/(a**2*d) - 5*atanh(sin(c + d*x))/(a**2*d) - 10*tan(c + d*x)*sec(c + d*x)**2/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_60():
    f = cos(c + d*x)**5/(a*cos(c + d*x) + a)**3
    F = -76*sin(c + d*x)*cos(c + d*x)**2/(15*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**4/(5*d*(a*cos(c + d*x) + a)**3) - 11*sin(c + d*x)*cos(c + d*x)**3/(15*a*d*(a*cos(c + d*x) + a)**2) + 13*x/(2*a**3) + 13*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) - 152*sin(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_61():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)**3
    F = 3*sin(c + d*x)/(d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**3/(5*d*(a*cos(c + d*x) + a)**3) - 3*sin(c + d*x)*cos(c + d*x)**2/(5*a*d*(a*cos(c + d*x) + a)**2) - 3*x/a**3 + 9*sin(c + d*x)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_62():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + a)**3
    F = -29*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**2/(5*d*(a*cos(c + d*x) + a)**3) + 7*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2) + x/a**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_63():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + a)**3
    F = 7*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) + sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) - 8*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_64():
    f = cos(c + d*x)/(a*cos(c + d*x) + a)**3
    F = sin(c + d*x)/(5*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) + sin(c + d*x)/(5*a*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_65():
    f = (a*cos(c + d*x) + a)**(-3)
    F = 2*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) + sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) + 2*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_66():
    f = sec(c + d*x)/(a*cos(c + d*x) + a)**3
    F = -22*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) - 7*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2) + atanh(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_67():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + a)**3
    F = -3*tan(c + d*x)/(d*(a**3*cos(c + d*x) + a**3)) - tan(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) - 3*tan(c + d*x)/(5*a*d*(a*cos(c + d*x) + a)**2) + 24*tan(c + d*x)/(5*a**3*d) - 3*atanh(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_68():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + a)**3
    F = -76*tan(c + d*x)*sec(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) - tan(c + d*x)*sec(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) - 11*tan(c + d*x)*sec(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2) + 13*tan(c + d*x)*sec(c + d*x)/(2*a**3*d) - 152*tan(c + d*x)/(15*a**3*d) + 13*atanh(sin(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_69():
    f = cos(c + d*x)**6/(a*cos(c + d*x) + a)**4
    F = -sin(c + d*x)*cos(c + d*x)**5/(7*d*(a*cos(c + d*x) + a)**4) - 2*sin(c + d*x)*cos(c + d*x)**4/(5*a*d*(a*cos(c + d*x) + a)**3) + 21*x/(2*a**4) + 21*sin(c + d*x)*cos(c + d*x)/(2*a**4*d) - 576*sin(c + d*x)/(35*a**4*d) - 288*sin(c + d*x)*cos(c + d*x)**2/(35*a**4*d*(cos(c + d*x) + 1)) - 43*sin(c + d*x)*cos(c + d*x)**3/(35*a**4*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_70():
    f = cos(c + d*x)**5/(a*cos(c + d*x) + a)**4
    F = -sin(c + d*x)*cos(c + d*x)**4/(7*d*(a*cos(c + d*x) + a)**4) - 12*sin(c + d*x)*cos(c + d*x)**3/(35*a*d*(a*cos(c + d*x) + a)**3) - 4*x/a**4 + 244*sin(c + d*x)/(105*a**4*d) + 4*sin(c + d*x)/(a**4*d*(cos(c + d*x) + 1)) - 88*sin(c + d*x)*cos(c + d*x)**2/(105*a**4*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_71():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)**4
    F = -sin(c + d*x)*cos(c + d*x)**3/(7*d*(a*cos(c + d*x) + a)**4) - 2*sin(c + d*x)*cos(c + d*x)**2/(7*a*d*(a*cos(c + d*x) + a)**3) + x/a**4 - 43*sin(c + d*x)/(21*a**4*d*(cos(c + d*x) + 1)) + 11*sin(c + d*x)/(21*a**4*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_72():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + a)**4
    F = -sin(c + d*x)*cos(c + d*x)**2/(7*d*(a*cos(c + d*x) + a)**4) + 8*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3) + 12*sin(c + d*x)/(35*a**4*d*(cos(c + d*x) + 1)) - 18*sin(c + d*x)/(35*a**4*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_73():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + a)**4
    F = 13*sin(c + d*x)/(105*d*(a**4*cos(c + d*x) + a**4)) + 13*sin(c + d*x)/(105*d*(a**2*cos(c + d*x) + a**2)**2) + sin(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) - 11*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_74():
    f = cos(c + d*x)/(a*cos(c + d*x) + a)**4
    F = 8*sin(c + d*x)/(105*d*(a**4*cos(c + d*x) + a**4)) + 8*sin(c + d*x)/(105*d*(a**2*cos(c + d*x) + a**2)**2) - sin(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) + 4*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_75():
    f = (a*cos(c + d*x) + a)**(-4)
    F = 2*sin(c + d*x)/(35*d*(a**4*cos(c + d*x) + a**4)) + 2*sin(c + d*x)/(35*d*(a**2*cos(c + d*x) + a**2)**2) + sin(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) + 3*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_76():
    f = sec(c + d*x)/(a*cos(c + d*x) + a)**4
    F = -sin(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) - 2*sin(c + d*x)/(7*a*d*(a*cos(c + d*x) + a)**3) + atanh(sin(c + d*x))/(a**4*d) - 32*sin(c + d*x)/(21*a**4*d*(cos(c + d*x) + 1)) - 11*sin(c + d*x)/(21*a**4*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_77():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + a)**4
    F = -tan(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) - 12*tan(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3) + 664*tan(c + d*x)/(105*a**4*d) - 4*atanh(sin(c + d*x))/(a**4*d) - 4*tan(c + d*x)/(a**4*d*(cos(c + d*x) + 1)) - 88*tan(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_78():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + a)**4
    F = -tan(c + d*x)*sec(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) - 2*tan(c + d*x)*sec(c + d*x)/(5*a*d*(a*cos(c + d*x) + a)**3) + 21*tan(c + d*x)*sec(c + d*x)/(2*a**4*d) - 576*tan(c + d*x)/(35*a**4*d) + 21*atanh(sin(c + d*x))/(2*a**4*d) - 288*tan(c + d*x)*sec(c + d*x)/(35*a**4*d*(cos(c + d*x) + 1)) - 43*tan(c + d*x)*sec(c + d*x)/(35*a**4*d*(cos(c + d*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_79():
    f = cos(c + d*x)**7/(a*cos(c + d*x) + a)**5
    F = -3832*sin(c + d*x)*cos(c + d*x)**2/(315*d*(a**5*cos(c + d*x) + a**5)) - sin(c + d*x)*cos(c + d*x)**6/(9*d*(a*cos(c + d*x) + a)**5) - 17*sin(c + d*x)*cos(c + d*x)**5/(63*a*d*(a*cos(c + d*x) + a)**4) - 28*sin(c + d*x)*cos(c + d*x)**4/(45*a**2*d*(a*cos(c + d*x) + a)**3) - 577*sin(c + d*x)*cos(c + d*x)**3/(315*a**3*d*(a*cos(c + d*x) + a)**2) + 31*x/(2*a**5) + 31*sin(c + d*x)*cos(c + d*x)/(2*a**5*d) - 7664*sin(c + d*x)/(315*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_80():
    f = cos(c + d*x)**6/(a*cos(c + d*x) + a)**5
    F = 5*sin(c + d*x)/(d*(a**5*cos(c + d*x) + a**5)) - sin(c + d*x)*cos(c + d*x)**5/(9*d*(a*cos(c + d*x) + a)**5) - 5*sin(c + d*x)*cos(c + d*x)**4/(21*a*d*(a*cos(c + d*x) + a)**4) - 29*sin(c + d*x)*cos(c + d*x)**3/(63*a**2*d*(a*cos(c + d*x) + a)**3) - 67*sin(c + d*x)*cos(c + d*x)**2/(63*a**3*d*(a*cos(c + d*x) + a)**2) - 5*x/a**5 + 181*sin(c + d*x)/(63*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_81():
    f = cos(c + d*x)**5/(a*cos(c + d*x) + a)**5
    F = -661*sin(c + d*x)/(315*d*(a**5*cos(c + d*x) + a**5)) - sin(c + d*x)*cos(c + d*x)**4/(9*d*(a*cos(c + d*x) + a)**5) - 13*sin(c + d*x)*cos(c + d*x)**3/(63*a*d*(a*cos(c + d*x) + a)**4) - 34*sin(c + d*x)*cos(c + d*x)**2/(105*a**2*d*(a*cos(c + d*x) + a)**3) + 173*sin(c + d*x)/(315*a**3*d*(a*cos(c + d*x) + a)**2) + x/a**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_82():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)**5
    F = 83*sin(c + d*x)/(315*d*(a**5*cos(c + d*x) + a**5)) - sin(c + d*x)*cos(c + d*x)**3/(9*d*(a*cos(c + d*x) + a)**5) - 11*sin(c + d*x)*cos(c + d*x)**2/(63*a*d*(a*cos(c + d*x) + a)**4) + 67*sin(c + d*x)/(315*a**2*d*(a*cos(c + d*x) + a)**3) - 142*sin(c + d*x)/(315*a**3*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_83():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + a)**5
    F = 5*sin(c + d*x)/(63*d*(a**5*cos(c + d*x) + a**5)) - sin(c + d*x)*cos(c + d*x)**2/(9*d*(a*cos(c + d*x) + a)**5) + sin(c + d*x)/(7*a*d*(a*cos(c + d*x) + a)**4) - 17*sin(c + d*x)/(63*a**2*d*(a*cos(c + d*x) + a)**3) + 5*sin(c + d*x)/(63*a**3*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_84():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + a)**5
    F = 2*sin(c + d*x)/(45*d*(a**5*cos(c + d*x) + a**5)) + sin(c + d*x)/(9*d*(a*cos(c + d*x) + a)**5) - 2*sin(c + d*x)/(9*a*d*(a*cos(c + d*x) + a)**4) + sin(c + d*x)/(15*a**2*d*(a*cos(c + d*x) + a)**3) + 2*sin(c + d*x)/(45*a**3*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_85():
    f = cos(c + d*x)/(a*cos(c + d*x) + a)**5
    F = 2*sin(c + d*x)/(63*d*(a**5*cos(c + d*x) + a**5)) - sin(c + d*x)/(9*d*(a*cos(c + d*x) + a)**5) + 2*sin(c + d*x)/(63*a*d*(a**2*cos(c + d*x) + a**2)**2) + 5*sin(c + d*x)/(63*a*d*(a*cos(c + d*x) + a)**4) + sin(c + d*x)/(21*a**2*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_86():
    f = (a*cos(c + d*x) + a)**(-5)
    F = 8*sin(c + d*x)/(315*d*(a**5*cos(c + d*x) + a**5)) + sin(c + d*x)/(9*d*(a*cos(c + d*x) + a)**5) + 8*sin(c + d*x)/(315*a*d*(a**2*cos(c + d*x) + a**2)**2) + 4*sin(c + d*x)/(63*a*d*(a*cos(c + d*x) + a)**4) + 4*sin(c + d*x)/(105*a**2*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_87():
    f = sec(c + d*x)/(a*cos(c + d*x) + a)**5
    F = -488*sin(c + d*x)/(315*d*(a**5*cos(c + d*x) + a**5)) - sin(c + d*x)/(9*d*(a*cos(c + d*x) + a)**5) - 13*sin(c + d*x)/(63*a*d*(a*cos(c + d*x) + a)**4) - 34*sin(c + d*x)/(105*a**2*d*(a*cos(c + d*x) + a)**3) - 173*sin(c + d*x)/(315*a**3*d*(a*cos(c + d*x) + a)**2) + atanh(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_88():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + a)**5
    F = -5*tan(c + d*x)/(d*(a**5*cos(c + d*x) + a**5)) - tan(c + d*x)/(9*d*(a*cos(c + d*x) + a)**5) - 5*tan(c + d*x)/(21*a*d*(a*cos(c + d*x) + a)**4) - 29*tan(c + d*x)/(63*a**2*d*(a*cos(c + d*x) + a)**3) - 67*tan(c + d*x)/(63*a**3*d*(a*cos(c + d*x) + a)**2) + 496*tan(c + d*x)/(63*a**5*d) - 5*atanh(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_89():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + a)**5
    F = -3832*tan(c + d*x)*sec(c + d*x)/(315*d*(a**5*cos(c + d*x) + a**5)) - tan(c + d*x)*sec(c + d*x)/(9*d*(a*cos(c + d*x) + a)**5) - 17*tan(c + d*x)*sec(c + d*x)/(63*a*d*(a*cos(c + d*x) + a)**4) - 28*tan(c + d*x)*sec(c + d*x)/(45*a**2*d*(a*cos(c + d*x) + a)**3) - 577*tan(c + d*x)*sec(c + d*x)/(315*a**3*d*(a*cos(c + d*x) + a)**2) + 31*tan(c + d*x)*sec(c + d*x)/(2*a**5*d) - 7664*tan(c + d*x)/(315*a**5*d) + 31*atanh(sin(c + d*x))/(2*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_90():
    f = cos(c + d*x)**5/(a*cos(c + d*x) + a)**6
    F = -sin(c + d*x)*cos(c + d*x)**4/(11*d*(a*cos(c + d*x) + a)**6) - 14*sin(c + d*x)*cos(c + d*x)**3/(99*a*d*(a*cos(c + d*x) + a)**5) - 118*sin(c + d*x)*cos(c + d*x)**2/(693*a**2*d*(a*cos(c + d*x) + a)**4) + 146*sin(c + d*x)/(693*a**6*d*(cos(c + d*x) + 1)) - 268*sin(c + d*x)/(693*a**6*d*(cos(c + d*x) + 1)**2) + 130*sin(c + d*x)/(693*a**6*d*(cos(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_91():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)**6
    F = -sin(c + d*x)*cos(c + d*x)**3/(11*d*(a*cos(c + d*x) + a)**6) - 4*sin(c + d*x)*cos(c + d*x)**2/(33*a*d*(a*cos(c + d*x) + a)**5) + 9*sin(c + d*x)/(77*a**2*d*(a*cos(c + d*x) + a)**4) + 61*sin(c + d*x)/(1155*a**6*d*(cos(c + d*x) + 1)) + 61*sin(c + d*x)/(1155*a**6*d*(cos(c + d*x) + 1)**2) - 241*sin(c + d*x)/(1155*a**6*d*(cos(c + d*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_92():
    f = sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**4
    F = 2*a*sin(c + d*x)*cos(c + d*x)**4/(9*d*sqrt(a*cos(c + d*x) + a)) + 16*a*sin(c + d*x)*cos(c + d*x)**3/(63*d*sqrt(a*cos(c + d*x) + a)) + 32*a*sin(c + d*x)/(45*d*sqrt(a*cos(c + d*x) + a)) - 64*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(315*d) + 32*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(105*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_93():
    f = sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**3
    F = 2*a*sin(c + d*x)*cos(c + d*x)**3/(7*d*sqrt(a*cos(c + d*x) + a)) + 4*a*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)) - 8*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(35*d) + 12*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_94():
    f = sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**2
    F = 14*a*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) - 4*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*d) + 2*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_95():
    f = sqrt(a*cos(c + d*x) + a)*cos(c + d*x)
    F = 2*a*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)) + 2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_96():
    f = sqrt(a*cos(c + d*x) + a)
    F = 2*a*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_97():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)
    F = 2*sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_98():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**2
    F = sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + a*tan(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_99():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**3
    F = 3*sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a*tan(c + d*x)*sec(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)) + 3*a*tan(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_100():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**4
    F = 5*sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a*tan(c + d*x)*sec(c + d*x)**2/(3*d*sqrt(a*cos(c + d*x) + a)) + 5*a*tan(c + d*x)*sec(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)) + 5*a*tan(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_101():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**4/(9*d*sqrt(a*cos(c + d*x) + a)) + 34*a**2*sin(c + d*x)*cos(c + d*x)**3/(63*d*sqrt(a*cos(c + d*x) + a)) + 68*a**2*sin(c + d*x)/(45*d*sqrt(a*cos(c + d*x) + a)) - 136*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(315*d) + 68*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_102():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2
    F = 152*a**2*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)) + 38*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(105*d) - 4*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(35*d) + 2*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_103():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)
    F = 8*a**2*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)) + 2*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(5*d) + 2*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_104():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = 8*a**2*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)) + 2*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_105():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)
    F = 2*a**(sympy.S(3)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 2*a**2*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_106():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2
    F = 3*a**(sympy.S(3)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + a**2*tan(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_107():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3
    F = 7*a**(sympy.S(3)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a**2*tan(c + d*x)*sec(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)) + 7*a**2*tan(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_108():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4
    F = 11*a**(sympy.S(3)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a**2*tan(c + d*x)*sec(c + d*x)**2/(3*d*sqrt(a*cos(c + d*x) + a)) + 11*a**2*tan(c + d*x)*sec(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)) + 11*a**2*tan(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_109():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3
    F = 46*a**3*sin(c + d*x)*cos(c + d*x)**4/(99*d*sqrt(a*cos(c + d*x) + a)) + 710*a**3*sin(c + d*x)*cos(c + d*x)**3/(693*d*sqrt(a*cos(c + d*x) + a)) + 284*a**3*sin(c + d*x)/(99*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**4/(11*d) - 568*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(693*d) + 284*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(231*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_110():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2
    F = 832*a**3*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)) + 208*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(315*d) + 26*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(105*d) - 4*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(63*d) + 2*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sin(c + d*x)/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_111():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)
    F = 64*a**3*sin(c + d*x)/(21*d*sqrt(a*cos(c + d*x) + a)) + 16*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(21*d) + 2*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(7*d) + 2*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_112():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = 64*a**3*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) + 16*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*d) + 2*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_113():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)
    F = 2*a**(sympy.S(5)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 14*a**3*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_114():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**2
    F = 5*a**(sympy.S(5)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + a**3*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)) + a**2*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_115():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**3
    F = 19*a**(sympy.S(5)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + 9*a**3*tan(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)) + a**2*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_116():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4
    F = 25*a**(sympy.S(5)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + 13*a**3*tan(c + d*x)*sec(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)) + 25*a**3*tan(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a)) + a**2*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**2/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_117():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**5
    F = 163*a**(sympy.S(5)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(64*d) + 17*a**3*tan(c + d*x)*sec(c + d*x)**2/(24*d*sqrt(a*cos(c + d*x) + a)) + 163*a**3*tan(c + d*x)*sec(c + d*x)/(96*d*sqrt(a*cos(c + d*x) + a)) + 163*a**3*tan(c + d*x)/(64*d*sqrt(a*cos(c + d*x) + a)) + a**2*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_118():
    f = (a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = 256*a**4*sin(c + d*x)/(35*d*sqrt(a*cos(c + d*x) + a)) + 64*a**3*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(35*d) + 24*a**2*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(35*d) + 2*a*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_119():
    f = cos(c + d*x)**4/sqrt(a*cos(c + d*x) + a)
    F = 2*sin(c + d*x)*cos(c + d*x)**3/(7*d*sqrt(a*cos(c + d*x) + a)) - 2*sin(c + d*x)*cos(c + d*x)**2/(35*d*sqrt(a*cos(c + d*x) + a)) - 148*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)) + 62*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(105*a*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_120():
    f = cos(c + d*x)**3/sqrt(a*cos(c + d*x) + a)
    F = 2*sin(c + d*x)*cos(c + d*x)**2/(5*d*sqrt(a*cos(c + d*x) + a)) + 28*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) - 2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*a*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_121():
    f = cos(c + d*x)**2/sqrt(a*cos(c + d*x) + a)
    F = -4*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)) + 2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*a*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_122():
    f = cos(c + d*x)/sqrt(a*cos(c + d*x) + a)
    F = 2*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_123():
    f = 1/sqrt(a*cos(c + d*x) + a)
    F = sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_124():
    f = sec(c + d*x)/sqrt(a*cos(c + d*x) + a)
    F = 2*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_125():
    f = sec(c + d*x)**2/sqrt(a*cos(c + d*x) + a)
    F = tan(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)) - atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_126():
    f = sec(c + d*x)**3/sqrt(a*cos(c + d*x) + a)
    F = tan(c + d*x)*sec(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)) - tan(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)) + 7*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_127():
    f = sec(c + d*x)**4/sqrt(a*cos(c + d*x) + a)
    F = tan(c + d*x)*sec(c + d*x)**2/(3*d*sqrt(a*cos(c + d*x) + a)) - tan(c + d*x)*sec(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)) + 7*tan(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a)) - 9*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*sqrt(a)*d) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_128():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*cos(c + d*x)**3/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 9*sin(c + d*x)*cos(c + d*x)**2/(10*a*d*sqrt(a*cos(c + d*x) + a)) + 31*sin(c + d*x)/(5*a*d*sqrt(a*cos(c + d*x) + a)) - 13*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(10*a**2*d) - 15*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_129():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*cos(c + d*x)**2/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - 13*sin(c + d*x)/(3*a*d*sqrt(a*cos(c + d*x) + a)) + 7*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(6*a**2*d) + 11*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_130():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 2*sin(c + d*x)/(a*d*sqrt(a*cos(c + d*x) + a)) - 7*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_131():
    f = cos(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_132():
    f = (a*cos(c + d*x) + a)**(sympy.S(-3)/2)
    F = sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_133():
    f = sec(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 2*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_134():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 3*tan(c + d*x)/(2*a*d*sqrt(a*cos(c + d*x) + a)) - 3*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 9*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_135():
    f = sec(c + d*x)**3/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -tan(c + d*x)*sec(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + tan(c + d*x)*sec(c + d*x)/(a*d*sqrt(a*cos(c + d*x) + a)) - 7*tan(c + d*x)/(4*a*d*sqrt(a*cos(c + d*x) + a)) + 19*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*a**(sympy.S(3)/2)*d) - 13*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_136():
    f = cos(c + d*x)**4/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*cos(c + d*x)**3/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 17*sin(c + d*x)*cos(c + d*x)**2/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - 197*sin(c + d*x)/(24*a**2*d*sqrt(a*cos(c + d*x) + a)) + 95*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(48*a**3*d) + 163*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_137():
    f = cos(c + d*x)**3/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*cos(c + d*x)**2/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 13*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 9*sin(c + d*x)/(4*a**2*d*sqrt(a*cos(c + d*x) + a)) - 75*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_138():
    f = cos(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 13*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 19*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_139():
    f = cos(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 5*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 5*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_140():
    f = (a*cos(c + d*x) + a)**(sympy.S(-5)/2)
    F = sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 3*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_141():
    f = sec(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 11*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 2*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 43*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_142():
    f = sec(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -tan(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 15*tan(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 35*tan(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) - 5*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 115*sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_143():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 10*a*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 10*a*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_144():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_145():
    f = (a*cos(c + d*x) + a)*sqrt(cos(c + d*x))
    F = 2*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_146():
    f = (a*cos(c + d*x) + a)/sqrt(cos(c + d*x))
    F = 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_147():
    f = (a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_148():
    f = (a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) - 2*a*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_149():
    f = (a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(7)/2)
    F = 6*a*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 2*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_150():
    f = (a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 4*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 32*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 20*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 32*a**2*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 20*a**2*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_151():
    f = (a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 4*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 8*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(7*d) + 12*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 8*a**2*elliptic_f(c/2 + d*x/2, 2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_152():
    f = (a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))
    F = 2*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 4*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 16*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_153():
    f = (a*cos(c + d*x) + a)**2/sqrt(cos(c + d*x))
    F = 2*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 4*a**2*elliptic_e(c/2 + d*x/2, 2)/d + 8*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_154():
    f = (a*cos(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 4*a**2*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_155():
    f = (a*cos(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(5)/2)
    F = 4*a**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) - 4*a**2*elliptic_e(c/2 + d*x/2, 2)/d + 8*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_156():
    f = (a*cos(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(7)/2)
    F = 16*a**2*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 4*a**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 16*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_157():
    f = (a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 6*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 68*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 44*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 68*a**3*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 44*a**3*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_158():
    f = (a*cos(c + d*x) + a)**3*sqrt(cos(c + d*x))
    F = 2*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 6*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 52*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 28*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 52*a**3*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_159():
    f = (a*cos(c + d*x) + a)**3/sqrt(cos(c + d*x))
    F = 2*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/d + 36*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**3*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_160():
    f = (a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a**3*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 4*a**3*elliptic_e(c/2 + d*x/2, 2)/d + 20*a**3*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_161():
    f = (a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(5)/2)
    F = 6*a**3*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a**3*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) - 4*a**3*elliptic_e(c/2 + d*x/2, 2)/d + 20*a**3*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_162():
    f = (a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(7)/2)
    F = 36*a**3*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 2*a**3*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**3*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 36*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**3*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_163():
    f = (a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(9)/2)
    F = 28*a**3*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 52*a**3*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 6*a**3*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**3*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) - 28*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 52*a**3*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_164():
    f = (a*cos(c + d*x) + a)**4*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)/(11*d) + 8*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 150*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(77*d) + 128*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 904*a**4*sin(c + d*x)*sqrt(cos(c + d*x))/(231*d) + 128*a**4*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 904*a**4*elliptic_f(c/2 + d*x/2, 2)/(231*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_165():
    f = (a*cos(c + d*x) + a)**4*sqrt(cos(c + d*x))
    F = 2*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 8*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 122*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 32*a**4*sin(c + d*x)*sqrt(cos(c + d*x))/(7*d) + 152*a**4*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 32*a**4*elliptic_f(c/2 + d*x/2, 2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_166():
    f = (a*cos(c + d*x) + a)**4/sqrt(cos(c + d*x))
    F = 2*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 8*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 94*a**4*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 64*a**4*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 136*a**4*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_167():
    f = (a*cos(c + d*x) + a)**4/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**4*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 8*a**4*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a**4*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 56*a**4*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 32*a**4*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_168():
    f = (a*cos(c + d*x) + a)**4/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**4*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 8*a**4*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*a**4*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 40*a**4*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_169():
    f = (a*cos(c + d*x) + a)**4/cos(c + d*x)**(sympy.S(7)/2)
    F = 66*a**4*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 8*a**4*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**4*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 56*a**4*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 32*a**4*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_170():
    f = (a*cos(c + d*x) + a)**4/cos(c + d*x)**(sympy.S(9)/2)
    F = 64*a**4*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 94*a**4*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 8*a**4*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**4*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) - 64*a**4*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 136*a**4*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_171():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(d*(a*cos(c + d*x) + a)) + 7*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a*d) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) + 21*elliptic_e(c/2 + d*x/2, 2)/(5*a*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_172():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(d*(a*cos(c + d*x) + a)) + 5*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) - 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) + 5*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_173():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*cos(c + d*x) + a)) + 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) - elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_174():
    f = sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*cos(c + d*x) + a)) - elliptic_e(c/2 + d*x/2, 2)/(a*d) + elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_175():
    f = 1/((a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*cos(c + d*x) + a)) + elliptic_e(c/2 + d*x/2, 2)/(a*d) + elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_176():
    f = 1/((a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(d*(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 3*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) - 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) - elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_177():
    f = 1/((a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(d*(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - 3*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) + 5*sin(c + d*x)/(3*a*d*cos(c + d*x)**(sympy.S(3)/2)) + 3*elliptic_e(c/2 + d*x/2, 2)/(a*d) + 5*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_178():
    f = cos(c + d*x)**(sympy.S(9)/2)/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(3*d*(a*cos(c + d*x) + a)**2) + 56*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a**2*d) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d) + 56*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(a**2*d) - 3*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_179():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*d*(a*cos(c + d*x) + a)**2) + 10*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d) - 7*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 10*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 7*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_180():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d*(a*cos(c + d*x) + a)**2) + 4*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_181():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*cos(c + d*x) + a)**2) - elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 2*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) + sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_182():
    f = sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**2
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*cos(c + d*x) + a)**2) + elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_183():
    f = 1/((a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*cos(c + d*x) + a)**2) + elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 2*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_184():
    f = 1/((a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))) + 4*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) - 4*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - 5*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 5*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_185():
    f = 1/((a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) - 7*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) + 10*sin(c + d*x)/(3*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + 7*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + 10*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - 7*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_186():
    f = cos(c + d*x)**(sympy.S(11)/2)/(a*cos(c + d*x) + a)**3
    F = -63*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(10*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)/(5*d*(a*cos(c + d*x) + a)**3) - 4*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(5*a*d*(a*cos(c + d*x) + a)**2) + 77*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(10*a**3*d) - 21*sin(c + d*x)*sqrt(cos(c + d*x))/(2*a**3*d) + 231*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) - 21*elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_187():
    f = cos(c + d*x)**(sympy.S(9)/2)/(a*cos(c + d*x) + a)**3
    F = -119*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(30*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(5*d*(a*cos(c + d*x) + a)**3) - 2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*a*d*(a*cos(c + d*x) + a)**2) + 11*sin(c + d*x)*sqrt(cos(c + d*x))/(2*a**3*d) - 119*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + 11*elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_188():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**3
    F = -13*sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(5*d*(a*cos(c + d*x) + a)**3) - 8*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a*d*(a*cos(c + d*x) + a)**2) + 49*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) - 13*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_189():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**3
    F = 9*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*(a*cos(c + d*x) + a)**3) - 2*sin(c + d*x)*sqrt(cos(c + d*x))/(5*a*d*(a*cos(c + d*x) + a)**2) - 9*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_190():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**3
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*cos(c + d*x) + a)**3) + 4*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) - elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_191():
    f = sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**3
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) + sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*cos(c + d*x) + a)**3) + sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) + elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_192():
    f = 1/((a*cos(c + d*x) + a)**3*sqrt(cos(c + d*x)))
    F = -9*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) - sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*cos(c + d*x) + a)**3) - 2*sin(c + d*x)*sqrt(cos(c + d*x))/(5*a*d*(a*cos(c + d*x) + a)**2) + 9*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_193():
    f = 1/((a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2))
    F = -13*sin(c + d*x)/(6*d*(a**3*cos(c + d*x) + a**3)*sqrt(cos(c + d*x))) - sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3*sqrt(cos(c + d*x))) - 8*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))) + 49*sin(c + d*x)/(10*a**3*d*sqrt(cos(c + d*x))) - 49*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) - 13*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_194():
    f = 1/((a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(5)/2))
    F = -119*sin(c + d*x)/(30*d*(a**3*cos(c + d*x) + a**3)*cos(c + d*x)**(sympy.S(3)/2)) - sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)) - 2*sin(c + d*x)/(3*a*d*(a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) - 119*sin(c + d*x)/(10*a**3*d*sqrt(cos(c + d*x))) + 11*sin(c + d*x)/(2*a**3*d*cos(c + d*x)**(sympy.S(3)/2)) + 119*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + 11*elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_195():
    f = sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = 5*sqrt(a)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) + 5*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(12*d*sqrt(a*cos(c + d*x) + a)) + 5*a*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_196():
    f = sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = 3*sqrt(a)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(a*cos(c + d*x) + a)) + 3*a*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_197():
    f = sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))
    F = sqrt(a)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + a*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_198():
    f = sqrt(a*cos(c + d*x) + a)/sqrt(cos(c + d*x))
    F = 2*sqrt(a)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_199():
    f = sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_200():
    f = sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = 4*a*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_201():
    f = sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(7)/2)
    F = 16*a*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 8*a*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_202():
    f = sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(9)/2)
    F = 32*a*sin(c + d*x)/(35*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 16*a*sin(c + d*x)/(35*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 12*a*sin(c + d*x)/(35*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 2*a*sin(c + d*x)/(7*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_203():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 11*a**(sympy.S(3)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) + 11*a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(12*d*sqrt(a*cos(c + d*x) + a)) + 11*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_204():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = 7*a**(sympy.S(3)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(a*cos(c + d*x) + a)) + 7*a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_205():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = 3*a**(sympy.S(3)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + a**2*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_206():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(3)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 2*a**2*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_207():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 10*a**2*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_208():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 12*a**2*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 6*a**2*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_209():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 208*a**2*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 104*a**2*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 26*a**2*sin(c + d*x)/(35*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**2*sin(c + d*x)/(7*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_210():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = 163*a**(sympy.S(5)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(64*d) + 17*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(24*d*sqrt(a*cos(c + d*x) + a)) + 163*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(96*d*sqrt(a*cos(c + d*x) + a)) + 163*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(64*d*sqrt(a*cos(c + d*x) + a)) + a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_211():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = 25*a**(sympy.S(5)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + 13*a**3*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(12*d*sqrt(a*cos(c + d*x) + a)) + 25*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d*sqrt(a*cos(c + d*x) + a)) + a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_212():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = 19*a**(sympy.S(5)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + 9*a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a)) + a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_213():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = 5*a**(sympy.S(5)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d - a**3*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_214():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(5)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 14*a**3*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_215():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 86*a**3*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 22*a**3*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_216():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 92*a**3*sin(c + d*x)/(21*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 46*a**3*sin(c + d*x)/(21*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 6*a**3*sin(c + d*x)/(7*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_217():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = 1168*a**3*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 584*a**3*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 146*a**3*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 38*a**3*sin(c + d*x)/(63*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_218():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/4)
    F = 4*a**2*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(1)/4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_219():
    f = sqrt(a*cos(e + f*x) + a)/sqrt(cos(e + f*x))
    F = 2*sqrt(a)*asin(sqrt(a)*sin(e + f*x)/sqrt(a*cos(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_220():
    f = sqrt(-a*cos(e + f*x) + a)/sqrt(-cos(e + f*x))
    F = -2*sqrt(a)*asin(sqrt(a)*sin(e + f*x)/sqrt(-a*cos(e + f*x) + a))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_221():
    f = cos(c + d*x)**(sympy.S(5)/2)/sqrt(a*cos(c + d*x) + a)
    F = sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(a*cos(c + d*x) + a)) - sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a)) + 7*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*sqrt(a)*d) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_222():
    f = cos(c + d*x)**(sympy.S(3)/2)/sqrt(a*cos(c + d*x) + a)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(a*cos(c + d*x) + a)) - asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_223():
    f = sqrt(cos(c + d*x))/sqrt(a*cos(c + d*x) + a)
    F = 2*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_224():
    f = 1/(sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_225():
    f = 1/(sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_226():
    f = 1/(sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = -2*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_227():
    f = 1/(sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    F = 26*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - 2*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_228():
    f = cos(c + d*x)**(sympy.S(5)/2)/sqrt(cos(c + d*x) + 1)
    F = -sqrt(2)*asin(sin(c + d*x)/(cos(c + d*x) + 1))/d + 7*asin(sin(c + d*x)/sqrt(cos(c + d*x) + 1))/(4*d) + sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(cos(c + d*x) + 1)) - sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_229():
    f = cos(c + d*x)**(sympy.S(3)/2)/sqrt(cos(c + d*x) + 1)
    F = sqrt(2)*asin(sin(c + d*x)/(cos(c + d*x) + 1))/d - asin(sin(c + d*x)/sqrt(cos(c + d*x) + 1))/d + sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_230():
    f = sqrt(cos(c + d*x))/sqrt(cos(c + d*x) + 1)
    F = -sqrt(2)*asin(sin(c + d*x)/(cos(c + d*x) + 1))/d + 2*asin(sin(c + d*x)/sqrt(cos(c + d*x) + 1))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_231():
    f = 1/(sqrt(cos(c + d*x) + 1)*sqrt(cos(c + d*x)))
    F = sqrt(2)*asin(sin(c + d*x)/(cos(c + d*x) + 1))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_232():
    f = 1/(sqrt(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*asin(sin(c + d*x)/(cos(c + d*x) + 1))/d + 2*sin(c + d*x)/(d*sqrt(cos(c + d*x) + 1)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_233():
    f = 1/(sqrt(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(5)/2))
    F = sqrt(2)*asin(sin(c + d*x)/(cos(c + d*x) + 1))/d - 2*sin(c + d*x)/(3*d*sqrt(cos(c + d*x) + 1)*sqrt(cos(c + d*x))) + 2*sin(c + d*x)/(3*d*sqrt(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_234():
    f = 1/(sqrt(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(7)/2))
    F = -sqrt(2)*asin(sin(c + d*x)/(cos(c + d*x) + 1))/d + 26*sin(c + d*x)/(15*d*sqrt(cos(c + d*x) + 1)*sqrt(cos(c + d*x))) - 2*sin(c + d*x)/(15*d*sqrt(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(3)/2)) + 2*sin(c + d*x)/(5*d*sqrt(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_235():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 3*sin(c + d*x)*sqrt(cos(c + d*x))/(2*a*d*sqrt(a*cos(c + d*x) + a)) - 3*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + 9*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_236():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 2*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_237():
    f = sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_238():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_239():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + 5*sin(c + d*x)/(2*a*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - 7*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_240():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) - 19*sin(c + d*x)/(6*a*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 7*sin(c + d*x)/(6*a*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 11*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_241():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 15*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 35*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) - 5*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + 115*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_242():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 11*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 2*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - 43*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_243():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 7*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_244():
    f = sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 5*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_245():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 9*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 19*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_246():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))) - 13*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + 49*sin(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - 75*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_247():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)) - 17*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) - 299*sin(c + d*x)/(48*a**2*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 95*sin(c + d*x)/(48*a**2*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 163*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_248():
    f = cos(c + d*x)**(sympy.S(9)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - 7*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 259*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 189*sin(c + d*x)*sqrt(cos(c + d*x))/(64*a**3*d*sqrt(a*cos(c + d*x) + a)) - 7*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(7)/2)*d) + 637*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_249():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - 17*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 49*sin(c + d*x)*sqrt(cos(c + d*x))/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 2*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(7)/2)*d) - 177*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_250():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - 13*sin(c + d*x)*sqrt(cos(c + d*x))/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 67*sin(c + d*x)*sqrt(cos(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 5*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_251():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) + 3*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 17*sin(c + d*x)*sqrt(cos(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 7*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_252():
    f = sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) + sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 13*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_253():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(cos(c + d*x)))
    F = -sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 103*sin(c + d*x)*sqrt(cos(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 63*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_254():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(cos(c + d*x))) - 19*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))) - 199*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + 691*sin(c + d*x)/(192*a**3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - 363*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_255():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**(sympy.S(3)/2)) - 23*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)) - 109*sin(c + d*x)/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) - 629*sin(c + d*x)/(64*a**3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 193*sin(c + d*x)/(64*a**3*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 1015*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_256():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**(sympy.S(9)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(8*d*(a*cos(c + d*x) + a)**(sympy.S(9)/2)) - 19*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(96*a*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - 187*sin(c + d*x)*sqrt(cos(c + d*x))/(768*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 853*sin(c + d*x)*sqrt(cos(c + d*x))/(3072*a**3*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 35*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(2048*a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_257():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(9)/2)
    F = -sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(8*d*(a*cos(c + d*x) + a)**(sympy.S(9)/2)) - 5*sin(c + d*x)*sqrt(cos(c + d*x))/(32*a*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) + 33*sin(c + d*x)*sqrt(cos(c + d*x))/(256*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + 73*sin(c + d*x)*sqrt(cos(c + d*x))/(1024*a**3*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 45*sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(2048*a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_258():
    f = 1/(sqrt(cos(x) + 1)*sqrt(cos(x)))
    F = sqrt(2)*asin(sin(x)/(cos(x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_259():
    f = 1/(sqrt(a*cos(x) + a)*sqrt(cos(x)))
    F = sqrt(2)*atan(sqrt(2)*sqrt(a)*sin(x)/(2*sqrt(a*cos(x) + a)*sqrt(cos(x))))/sqrt(a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_260():
    f = sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = -3*sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/(sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*d) - a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(-a*cos(c + d*x) + a)) + 3*a*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(-a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_261():
    f = sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))
    F = sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/(sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/d - a*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(-a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_262():
    f = sqrt(-a*cos(c + d*x) + a)/sqrt(cos(c + d*x))
    F = -2*sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/(sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_263():
    f = sqrt(-a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_264():
    f = sqrt(-a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = -4*a*sin(c + d*x)/(3*d*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a*sin(c + d*x)/(3*d*sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_265():
    f = sqrt(-a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(7)/2)
    F = 16*a*sin(c + d*x)/(15*d*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - 8*a*sin(c + d*x)/(15*d*sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(5*d*sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_266():
    f = sqrt(1 - cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = -3*atanh(sin(c + d*x)/(sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/(4*d) - sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(1 - cos(c + d*x))) + 3*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_267():
    f = sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))
    F = atanh(sin(c + d*x)/(sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d - sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_268():
    f = sqrt(1 - cos(c + d*x))/sqrt(cos(c + d*x))
    F = -2*atanh(sin(c + d*x)/(sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_269():
    f = sqrt(1 - cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*sin(c + d*x)/(d*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_270():
    f = sqrt(1 - cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = -4*sin(c + d*x)/(3*d*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))) + 2*sin(c + d*x)/(3*d*sqrt(1 - cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_271():
    f = sqrt(1 - cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = 16*sin(c + d*x)/(15*d*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))) - 8*sin(c + d*x)/(15*d*sqrt(1 - cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + 2*sin(c + d*x)/(5*d*sqrt(1 - cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_272():
    f = cos(c + d*x)**(sympy.S(5)/2)/sqrt(-a*cos(c + d*x) + a)
    F = sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(-a*cos(c + d*x) + a)) + sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(-a*cos(c + d*x) + a)) + 7*atanh(sqrt(a)*sin(c + d*x)/(sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_273():
    f = cos(c + d*x)**(sympy.S(3)/2)/sqrt(-a*cos(c + d*x) + a)
    F = sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(-a*cos(c + d*x) + a)) + atanh(sqrt(a)*sin(c + d*x)/(sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_274():
    f = sqrt(cos(c + d*x))/sqrt(-a*cos(c + d*x) + a)
    F = 2*atanh(sqrt(a)*sin(c + d*x)/(sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_275():
    f = 1/(sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = -sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_276():
    f = 1/(sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*sin(c + d*x)/(d*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_277():
    f = 1/(sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*sin(c + d*x)/(3*d*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*sin(c + d*x)/(3*d*sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_278():
    f = 1/(sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    F = 26*sin(c + d*x)/(15*d*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*sin(c + d*x)/(15*d*sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*sin(c + d*x)/(5*d*sqrt(-a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) - sqrt(2)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(-a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_279():
    f = cos(c + d*x)**(sympy.S(5)/2)/sqrt(1 - cos(c + d*x))
    F = 7*atanh(sin(c + d*x)/(sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/(4*d) - sqrt(2)*atanh(sqrt(2)*sin(c + d*x)/(2*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d + sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(1 - cos(c + d*x))) + sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_280():
    f = cos(c + d*x)**(sympy.S(3)/2)/sqrt(1 - cos(c + d*x))
    F = atanh(sin(c + d*x)/(sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d - sqrt(2)*atanh(sqrt(2)*sin(c + d*x)/(2*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d + sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(1 - cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_281():
    f = sqrt(cos(c + d*x))/sqrt(1 - cos(c + d*x))
    F = 2*atanh(sin(c + d*x)/(sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d - sqrt(2)*atanh(sqrt(2)*sin(c + d*x)/(2*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_282():
    f = 1/(sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x)))
    F = -sqrt(2)*atanh(sqrt(2)*sin(c + d*x)/(2*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_283():
    f = 1/(sqrt(1 - cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(2)*atanh(sqrt(2)*sin(c + d*x)/(2*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d + 2*sin(c + d*x)/(d*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_284():
    f = 1/(sqrt(1 - cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = -sqrt(2)*atanh(sqrt(2)*sin(c + d*x)/(2*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))))/d + 2*sin(c + d*x)/(3*d*sqrt(1 - cos(c + d*x))*sqrt(cos(c + d*x))) + 2*sin(c + d*x)/(3*d*sqrt(1 - cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_285():
    f = (a*cos(c + d*x) + a)**(sympy.S(1)/3)*cos(c + d*x)**(sympy.S(4)/3)
    F = 2**(sympy.S(5)/6)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S(1)/6, sympy.S(3)/2, 1 - cos(c + d*x), sympy.S.Half - cos(c + d*x)/2)/(d*(cos(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_286():
    f = (a*cos(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)**(sympy.S(4)/3)
    F = 2*2**(sympy.S(1)/6)*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S(-1)/6, sympy.S(3)/2, 1 - cos(c + d*x), sympy.S.Half - cos(c + d*x)/2)/(d*(cos(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_287():
    f = (a*cos(c + d*x) + a)**(sympy.S(2)/3)*cos(c + d*x)**(sympy.S(5)/3)
    F = 2*2**(sympy.S(1)/6)*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/3, sympy.S(-1)/6, sympy.S(3)/2, 1 - cos(c + d*x), sympy.S.Half - cos(c + d*x)/2)/(d*(cos(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_288():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 6*a*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 6*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_289():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_290():
    f = (a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_291():
    f = (a*cos(c + d*x) + a)*sqrt(sec(c + d*x))
    F = 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_292():
    f = (a*cos(c + d*x) + a)/sqrt(sec(c + d*x))
    F = 2*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_293():
    f = (a*cos(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_294():
    f = (a*cos(c + d*x) + a)/sec(c + d*x)**(sympy.S(5)/2)
    F = 10*a*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 6*a*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 10*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_295():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 4*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 16*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 16*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_296():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 4*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/d - 4*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 8*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_297():
    f = (a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/d + 4*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_298():
    f = (a*cos(c + d*x) + a)**2*sqrt(sec(c + d*x))
    F = 2*a**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 8*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_299():
    f = (a*cos(c + d*x) + a)**2/sqrt(sec(c + d*x))
    F = 4*a**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 16*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_300():
    f = (a*cos(c + d*x) + a)**2/sec(c + d*x)**(sympy.S(3)/2)
    F = 8*a**2*sin(c + d*x)/(7*d*sqrt(sec(c + d*x))) + 4*a**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 12*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 8*a**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_301():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 6*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 52*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 28*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 28*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 52*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_302():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/d + 36*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 36*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_303():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 6*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/d - 4*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 20*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_304():
    f = (a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/d + 2*a**3*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 20*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_305():
    f = (a*cos(c + d*x) + a)**3*sqrt(sec(c + d*x))
    F = 2*a**3*sin(c + d*x)/(d*sqrt(sec(c + d*x))) + 2*a**3*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 36*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_306():
    f = (a*cos(c + d*x) + a)**3/sqrt(sec(c + d*x))
    F = 52*a**3*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 6*a**3*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**3*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 28*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 52*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_307():
    f = (a*cos(c + d*x) + a)**3/sec(c + d*x)**(sympy.S(3)/2)
    F = 44*a**3*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 68*a**3*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a**3*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**3*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 68*a**3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 44*a**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_308():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 8*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 94*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 64*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 64*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 136*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_309():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 8*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 66*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 56*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 32*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_310():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**4*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 8*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/d + 2*a**4*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 40*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_311():
    f = (a*cos(c + d*x) + a)**4*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**4*sin(c + d*x)*sqrt(sec(c + d*x))/d + 8*a**4*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a**4*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 56*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 32*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_312():
    f = (a*cos(c + d*x) + a)**4*sqrt(sec(c + d*x))
    F = 94*a**4*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 8*a**4*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**4*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 64*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 136*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_313():
    f = (a*cos(c + d*x) + a)**4/sqrt(sec(c + d*x))
    F = 32*a**4*sin(c + d*x)/(7*d*sqrt(sec(c + d*x))) + 122*a**4*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 8*a**4*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a**4*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 152*a**4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 32*a**4*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_314():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(d*(a*sec(c + d*x) + a)) + 5*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d) - 3*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d) + 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_315():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*(a*sec(c + d*x) + a)) + 3*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d) - 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) - sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_316():
    f = sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) + sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_317():
    f = 1/((a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) - sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_318():
    f = 1/((a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) + 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) - sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_319():
    f = 1/((a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) + 5*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) - 3*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_320():
    f = 1/((a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - 5*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) + 7*sin(c + d*x)/(5*a*d*sec(c + d*x)**(sympy.S(3)/2)) + 21*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_321():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(3*d*(a*sec(c + d*x) + a)**2) + 10*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d) - 7*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d) + 7*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 7*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_322():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*d*(a*sec(c + d*x) + a)**2) + 4*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d) - 4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 5*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_323():
    f = sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) + sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_324():
    f = 1/((a*cos(c + d*x) + a)**2*sqrt(sec(c + d*x)))
    F = sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_325():
    f = 1/((a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) - sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) + sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_326():
    f = 1/((a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*(a*sec(c + d*x) + a)**2) + 4*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 5*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_327():
    f = 1/((a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) + 10*sin(c + d*x)/(3*a**2*d*sqrt(sec(c + d*x))) - 7*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + 10*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - 7*sin(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_328():
    f = 1/((a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(9)/2))
    F = -sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)) - 5*sin(c + d*x)/(a**2*d*sqrt(sec(c + d*x))) + 56*sin(c + d*x)/(15*a**2*d*sec(c + d*x)**(sympy.S(3)/2)) + 56*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a**2*d) - 5*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - 3*sin(c + d*x)/(a**2*d*(sec(c + d*x) + 1)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_329():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**3
    F = -13*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 8*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(15*a*d*(a*sec(c + d*x) + a)**2) + 49*sin(c + d*x)*sqrt(sec(c + d*x))/(10*a**3*d) - 49*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) - 13*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_330():
    f = sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**3
    F = -9*sin(c + d*x)*sqrt(sec(c + d*x))/(10*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*a*d*(a*sec(c + d*x) + a)**2) + 9*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_331():
    f = 1/((a*cos(c + d*x) + a)**3*sqrt(sec(c + d*x)))
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) - 4*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_332():
    f = 1/((a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) + sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) - sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) - sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_333():
    f = 1/((a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(5)/2))
    F = sin(c + d*x)*sqrt(sec(c + d*x))/(2*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) + 2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*a*d*(a*sec(c + d*x) + a)**2) - 9*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_334():
    f = 1/((a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(7)/2))
    F = -13*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) - sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*(a*sec(c + d*x) + a)**3) - 8*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + 49*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) - 13*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_335():
    f = 1/((a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(9)/2))
    F = -119*sin(c + d*x)/(30*d*(a**3*sec(c + d*x) + a**3)*sqrt(sec(c + d*x))) - sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x))) - 2*sin(c + d*x)/(3*a*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) + 11*sin(c + d*x)/(2*a**3*d*sqrt(sec(c + d*x))) - 119*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + 11*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_336():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d*sqrt(a*cos(c + d*x) + a)) + 12*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d*sqrt(a*cos(c + d*x) + a)) + 16*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(35*d*sqrt(a*cos(c + d*x) + a)) + 32*a*sin(c + d*x)*sqrt(sec(c + d*x))/(35*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_337():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*sqrt(a*cos(c + d*x) + a)) + 8*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(a*cos(c + d*x) + a)) + 16*a*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_338():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) + 4*a*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_339():
    f = sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_340():
    f = sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))
    F = 2*sqrt(a)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_341():
    f = sqrt(a*cos(c + d*x) + a)/sqrt(sec(c + d*x))
    F = sqrt(a)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d + a*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_342():
    f = sqrt(a*cos(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = 3*sqrt(a)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + 3*a*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a*sin(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_343():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d*sqrt(a*cos(c + d*x) + a)) + 26*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d*sqrt(a*cos(c + d*x) + a)) + 104*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) + 208*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_344():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*sqrt(a*cos(c + d*x) + a)) + 6*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*sqrt(a*cos(c + d*x) + a)) + 12*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_345():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) + 10*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_346():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d + 2*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_347():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = 3*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d + a**2*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_348():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = 7*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + 7*a**2*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a**2*sin(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_349():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = 11*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + 11*a**2*sin(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + 11*a**2*sin(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + a**2*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_350():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(11)/2)
    F = 38*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d*sqrt(a*cos(c + d*x) + a)) + 146*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) + 584*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*d*sqrt(a*cos(c + d*x) + a)) + 1168*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_351():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 6*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(7*d*sqrt(a*cos(c + d*x) + a)) + 46*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d*sqrt(a*cos(c + d*x) + a)) + 92*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(21*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_352():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = 22*a**3*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(a*cos(c + d*x) + a)) + 86*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_353():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d + 14*a**3*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_354():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = 5*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d - a**3*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_355():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = 19*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + 9*a**3*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_356():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = 25*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + 25*a**3*sin(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + 13*a**3*sin(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_357():
    f = (a*cos(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = 163*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(64*d) + 163*a**3*sin(c + d*x)/(64*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + 163*a**3*sin(c + d*x)/(96*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + 17*a**3*sin(c + d*x)/(24*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + a**2*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(4*d*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_358():
    f = sec(c + d*x)**(sympy.S(7)/2)/sqrt(cos(c + d*x) + 1)
    F = -sqrt(2)*sqrt(cos(c + d*x))*asin(sin(c + d*x)/(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d + 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*sqrt(cos(c + d*x) + 1)) - 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(cos(c + d*x) + 1)) + 26*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_359():
    f = sec(c + d*x)**(sympy.S(5)/2)/sqrt(cos(c + d*x) + 1)
    F = sqrt(2)*sqrt(cos(c + d*x))*asin(sin(c + d*x)/(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d + 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*sqrt(cos(c + d*x) + 1)) - 2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_360():
    f = sec(c + d*x)**(sympy.S(3)/2)/sqrt(cos(c + d*x) + 1)
    F = -sqrt(2)*sqrt(cos(c + d*x))*asin(sin(c + d*x)/(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d + 2*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_361():
    f = sqrt(sec(c + d*x))/sqrt(cos(c + d*x) + 1)
    F = sqrt(2)*sqrt(cos(c + d*x))*asin(sin(c + d*x)/(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_362():
    f = 1/(sqrt(cos(c + d*x) + 1)*sqrt(sec(c + d*x)))
    F = -sqrt(2)*sqrt(cos(c + d*x))*asin(sin(c + d*x)/(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d + 2*sqrt(cos(c + d*x))*asin(sin(c + d*x)/sqrt(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_363():
    f = 1/(sqrt(cos(c + d*x) + 1)*sec(c + d*x)**(sympy.S(3)/2))
    F = sqrt(2)*sqrt(cos(c + d*x))*asin(sin(c + d*x)/(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d - sqrt(cos(c + d*x))*asin(sin(c + d*x)/sqrt(cos(c + d*x) + 1))*sqrt(sec(c + d*x))/d + sin(c + d*x)/(d*sqrt(cos(c + d*x) + 1)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_364():
    f = sec(c + d*x)**(sympy.S(7)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*sqrt(a*cos(c + d*x) + a)) - 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(a*cos(c + d*x) + a)) + 26*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_365():
    f = sec(c + d*x)**(sympy.S(5)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) - 2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_366():
    f = sec(c + d*x)**(sympy.S(3)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_367():
    f = 1/(sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) - sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d) + sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_368():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 7*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*a*d*sqrt(a*cos(c + d*x) + a)) - 19*sin(c + d*x)*sqrt(sec(c + d*x))/(6*a*d*sqrt(a*cos(c + d*x) + a)) + 11*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_369():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 5*sin(c + d*x)*sqrt(sec(c + d*x))/(2*a*d*sqrt(a*cos(c + d*x) + a)) - 7*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_370():
    f = sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 3*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_371():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_372():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 2*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) - 5*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_373():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) + 3*sin(c + d*x)/(2*a*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) - 3*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) + 9*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_374():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 17*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 95*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(48*a**2*d*sqrt(a*cos(c + d*x) + a)) - 299*sin(c + d*x)*sqrt(sec(c + d*x))/(48*a**2*d*sqrt(a*cos(c + d*x) + a)) + 163*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_375():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 13*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 49*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) - 75*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_376():
    f = sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - 9*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 19*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_377():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 5*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_378():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + 7*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 3*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_379():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)) - 11*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 2*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) - 43*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_380():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)) - 15*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) + 35*sin(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) - 5*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) + 115*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_381():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - 23*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 109*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 193*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(64*a**3*d*sqrt(a*cos(c + d*x) + a)) - 629*sin(c + d*x)*sqrt(sec(c + d*x))/(64*a**3*d*sqrt(a*cos(c + d*x) + a)) + 1015*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_382():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - 19*sin(c + d*x)*sqrt(sec(c + d*x))/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - 199*sin(c + d*x)*sqrt(sec(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + 691*sin(c + d*x)*sqrt(sec(c + d*x))/(192*a**3*d*sqrt(a*cos(c + d*x) + a)) - 363*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_383():
    f = sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x))) - 5*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - 103*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 63*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_384():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x)))
    F = sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x))) + sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - 5*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 13*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_385():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x))) + 3*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + 17*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 7*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_386():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(3)/2)) - 13*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + 67*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 5*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_387():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(5)/2)) - 17*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)) - 49*sin(c + d*x)/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 2*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(7)/2)*d) - 177*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_388():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(9)/2))
    F = -sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(7)/2)) - 7*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)) - 259*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) + 189*sin(c + d*x)/(64*a**3*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) - 7*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(7)/2)*d) + 637*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_389():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(9)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = -sin(c + d*x)/(8*d*(a*cos(c + d*x) + a)**(sympy.S(9)/2)*sec(c + d*x)**(sympy.S(3)/2)) - 5*sin(c + d*x)/(32*a*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x))) + 33*sin(c + d*x)/(256*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + 73*sin(c + d*x)/(1024*a**3*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 45*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(2048*a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_390():
    f = 1/((a*cos(c + d*x) + a)**(sympy.S(9)/2)*sec(c + d*x)**(sympy.S(7)/2))
    F = -sin(c + d*x)/(8*d*(a*cos(c + d*x) + a)**(sympy.S(9)/2)*sec(c + d*x)**(sympy.S(5)/2)) - 19*sin(c + d*x)/(96*a*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(3)/2)) - 187*sin(c + d*x)/(768*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + 853*sin(c + d*x)/(3072*a**3*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 35*sqrt(2)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(2048*a**(sympy.S(9)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_391():
    f = (a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/4)
    F = 4*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(1)/4)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_392():
    f = (a*cos(c + d*x) + a)**4*cos(c + d*x)**m
    F = -4*a**4*(2*m + 5)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*(m + 3)*sqrt(sin(c + d*x)**2)) + a**4*(4*m**2 + 29*m + 55)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 2)*(m + 3)*(m + 4)) - a**4*(8*m**2 + 40*m + 35)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*(m + 2)*(m + 4)*sqrt(sin(c + d*x)**2)) + (a**2*cos(c + d*x) + a**2)**2*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 4)) + (2*m + 10)*(a**4*cos(c + d*x) + a**4)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 3)*(m + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_393():
    f = (a*cos(c + d*x) + a)**3*cos(c + d*x)**m
    F = a**3*(2*m + 7)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 2)*(m + 3)) - a**3*(4*m + 11)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*(m + 3)*sqrt(sin(c + d*x)**2)) - a**3*(4*m + 5)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*(m + 2)*sqrt(sin(c + d*x)**2)) + (a**3*cos(c + d*x) + a**3)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_394():
    f = (a*cos(c + d*x) + a)**2*cos(c + d*x)**m
    F = a**2*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 2)) - 2*a**2*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*sqrt(sin(c + d*x)**2)) - a**2*(2*m + 3)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*(m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_395():
    f = (a*cos(c + d*x) + a)*cos(c + d*x)**m
    F = -a*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*sqrt(sin(c + d*x)**2)) - a*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_396():
    f = cos(c + d*x)**m/(a*cos(c + d*x) + a)
    F = sin(c + d*x)*cos(c + d*x)**m/(d*(a*cos(c + d*x) + a)) + m*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(a*d*(m + 1)*sqrt(sin(c + d*x)**2)) - sin(c + d*x)*cos(c + d*x)**m*hyper((sympy.S.Half, m/2), (m/2 + 1,), cos(c + d*x)**2)/(a*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_397():
    f = cos(c + d*x)**m/(a*cos(c + d*x) + a)**2
    F = -sin(c + d*x)*cos(c + d*x)**(m + 1)/(3*d*(a*cos(c + d*x) + a)**2) + m*(1 - 2*m)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(3*a**2*d*(m + 1)*sqrt(sin(c + d*x)**2)) - (2 - 2*m)*(m + 1)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(3*a**2*d*(m + 2)*sqrt(sin(c + d*x)**2)) - (2 - 2*m)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_398():
    f = (a + b*cos(c + d*x))*cos(c + d*x)**7
    F = -a*sin(c + d*x)**7/(7*d) + 3*a*sin(c + d*x)**5/(5*d) - a*sin(c + d*x)**3/d + a*sin(c + d*x)/d + 35*b*x/128 + b*sin(c + d*x)*cos(c + d*x)**7/(8*d) + 7*b*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 35*b*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 35*b*sin(c + d*x)*cos(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_399():
    f = (a + b*cos(c + d*x))*cos(c + d*x)**6
    F = 5*a*x/16 + a*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a*sin(c + d*x)*cos(c + d*x)/(16*d) - b*sin(c + d*x)**7/(7*d) + 3*b*sin(c + d*x)**5/(5*d) - b*sin(c + d*x)**3/d + b*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_400():
    f = (a + b*cos(c + d*x))*cos(c + d*x)**5
    F = a*sin(c + d*x)**5/(5*d) - 2*a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d + 5*b*x/16 + b*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*b*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*b*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_401():
    f = (a + b*cos(c + d*x))*cos(c + d*x)**4
    F = 3*a*x/8 + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) + b*sin(c + d*x)**5/(5*d) - 2*b*sin(c + d*x)**3/(3*d) + b*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_402():
    f = (a + b*cos(c + d*x))*cos(c + d*x)**3
    F = -a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d + 3*b*x/8 + b*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*b*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_403():
    f = (a + b*cos(c + d*x))*cos(c + d*x)**2
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) - b*sin(c + d*x)**3/(3*d) + b*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_404():
    f = (a + b*cos(c + d*x))*cos(c + d*x)
    F = a*sin(c + d*x)/d + b*x/2 + b*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_405():
    f = a + b*cos(c + d*x)
    F = a*x + b*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_406():
    f = (a + b*cos(c + d*x))*sec(c + d*x)
    F = a*atanh(sin(c + d*x))/d + b*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_407():
    f = (a + b*cos(c + d*x))*sec(c + d*x)**2
    F = a*tan(c + d*x)/d + b*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_408():
    f = (a + b*cos(c + d*x))*sec(c + d*x)**3
    F = a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*atanh(sin(c + d*x))/(2*d) + b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_409():
    f = (a + b*cos(c + d*x))*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*tan(c + d*x)*sec(c + d*x)/(2*d) + b*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_410():
    f = (a + b*cos(c + d*x))*sec(c + d*x)**5
    F = a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a*atanh(sin(c + d*x))/(8*d) + b*tan(c + d*x)**3/(3*d) + b*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_411():
    f = (a + b*cos(c + d*x))*sec(c + d*x)**6
    F = a*tan(c + d*x)**5/(5*d) + 2*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*b*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*b*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_412():
    f = (a + b*cos(c + d*x))**2*cos(c + d*x)**4
    F = 2*a*b*sin(c + d*x)**5/(5*d) - 4*a*b*sin(c + d*x)**3/(3*d) + 2*a*b*sin(c + d*x)/d + b**2*sin(c + d*x)*cos(c + d*x)**5/(6*d) + x*(3*a**2/8 + 5*b**2/16) + (6*a**2 + 5*b**2)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + (6*a**2 + 5*b**2)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_413():
    f = (a + b*cos(c + d*x))**2*cos(c + d*x)**3
    F = 3*a*b*x/4 + a*b*sin(c + d*x)*cos(c + d*x)**3/(2*d) + 3*a*b*sin(c + d*x)*cos(c + d*x)/(4*d) + b**2*sin(c + d*x)**5/(5*d) + (a**2 + b**2)*sin(c + d*x)/d - (a**2 + 2*b**2)*sin(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_414():
    f = (a + b*cos(c + d*x))**2*cos(c + d*x)**2
    F = -2*a*b*sin(c + d*x)**3/(3*d) + 2*a*b*sin(c + d*x)/d + b**2*sin(c + d*x)*cos(c + d*x)**3/(4*d) + x*(a**2/2 + 3*b**2/8) + (4*a**2 + 3*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_415():
    f = (a + b*cos(c + d*x))**2*cos(c + d*x)
    F = a*b*x + a*b*sin(c + d*x)*cos(c + d*x)/(3*d) + (a + b*cos(c + d*x))**2*sin(c + d*x)/(3*d) + (2*a**2 + 2*b**2)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_416():
    f = (a + b*cos(c + d*x))**2
    F = 2*a*b*sin(c + d*x)/d + b**2*sin(c + d*x)*cos(c + d*x)/(2*d) + x*(a**2 + b**2/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_417():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)
    F = a**2*atanh(sin(c + d*x))/d + 2*a*b*x + b**2*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_418():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**2
    F = a**2*tan(c + d*x)/d + 2*a*b*atanh(sin(c + d*x))/d + b**2*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_419():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**3
    F = a**2*tan(c + d*x)*sec(c + d*x)/(2*d) + 2*a*b*tan(c + d*x)/d + (a**2 + 2*b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_420():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**4
    F = a**2*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a*b*tan(c + d*x)*sec(c + d*x)/d + a*b*atanh(sin(c + d*x))/d + (2*a**2 + 3*b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_421():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**5
    F = a**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 2*a*b*tan(c + d*x)**3/(3*d) + 2*a*b*tan(c + d*x)/d + (3*a**2 + 4*b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + (3*a**2 + 4*b**2)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_422():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**6
    F = a**2*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a*b*tan(c + d*x)*sec(c + d*x)**3/(2*d) + 3*a*b*tan(c + d*x)*sec(c + d*x)/(4*d) + 3*a*b*atanh(sin(c + d*x))/(4*d) + (4*a**2 + 5*b**2)*tan(c + d*x)**3/(15*d) + (4*a**2 + 5*b**2)*tan(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_423():
    f = (a + b*cos(c + d*x))**3*cos(c + d*x)**2
    F = a*x*(4*a**2 + 9*b**2)/8 - a*(6*a**2 - 71*b**2)*sin(c + d*x)*cos(c + d*x)/(120*d) - a*(a + b*cos(c + d*x))**3*sin(c + d*x)/(20*b*d) + (a + b*cos(c + d*x))**4*sin(c + d*x)/(5*b*d) - (a + b*cos(c + d*x))**2*(3*a**2 - 16*b**2)*sin(c + d*x)/(60*b*d) - (3*a**4 - 52*a**2*b**2 - 16*b**4)*sin(c + d*x)/(30*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_424():
    f = (a + b*cos(c + d*x))**3*cos(c + d*x)
    F = a*(a + b*cos(c + d*x))**2*sin(c + d*x)/(4*d) + a*(a**2 + 4*b**2)*sin(c + d*x)/(2*d) + 3*b*x*(4*a**2 + b**2)/8 + b*(2*a**2 + 3*b**2)*sin(c + d*x)*cos(c + d*x)/(8*d) + (a + b*cos(c + d*x))**3*sin(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_425():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)
    F = a**3*atanh(sin(c + d*x))/d + 5*a*b**2*sin(c + d*x)/(2*d) + b**2*(a + b*cos(c + d*x))*sin(c + d*x)/(2*d) + b*x*(6*a**2 + b**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_426():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**2
    F = 3*a**2*b*atanh(sin(c + d*x))/d + a**2*(a + b*cos(c + d*x))*tan(c + d*x)/d + 3*a*b**2*x - b*(a**2 - b**2)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_427():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**3
    F = 5*a**2*b*tan(c + d*x)/(2*d) + a**2*(a + b*cos(c + d*x))*tan(c + d*x)*sec(c + d*x)/(2*d) + a*(a**2 + 6*b**2)*atanh(sin(c + d*x))/(2*d) + b**3*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_428():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**4
    F = 7*a**2*b*tan(c + d*x)*sec(c + d*x)/(6*d) + a**2*(a + b*cos(c + d*x))*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a*(2*a**2 + 9*b**2)*tan(c + d*x)/(3*d) + b*(3*a**2 + 2*b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_429():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**5
    F = 3*a**2*b*tan(c + d*x)*sec(c + d*x)**2/(4*d) + a**2*(a + b*cos(c + d*x))*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*(a**2 + 4*b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a*(a**2 + 4*b**2)*atanh(sin(c + d*x))/(8*d) + b*(2*a**2 + b**2)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_430():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**6
    F = 11*a**2*b*tan(c + d*x)*sec(c + d*x)**3/(20*d) + a**2*(a + b*cos(c + d*x))*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a*(4*a**2 + 15*b**2)*tan(c + d*x)**3/(15*d) + a*(4*a**2 + 15*b**2)*tan(c + d*x)/(5*d) + b*(9*a**2 + 4*b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + b*(9*a**2 + 4*b**2)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_431():
    f = (a + b*cos(c + d*x))**4*cos(c + d*x)**3
    F = 8*a*b**3*sin(c + d*x)*cos(c + d*x)**5/(21*d) + a*b*x*(6*a**2 + 5*b**2)/4 + a*b*(6*a**2 + 5*b**2)*sin(c + d*x)*cos(c + d*x)**3/(6*d) + a*b*(6*a**2 + 5*b**2)*sin(c + d*x)*cos(c + d*x)/(4*d) + b**2*(a + b*cos(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**4/(7*d) + b**2*(37*a**2 + 6*b**2)*sin(c + d*x)*cos(c + d*x)**4/(35*d) - (35*a**4 + 168*a**2*b**2 + 24*b**4)*sin(c + d*x)**3/(105*d) + (35*a**4 + 168*a**2*b**2 + 24*b**4)*sin(c + d*x)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_432():
    f = (a + b*cos(c + d*x))**4*cos(c + d*x)**2
    F = -a*(a + b*cos(c + d*x))**4*sin(c + d*x)/(30*b*d) - a*(a + b*cos(c + d*x))**2*(4*a**2 - 53*b**2)*sin(c + d*x)/(120*b*d) - a*(4*a**4 - 121*a**2*b**2 - 128*b**4)*sin(c + d*x)/(60*b*d) + x*(a**4/2 + 9*a**2*b**2/4 + 5*b**4/16) - (8*a**4 - 178*a**2*b**2 - 75*b**4)*sin(c + d*x)*cos(c + d*x)/(240*d) + (a + b*cos(c + d*x))**5*sin(c + d*x)/(6*b*d) - (a + b*cos(c + d*x))**3*(4*a**2 - 25*b**2)*sin(c + d*x)/(120*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_433():
    f = (a + b*cos(c + d*x))**4*cos(c + d*x)
    F = a*b*x*(4*a**2 + 3*b**2)/2 + a*b*(6*a**2 + 29*b**2)*sin(c + d*x)*cos(c + d*x)/(30*d) + a*(a + b*cos(c + d*x))**3*sin(c + d*x)/(5*d) + (a + b*cos(c + d*x))**4*sin(c + d*x)/(5*d) + (a + b*cos(c + d*x))**2*(3*a**2 + 4*b**2)*sin(c + d*x)/(15*d) + (6*a**4 + 56*a**2*b**2 + 8*b**4)*sin(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_434():
    f = (a + b*cos(c + d*x))**4
    F = 7*a*b*(a + b*cos(c + d*x))**2*sin(c + d*x)/(12*d) + a*b*(19*a**2 + 16*b**2)*sin(c + d*x)/(6*d) + b**2*(26*a**2 + 9*b**2)*sin(c + d*x)*cos(c + d*x)/(24*d) + b*(a + b*cos(c + d*x))**3*sin(c + d*x)/(4*d) + x*(a**4 + 3*a**2*b**2 + 3*b**4/8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_435():
    f = (a + b*cos(c + d*x))**4*sec(c + d*x)
    F = a**4*atanh(sin(c + d*x))/d + 4*a*b**3*sin(c + d*x)*cos(c + d*x)/(3*d) + 2*a*b*x*(2*a**2 + b**2) + b**2*(a + b*cos(c + d*x))**2*sin(c + d*x)/(3*d) + b**2*(17*a**2 + 2*b**2)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_436():
    f = (a + b*cos(c + d*x))**4*sec(c + d*x)**2
    F = 4*a**3*b*atanh(sin(c + d*x))/d + a**2*(a + b*cos(c + d*x))**2*tan(c + d*x)/d - 2*a*b*(a**2 - 2*b**2)*sin(c + d*x)/d + b**2*x*(12*a**2 + b**2)/2 - b**2*(2*a**2 - b**2)*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_437():
    f = (a + b*cos(c + d*x))**4*sec(c + d*x)**3
    F = 3*a**3*b*tan(c + d*x)/d + a**2*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)/(2*d) + a**2*(a**2 + 12*b**2)*atanh(sin(c + d*x))/(2*d) + 4*a*b**3*x - b**2*(a**2 - 2*b**2)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_438():
    f = (a + b*cos(c + d*x))**4*sec(c + d*x)**4
    F = 4*a**3*b*tan(c + d*x)*sec(c + d*x)/(3*d) + a**2*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a**2*(2*a**2 + 17*b**2)*tan(c + d*x)/(3*d) + 2*a*b*(a**2 + 2*b**2)*atanh(sin(c + d*x))/d + b**4*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_439():
    f = (a + b*cos(c + d*x))**4*sec(c + d*x)**5
    F = 5*a**3*b*tan(c + d*x)*sec(c + d*x)**2/(6*d) + a**2*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a**2*(3*a**2 + 22*b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + 4*a*b*(2*a**2 + 3*b**2)*tan(c + d*x)/(3*d) + (3*a**4 + 24*a**2*b**2 + 8*b**4)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_440():
    f = (a + b*cos(c + d*x))**4*sec(c + d*x)**6
    F = 3*a**3*b*tan(c + d*x)*sec(c + d*x)**3/(5*d) + a**2*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a**2*(4*a**2 + 27*b**2)*tan(c + d*x)*sec(c + d*x)**2/(15*d) + a*b*(3*a**2 + 4*b**2)*tan(c + d*x)*sec(c + d*x)/(2*d) + a*b*(3*a**2 + 4*b**2)*atanh(sin(c + d*x))/(2*d) + (8*a**4 + 60*a**2*b**2 + 15*b**4)*tan(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_441():
    f = (a + b*cos(c + d*x))**4*sec(c + d*x)**7
    F = 7*a**3*b*tan(c + d*x)*sec(c + d*x)**4/(15*d) + a**2*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)**5/(6*d) + a**2*(5*a**2 + 32*b**2)*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 4*a*b*(4*a**2 + 5*b**2)*tan(c + d*x)**3/(15*d) + 4*a*b*(4*a**2 + 5*b**2)*tan(c + d*x)/(5*d) + (5*a**4 + 36*a**2*b**2 + 8*b**4)*tan(c + d*x)*sec(c + d*x)/(16*d) + (5*a**4 + 36*a**2*b**2 + 8*b**4)*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_442():
    f = cos(c + d*x)**5/(a + b*cos(c + d*x))
    F = -2*a**5*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*sqrt(a - b)*sqrt(a + b)) - a*sin(c + d*x)*cos(c + d*x)**2/(3*b**2*d) - a*(3*a**2 + 2*b**2)*sin(c + d*x)/(3*b**4*d) + sin(c + d*x)*cos(c + d*x)**3/(4*b*d) + (4*a**2 + 3*b**2)*sin(c + d*x)*cos(c + d*x)/(8*b**3*d) + x*(8*a**4 + 4*a**2*b**2 + 3*b**4)/(8*b**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_443():
    f = cos(c + d*x)**4/(a + b*cos(c + d*x))
    F = 2*a**4*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*sqrt(a - b)*sqrt(a + b)) - a*sin(c + d*x)*cos(c + d*x)/(2*b**2*d) - a*x*(2*a**2 + b**2)/(2*b**4) + sin(c + d*x)*cos(c + d*x)**2/(3*b*d) + (3*a**2 + 2*b**2)*sin(c + d*x)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_444():
    f = cos(c + d*x)**3/(a + b*cos(c + d*x))
    F = -2*a**3*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*sqrt(a - b)*sqrt(a + b)) - a*sin(c + d*x)/(b**2*d) + sin(c + d*x)*cos(c + d*x)/(2*b*d) + x*(2*a**2 + b**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_445():
    f = cos(c + d*x)**2/(a + b*cos(c + d*x))
    F = 2*a**2*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*sqrt(a - b)*sqrt(a + b)) - a*x/b**2 + sin(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_446():
    f = cos(c + d*x)/(a + b*cos(c + d*x))
    F = -2*a*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b*d*sqrt(a - b)*sqrt(a + b)) + x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_447():
    f = 1/(a + b*cos(c + d*x))
    F = 2*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_448():
    f = sec(c + d*x)/(a + b*cos(c + d*x))
    F = -2*b*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*sqrt(a - b)*sqrt(a + b)) + atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_449():
    f = sec(c + d*x)**2/(a + b*cos(c + d*x))
    F = tan(c + d*x)/(a*d) + 2*b**2*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*sqrt(a - b)*sqrt(a + b)) - b*atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_450():
    f = sec(c + d*x)**3/(a + b*cos(c + d*x))
    F = tan(c + d*x)*sec(c + d*x)/(2*a*d) - b*tan(c + d*x)/(a**2*d) - 2*b**3*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*sqrt(a - b)*sqrt(a + b)) + (a**2 + 2*b**2)*atanh(sin(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_451():
    f = sec(c + d*x)**4/(a + b*cos(c + d*x))
    F = tan(c + d*x)*sec(c + d*x)**2/(3*a*d) - b*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) + (2*a**2 + 3*b**2)*tan(c + d*x)/(3*a**3*d) + 2*b**4*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*sqrt(a - b)*sqrt(a + b)) - b*(a**2 + 2*b**2)*atanh(sin(c + d*x))/(2*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_452():
    f = cos(c + d*x)**5/(a + b*cos(c + d*x))**2
    F = 2*a**4*(4*a**2 - 5*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - a**2*sin(c + d*x)*cos(c + d*x)**3/(b*d*(a + b*cos(c + d*x))*(a**2 - b**2)) - a*(2*a**2 - b**2)*sin(c + d*x)*cos(c + d*x)/(b**3*d*(a**2 - b**2)) - a*x*(4*a**2 + b**2)/b**5 + (4*a**2 - b**2)*sin(c + d*x)*cos(c + d*x)**2/(3*b**2*d*(a**2 - b**2)) + (12*a**4 - 7*a**2*b**2 - 2*b**4)*sin(c + d*x)/(3*b**4*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_453():
    f = cos(c + d*x)**3/(a + b*cos(c + d*x))**2
    F = -a**2*sin(c + d*x)*cos(c + d*x)/(b*d*(a + b*cos(c + d*x))*(a**2 - b**2)) + 2*a**2*(2*a**2 - 3*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - 2*a*x/b**3 + (2*a**2 - b**2)*sin(c + d*x)/(b**2*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_454():
    f = cos(c + d*x)**2/(a + b*cos(c + d*x))**2
    F = -a**2*sin(c + d*x)/(b*d*(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*a*(a**2 - 2*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_455():
    f = cos(c + d*x)/(a + b*cos(c + d*x))**2
    F = a*sin(c + d*x)/(d*(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*b*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_456():
    f = (a + b*cos(c + d*x))**(-2)
    F = 2*a*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - b*sin(c + d*x)/(d*(a + b*cos(c + d*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_457():
    f = sec(c + d*x)/(a + b*cos(c + d*x))**2
    F = b**2*sin(c + d*x)/(a*d*(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*b*(2*a**2 - b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_458():
    f = sec(c + d*x)**2/(a + b*cos(c + d*x))**2
    F = b**2*tan(c + d*x)/(a*d*(a + b*cos(c + d*x))*(a**2 - b**2)) + (a**2 - 2*b**2)*tan(c + d*x)/(a**2*d*(a**2 - b**2)) + 2*b**2*(3*a**2 - 2*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - 2*b*atanh(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_459():
    f = sec(c + d*x)**3/(a + b*cos(c + d*x))**2
    F = b**2*tan(c + d*x)*sec(c + d*x)/(a*d*(a + b*cos(c + d*x))*(a**2 - b**2)) + (a**2 - 3*b**2)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d*(a**2 - b**2)) - b*(2*a**2 - 3*b**2)*tan(c + d*x)/(a**3*d*(a**2 - b**2)) - 2*b**3*(4*a**2 - 3*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + (a**2 + 6*b**2)*atanh(sin(c + d*x))/(2*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_460():
    f = sec(c + d*x)**4/(a + b*cos(c + d*x))**2
    F = b**2*tan(c + d*x)*sec(c + d*x)**2/(a*d*(a + b*cos(c + d*x))*(a**2 - b**2)) + (a**2 - 4*b**2)*tan(c + d*x)*sec(c + d*x)**2/(3*a**2*d*(a**2 - b**2)) - b*(a**2 - 2*b**2)*tan(c + d*x)*sec(c + d*x)/(a**3*d*(a**2 - b**2)) + (2*a**4 + 7*a**2*b**2 - 12*b**4)*tan(c + d*x)/(3*a**4*d*(a**2 - b**2)) + 2*b**4*(5*a**2 - 4*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - b*(a**2 + 4*b**2)*atanh(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_461():
    f = cos(c + d*x)**5/(a + b*cos(c + d*x))**3
    F = -a**3*(12*a**4 - 29*a**2*b**2 + 20*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - a**2*sin(c + d*x)*cos(c + d*x)**3/(2*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) - a**2*(4*a**2 - 7*b**2)*sin(c + d*x)*cos(c + d*x)**2/(2*b**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - 3*a*(4*a**4 - 7*a**2*b**2 + 2*b**4)*sin(c + d*x)/(2*b**4*d*(a**2 - b**2)**2) + (6*a**4 - 10*a**2*b**2 + b**4)*sin(c + d*x)*cos(c + d*x)/(2*b**3*d*(a**2 - b**2)**2) + x*(12*a**2 + b**2)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_462():
    f = cos(c + d*x)**4/(a + b*cos(c + d*x))**3
    F = 3*a**3*(a**2 - 2*b**2)*sin(c + d*x)/(2*b**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - a**2*sin(c + d*x)*cos(c + d*x)**2/(2*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + 3*a**2*(2*a**4 - 5*a**2*b**2 + 4*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - 3*a*x/b**4 + (3*a**2 - 2*b**2)*sin(c + d*x)/(2*b**3*d*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_463():
    f = cos(c + d*x)**3/(a + b*cos(c + d*x))**3
    F = -a**2*sin(c + d*x)*cos(c + d*x)/(2*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) - a**2*(2*a**2 - 5*b**2)*sin(c + d*x)/(2*b**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - a*(2*a**4 - 5*a**2*b**2 + 6*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + x/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_464():
    f = cos(c + d*x)**2/(a + b*cos(c + d*x))**3
    F = -a**2*sin(c + d*x)/(2*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + a*(a**2 - 4*b**2)*sin(c + d*x)/(2*b*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) + (a**2 + 2*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_465():
    f = cos(c + d*x)/(a + b*cos(c + d*x))**3
    F = -3*a*b*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a*sin(c + d*x)/(d*(a + b*cos(c + d*x))**2*(2*a**2 - 2*b**2)) + (a**2 + 2*b**2)*sin(c + d*x)/(2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_466():
    f = (a + b*cos(c + d*x))**(-3)
    F = -3*a*b*sin(c + d*x)/(2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - b*sin(c + d*x)/(d*(a + b*cos(c + d*x))**2*(2*a**2 - 2*b**2)) + (2*a**2 + b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_467():
    f = sec(c + d*x)/(a + b*cos(c + d*x))**3
    F = b**2*sin(c + d*x)/(2*a*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + b**2*(5*a**2 - 2*b**2)*sin(c + d*x)/(2*a**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - b*(6*a**4 - 5*a**2*b**2 + 2*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + atanh(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_468():
    f = sec(c + d*x)**2/(a + b*cos(c + d*x))**3
    F = b**2*tan(c + d*x)/(2*a*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + 3*b**2*(2*a**2 - b**2)*tan(c + d*x)/(2*a**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) + (2*a**4 - 11*a**2*b**2 + 6*b**4)*tan(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) + 3*b**2*(4*a**4 - 5*a**2*b**2 + 2*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - 3*b*atanh(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_469():
    f = sec(c + d*x)**3/(a + b*cos(c + d*x))**3
    F = b**2*tan(c + d*x)*sec(c + d*x)/(2*a*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + b**2*(7*a**2 - 4*b**2)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) + (a**4 - 10*a**2*b**2 + 6*b**4)*tan(c + d*x)*sec(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) - 3*b*(2*a**4 - 7*a**2*b**2 + 4*b**4)*tan(c + d*x)/(2*a**4*d*(a**2 - b**2)**2) - b**3*(20*a**4 - 29*a**2*b**2 + 12*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + (a**2 + 12*b**2)*atanh(sin(c + d*x))/(2*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_470():
    f = cos(c + d*x)**5/(a + b*cos(c + d*x))**4
    F = a**3*(4*a**4 - 11*a**2*b**2 + 12*b**4)*sin(c + d*x)/(2*b**4*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - a**2*sin(c + d*x)*cos(c + d*x)**3/(3*b*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) - a**2*(4*a**2 - 9*b**2)*sin(c + d*x)*cos(c + d*x)**2/(6*b**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + a**2*(8*a**6 - 28*a**4*b**2 + 35*a**2*b**4 - 20*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - 4*a*x/b**5 + (12*a**4 - 23*a**2*b**2 + 6*b**4)*sin(c + d*x)/(6*b**4*d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_471():
    f = cos(c + d*x)**4/(a + b*cos(c + d*x))**4
    F = a**3*(3*a**2 - 8*b**2)*sin(c + d*x)/(6*b**3*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) - a**2*sin(c + d*x)*cos(c + d*x)**2/(3*b*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) - a**2*(9*a**4 - 28*a**2*b**2 + 34*b**4)*sin(c + d*x)/(6*b**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - a*(2*a**6 - 7*a**4*b**2 + 8*a**2*b**4 - 8*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + x/b**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_472():
    f = cos(c + d*x)**3/(a + b*cos(c + d*x))**4
    F = -a**2*sin(c + d*x)*cos(c + d*x)/(3*b*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) - a**2*(2*a**2 - 7*b**2)*sin(c + d*x)/(6*b**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + a*(2*a**4 - 5*a**2*b**2 + 18*b**4)*sin(c + d*x)/(6*b**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - b*(3*a**2 + 2*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_473():
    f = cos(c + d*x)**2/(a + b*cos(c + d*x))**4
    F = -a**2*sin(c + d*x)/(3*b*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + a*(a**2 + 4*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + a*(a**2 - 6*b**2)*sin(c + d*x)/(6*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + (a**4 - 10*a**2*b**2 - 6*b**4)*sin(c + d*x)/(6*b*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_474():
    f = cos(c + d*x)/(a + b*cos(c + d*x))**4
    F = a*(2*a**2 + 13*b**2)*sin(c + d*x)/(6*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) + a*sin(c + d*x)/(d*(a + b*cos(c + d*x))**3*(3*a**2 - 3*b**2)) - b*(4*a**2 + b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (2*a**2 + 3*b**2)*sin(c + d*x)/(6*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_475():
    f = (a + b*cos(c + d*x))**(-4)
    F = -5*a*b*sin(c + d*x)/(6*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + a*(2*a**2 + 3*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - b*(11*a**2 + 4*b**2)*sin(c + d*x)/(6*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - b*sin(c + d*x)/(d*(a + b*cos(c + d*x))**3*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_476():
    f = sec(c + d*x)/(a + b*cos(c + d*x))**4
    F = b**2*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + b**2*(8*a**2 - 3*b**2)*sin(c + d*x)/(6*a**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + b**2*(26*a**4 - 17*a**2*b**2 + 6*b**4)*sin(c + d*x)/(6*a**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - b*(8*a**6 - 8*a**4*b**2 + 7*a**2*b**4 - 2*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + atanh(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_477():
    f = sec(c + d*x)**2/(a + b*cos(c + d*x))**4
    F = b**2*tan(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + b**2*(9*a**2 - 4*b**2)*tan(c + d*x)/(6*a**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + b**2*(12*a**4 - 11*a**2*b**2 + 4*b**4)*tan(c + d*x)/(2*a**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) + (6*a**6 - 65*a**4*b**2 + 68*a**2*b**4 - 24*b**6)*tan(c + d*x)/(6*a**4*d*(a**2 - b**2)**3) + b**2*(20*a**6 - 35*a**4*b**2 + 28*a**2*b**4 - 8*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - 4*b*atanh(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_478():
    f = sqrt(a + b*cos(c + d*x))*cos(c + d*x)**3
    F = -8*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(35*b**2*d) + 2*a*sqrt(a + b*cos(c + d*x))*(8*a**2 + 19*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(105*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)/(7*b*d) + sqrt(a + b*cos(c + d*x))*(16*a**2 + 50*b**2)*sin(c + d*x)/(105*b**2*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(16*a**4 + 34*a**2*b**2 - 50*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(105*b**3*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_479():
    f = sqrt(a + b*cos(c + d*x))*cos(c + d*x)**2
    F = -4*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(15*b*d) + 4*a*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b**2*d*sqrt(a + b*cos(c + d*x))) + 2*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d) - sqrt(a + b*cos(c + d*x))*(4*a**2 - 18*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_480():
    f = sqrt(a + b*cos(c + d*x))*cos(c + d*x)
    F = 2*a*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_481():
    f = sqrt(a + b*cos(c + d*x))
    F = 2*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_482():
    f = sqrt(a + b*cos(c + d*x))*sec(c + d*x)
    F = ((Integer(2) * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_483():
    f = sqrt(a + b*cos(c + d*x))*sec(c + d*x)**2
    F = (Integer(-1) * ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_484():
    f = sqrt(a + b*cos(c + d*x))*sec(c + d*x)**3
    F = ((Integer(-1) * (Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_485():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**3
    F = -8*a*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(63*b**2*d) + 2*a*sqrt(a + b*cos(c + d*x))*(8*a**2 + 39*b**2)*sin(c + d*x)/(315*b**2*d) - 2*a*sqrt((a + b*cos(c + d*x))/(a + b))*(8*a**4 + 31*a**2*b**2 - 39*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(315*b**3*d*sqrt(a + b*cos(c + d*x))) + 2*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)*cos(c + d*x)/(9*b*d) + (a + b*cos(c + d*x))**(sympy.S(3)/2)*(16*a**2 + 98*b**2)*sin(c + d*x)/(315*b**2*d) + sqrt(a + b*cos(c + d*x))*(16*a**4 + 66*a**2*b**2 + 294*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(315*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_486():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2
    F = -4*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(35*b*d) - 4*a*sqrt(a + b*cos(c + d*x))*(3*a**2 - 41*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(105*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d) - sqrt(a + b*cos(c + d*x))*(12*a**2 - 50*b**2)*sin(c + d*x)/(105*b*d) + sqrt((a + b*cos(c + d*x))/(a + b))*(12*a**4 - 62*a**2*b**2 + 50*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(105*b**2*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_487():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*d) - 2*a*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(5*b*d*sqrt(a + b*cos(c + d*x))) + 2*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + sqrt(a + b*cos(c + d*x))*(2*a**2 + 6*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(5*b*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_488():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 8*a*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_489():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)
    F = ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_490():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_491():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**3
    F = ((Integer(-5) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(7) * Symbol('a') * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(5) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_492():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**3
    F = -8*a*(a + b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(99*b**2*d) + 2*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(8*a**2 + 67*b**2)*sin(c + d*x)/(693*b**2*d) + 2*a*sqrt(a + b*cos(c + d*x))*(8*a**4 + 51*a**2*b**2 + 741*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(693*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*(a + b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)*cos(c + d*x)/(11*b*d) + (a + b*cos(c + d*x))**(sympy.S(5)/2)*(16*a**2 + 162*b**2)*sin(c + d*x)/(693*b**2*d) + sqrt(a + b*cos(c + d*x))*(16*a**4 + 114*a**2*b**2 + 270*b**4)*sin(c + d*x)/(693*b**2*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(16*a**6 + 98*a**4*b**2 + 156*a**2*b**4 - 270*b**6)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(693*b**3*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_493():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2
    F = -4*a*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(63*b*d) - 4*a*sqrt(a + b*cos(c + d*x))*(5*a**2 - 57*b**2)*sin(c + d*x)/(315*b*d) + 4*a*sqrt((a + b*cos(c + d*x))/(a + b))*(5*a**4 - 62*a**2*b**2 + 57*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(315*b**2*d*sqrt(a + b*cos(c + d*x))) + 2*(a + b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b*d) - (a + b*cos(c + d*x))**(sympy.S(3)/2)*(20*a**2 - 98*b**2)*sin(c + d*x)/(315*b*d) - sqrt(a + b*cos(c + d*x))*(20*a**4 - 558*a**2*b**2 - 294*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(315*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_494():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)
    F = 2*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(7*d) + 2*a*sqrt(a + b*cos(c + d*x))*(3*a**2 + 29*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(21*b*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) + sqrt(a + b*cos(c + d*x))*(6*a**2 + 10*b**2)*sin(c + d*x)/(21*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(6*a**4 + 4*a**2*b**2 - 10*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(21*b*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_495():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 16*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(15*d) - 16*a*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*d*sqrt(a + b*cos(c + d*x))) + 2*b*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + sqrt(a + b*cos(c + d*x))*(46*a**2 + 18*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_496():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)
    F = ((Integer(14) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_497():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2
    F = (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_498():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**3
    F = ((Integer(-9) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(8) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(15) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(9) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_499():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**4
    F = ((Integer(-1) * (((Integer(16) * (Symbol('a'))**(Integer(2))) + (Integer(33) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(16) * (Symbol('a'))**(Integer(2))) + (Integer(59) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(5) * Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(16) * (Symbol('a'))**(Integer(2))) + (Integer(33) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_500():
    f = (a + b*cos(c + d*x))**(sympy.S(7)/2)
    F = 24*a*b*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(35*d) + 32*a*sqrt(a + b*cos(c + d*x))*(11*a**2 + 13*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(105*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*b*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) + 2*b*sqrt(a + b*cos(c + d*x))*(71*a**2 + 25*b**2)*sin(c + d*x)/(105*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(142*a**4 - 92*a**2*b**2 - 50*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(105*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_501():
    f = sqrt(4*cos(c + d*x) + 3)*cos(c + d*x)**3
    F = (4*cos(c + d*x) + 3)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)/(14*d) - 3*(4*cos(c + d*x) + 3)**(sympy.S(3)/2)*sin(c + d*x)/(70*d) + 59*sqrt(4*cos(c + d*x) + 3)*sin(c + d*x)/(105*d) + 47*sqrt(7)*elliptic_e(c/2 + d*x/2, sympy.S(8)/7)/(140*d) + 59*sqrt(7)*elliptic_f(c/2 + d*x/2, sympy.S(8)/7)/(420*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_502():
    f = sqrt(4*cos(c + d*x) + 3)*cos(c + d*x)**2
    F = (4*cos(c + d*x) + 3)**(sympy.S(3)/2)*sin(c + d*x)/(10*d) - sqrt(4*cos(c + d*x) + 3)*sin(c + d*x)/(5*d) + 21*sqrt(7)*elliptic_e(c/2 + d*x/2, sympy.S(8)/7)/(20*d) - sqrt(7)*elliptic_f(c/2 + d*x/2, sympy.S(8)/7)/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_503():
    f = sqrt(4*cos(c + d*x) + 3)*cos(c + d*x)
    F = 2*sqrt(4*cos(c + d*x) + 3)*sin(c + d*x)/(3*d) + sqrt(7)*elliptic_e(c/2 + d*x/2, sympy.S(8)/7)/(2*d) + sqrt(7)*elliptic_f(c/2 + d*x/2, sympy.S(8)/7)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_504():
    f = sqrt(4*cos(c + d*x) + 3)
    F = 2*sqrt(7)*elliptic_e(c/2 + d*x/2, sympy.S(8)/7)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_505():
    f = sqrt(4*cos(c + d*x) + 3)*sec(c + d*x)
    F = ((Integer(8) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + ((Integer(6) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_506():
    f = sqrt(4*cos(c + d*x) + 3)*sec(c + d*x)**2
    F = (Integer(-1) * ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(3) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + ((Integer(4) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Integer(3) + (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_507():
    f = sqrt(4*cos(c + d*x) + 3)*sec(c + d*x)**3
    F = (Integer(-1) * ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + ((Integer(5) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Integer(3) + (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Integer(3) + (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_508():
    f = sqrt(3 - 4*cos(c + d*x))*cos(c + d*x)**3
    F = -(3 - 4*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)/(14*d) - 3*(3 - 4*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(70*d) + 59*sqrt(3 - 4*cos(c + d*x))*sin(c + d*x)/(105*d) - 47*sqrt(7)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(140*d) - 59*sqrt(7)*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(420*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_509():
    f = sqrt(3 - 4*cos(c + d*x))*cos(c + d*x)**2
    F = -(3 - 4*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(10*d) + sqrt(3 - 4*cos(c + d*x))*sin(c + d*x)/(5*d) + 21*sqrt(7)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(20*d) - sqrt(7)*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_510():
    f = sqrt(3 - 4*cos(c + d*x))*cos(c + d*x)
    F = 2*sqrt(3 - 4*cos(c + d*x))*sin(c + d*x)/(3*d) - sqrt(7)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(2*d) - sqrt(7)*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_511():
    f = sqrt(3 - 4*cos(c + d*x))
    F = 2*sqrt(7)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_512():
    f = sqrt(3 - 4*cos(c + d*x))*sec(c + d*x)
    F = (Integer(-1) * ((Integer(8) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(6) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_513():
    f = sqrt(3 - 4*cos(c + d*x))*sec(c + d*x)**2
    F = (Integer(-1) * ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(3) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + ((Integer(4) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_514():
    f = sqrt(3 - 4*cos(c + d*x))*sec(c + d*x)**3
    F = ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_515():
    f = cos(c + d*x)**3/sqrt(a + b*cos(c + d*x))
    F = -8*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(15*b**2*d) - 2*a*sqrt((a + b*cos(c + d*x))/(a + b))*(8*a**2 + 7*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b**3*d*sqrt(a + b*cos(c + d*x))) + 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)/(5*b*d) + sqrt(a + b*cos(c + d*x))*(16*a**2 + 18*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_516():
    f = cos(c + d*x)**2/sqrt(a + b*cos(c + d*x))
    F = -4*a*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b))) + 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*b*d) + sqrt((a + b*cos(c + d*x))/(a + b))*(4*a**2 + 2*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_517():
    f = cos(c + d*x)/sqrt(a + b*cos(c + d*x))
    F = -2*a*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt(a + b*cos(c + d*x))) + 2*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_518():
    f = 1/sqrt(a + b*cos(c + d*x))
    F = 2*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_519():
    f = sec(c + d*x)/sqrt(a + b*cos(c + d*x))
    F = (Integer(2) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_520():
    f = sec(c + d*x)**2/sqrt(a + b*cos(c + d*x))
    F = (Integer(-1) * ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_521():
    f = sec(c + d*x)**3/sqrt(a + b*cos(c + d*x))
    F = ((Integer(3) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_522():
    f = cos(c + d*x)**4/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*sin(c + d*x)*cos(c + d*x)**2/(b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*a*sqrt(a + b*cos(c + d*x))*(8*a**2 - 3*b**2)*sin(c + d*x)/(5*b**3*d*(a**2 - b**2)) - 8*a*sqrt((a + b*cos(c + d*x))/(a + b))*(4*a**2 + b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(5*b**4*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(12*a**2 - 2*b**2)*sin(c + d*x)*cos(c + d*x)/(5*b**2*d*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(32*a**4 - 16*a**2*b**2 - 6*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(5*b**4*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_523():
    f = cos(c + d*x)**3/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*sin(c + d*x)*cos(c + d*x)/(b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*a*sqrt(a + b*cos(c + d*x))*(8*a**2 - 5*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(8*a**2 - 2*b**2)*sin(c + d*x)/(3*b**2*d*(a**2 - b**2)) + sqrt((a + b*cos(c + d*x))/(a + b))*(16*a**2 + 2*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_524():
    f = cos(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*sin(c + d*x)/(b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - 4*a*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(b**2*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(4*a**2 - 2*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(b**2*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_525():
    f = cos(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*a*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)) + 2*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_526():
    f = (a + b*cos(c + d*x))**(sympy.S(-3)/2)
    F = -2*b*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + 2*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_527():
    f = sec(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(-2) * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_528():
    f = sec(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (sympy.tan((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_529():
    f = sec(c + d*x)**3/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(15) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_530():
    f = cos(c + d*x)**5/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*sin(c + d*x)*cos(c + d*x)**3/(3*b*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - 8*a**2*(2*a**2 - 3*b**2)*sin(c + d*x)*cos(c + d*x)**2/(3*b**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - 4*a*sqrt(a + b*cos(c + d*x))*(32*a**4 - 49*a**2*b**2 + 7*b**4)*sin(c + d*x)/(15*b**4*d*(a**2 - b**2)**2) - 2*a*sqrt((a + b*cos(c + d*x))/(a + b))*(128*a**4 - 116*a**2*b**2 - 17*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b**5*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(96*a**4 - 142*a**2*b**2 + 6*b**4)*sin(c + d*x)*cos(c + d*x)/(15*b**3*d*(a**2 - b**2)**2) + sqrt(a + b*cos(c + d*x))*(256*a**6 - 424*a**4*b**2 + 110*a**2*b**4 + 18*b**6)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b**5*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_531():
    f = cos(c + d*x)**4/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 4*a**3*(3*a**2 - 5*b**2)*sin(c + d*x)/(3*b**3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - 2*a**2*sin(c + d*x)*cos(c + d*x)**2/(3*b*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - 8*a*sqrt(a + b*cos(c + d*x))*(4*a**4 - 7*a**2*b**2 + 2*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**4*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2) + sqrt(a + b*cos(c + d*x))*(4*a**2 - 2*b**2)*sin(c + d*x)/(3*b**3*d*(a**2 - b**2)) + sqrt((a + b*cos(c + d*x))/(a + b))*(32*a**4 - 32*a**2*b**2 - 2*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**4*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_532():
    f = cos(c + d*x)**3/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*sin(c + d*x)*cos(c + d*x)/(3*b*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - 8*a**2*(a**2 - 2*b**2)*sin(c + d*x)/(3*b**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - 2*a*sqrt((a + b*cos(c + d*x))/(a + b))*(8*a**2 - 9*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(16*a**4 - 30*a**2*b**2 + 6*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_533():
    f = cos(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*sin(c + d*x)/(3*b*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 4*a*(a**2 - 3*b**2)*sin(c + d*x)/(3*b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - 4*a*sqrt(a + b*cos(c + d*x))*(a**2 - 3*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2) + sqrt((a + b*cos(c + d*x))/(a + b))*(4*a**2 - 6*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_534():
    f = cos(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + 2*a*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + (2*a**2 + 6*b**2)*sin(c + d*x)/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt(a + b*cos(c + d*x))*(-2*a**2 - 6*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_535():
    f = (a + b*cos(c + d*x))**(sympy.S(-5)/2)
    F = -8*a*b*sin(c + d*x)/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + 8*a*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2) - 2*b*sin(c + d*x)/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) - 2*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(d*sqrt(a + b*cos(c + d*x))*(3*a**2 - 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_536():
    f = sec(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(-2) * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_537():
    f = sec(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (sympy.tan((Symbol('c') + (Symbol('d') * x))) * ((Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_538():
    f = (a + b*cos(c + d*x))**(sympy.S(-7)/2)
    F = -16*a*b*sin(c + d*x)/(15*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)**2) - 16*a*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - 2*b*(23*a**2 + 9*b**2)*sin(c + d*x)/(15*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**3) - 2*b*sin(c + d*x)/(d*(a + b*cos(c + d*x))**(sympy.S(5)/2)*(5*a**2 - 5*b**2)) + sqrt(a + b*cos(c + d*x))*(46*a**2 + 18*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_539():
    f = cos(c + d*x)**3/sqrt(4*cos(c + d*x) + 3)
    F = sqrt(4*cos(c + d*x) + 3)*sin(c + d*x)*cos(c + d*x)/(10*d) - sqrt(4*cos(c + d*x) + 3)*sin(c + d*x)/(10*d) + 9*sqrt(7)*elliptic_e(c/2 + d*x/2, sympy.S(8)/7)/(20*d) - 23*sqrt(7)*elliptic_f(c/2 + d*x/2, sympy.S(8)/7)/(140*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_540():
    f = cos(c + d*x)**2/sqrt(4*cos(c + d*x) + 3)
    F = sqrt(4*cos(c + d*x) + 3)*sin(c + d*x)/(6*d) - sqrt(7)*elliptic_e(c/2 + d*x/2, sympy.S(8)/7)/(4*d) + 17*sqrt(7)*elliptic_f(c/2 + d*x/2, sympy.S(8)/7)/(84*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_541():
    f = cos(c + d*x)/sqrt(4*cos(c + d*x) + 3)
    F = sqrt(7)*elliptic_e(c/2 + d*x/2, sympy.S(8)/7)/(2*d) - 3*sqrt(7)*elliptic_f(c/2 + d*x/2, sympy.S(8)/7)/(14*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_542():
    f = 1/sqrt(4*cos(c + d*x) + 3)
    F = 2*sqrt(7)*elliptic_f(c/2 + d*x/2, sympy.S(8)/7)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_543():
    f = sec(c + d*x)/sqrt(4*cos(c + d*x) + 3)
    F = (Integer(2) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_544():
    f = sec(c + d*x)**2/sqrt(4*cos(c + d*x) + 3)
    F = (Integer(-1) * ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + (sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1)))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Integer(3) + (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_545():
    f = sec(c + d*x)**3/sqrt(4*cos(c + d*x) + 3)
    F = ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1)))) * ((Integer(3) * sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt(Integer(7)) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Integer(3) + (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Integer(3) + (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(6) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_546():
    f = cos(c + d*x)**3/sqrt(3 - 4*cos(c + d*x))
    F = -sqrt(3 - 4*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)/(10*d) - sqrt(3 - 4*cos(c + d*x))*sin(c + d*x)/(10*d) - 9*sqrt(7)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(20*d) + 23*sqrt(7)*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(140*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_547():
    f = cos(c + d*x)**2/sqrt(3 - 4*cos(c + d*x))
    F = -sqrt(3 - 4*cos(c + d*x))*sin(c + d*x)/(6*d) - sqrt(7)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(4*d) + 17*sqrt(7)*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(84*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_548():
    f = cos(c + d*x)/sqrt(3 - 4*cos(c + d*x))
    F = -sqrt(7)*elliptic_e(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(2*d) + 3*sqrt(7)*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(14*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_549():
    f = 1/sqrt(3 - 4*cos(c + d*x))
    F = 2*sqrt(7)*elliptic_f(c/2 + d*x/2 + pi/2, sympy.S(8)/7)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_550():
    f = sec(c + d*x)/sqrt(3 - 4*cos(c + d*x))
    F = Integer(-1) * ((Integer(2) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_551():
    f = sec(c + d*x)**2/sqrt(3 - 4*cos(c + d*x))
    F = (Integer(-1) * ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + (sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1)))) * ((sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(4) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_552():
    f = sec(c + d*x)**3/sqrt(3 - 4*cos(c + d*x))
    F = (Integer(-1) * ((sympy.sqrt(Integer(7)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + (sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1)))) * ((Integer(3) * sympy.sqrt(Integer(7)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt(Integer(7)) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + sympy.pi + (Symbol('d') * x))), (Integer(8) * (Integer(7))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(4) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(6) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_553():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*A*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*B*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 10*B*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 10*B*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_554():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*A*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*B*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_555():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))
    F = 2*A*elliptic_e(c/2 + d*x/2, 2)/d + 2*B*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*B*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_556():
    f = (A + B*cos(c + d*x))/sqrt(cos(c + d*x))
    F = 2*A*elliptic_f(c/2 + d*x/2, 2)/d + 2*B*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_557():
    f = (A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*A*elliptic_e(c/2 + d*x/2, 2)/d + 2*B*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_558():
    f = (A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*A*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*B*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*B*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_559():
    f = (A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = 6*A*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) + 2*A*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) - 6*A*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*B*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*B*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_560():
    f = (a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2)
    F = 4*a*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 20*a*b*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 20*a*b*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*b**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + (18*a**2 + 14*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + (18*a**2 + 14*b**2)*elliptic_e(c/2 + d*x/2, 2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_561():
    f = (a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 4*a*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 12*a*b*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*b**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + (14*a**2 + 10*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (14*a**2 + 10*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_562():
    f = (a + b*cos(c + d*x))**2*sqrt(cos(c + d*x))
    F = 4*a*b*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 4*a*b*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*b**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + (10*a**2 + 6*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_563():
    f = (a + b*cos(c + d*x))**2/sqrt(cos(c + d*x))
    F = 4*a*b*elliptic_e(c/2 + d*x/2, 2)/d + 2*b**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + (6*a**2 + 2*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_564():
    f = (a + b*cos(c + d*x))**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 4*a*b*elliptic_f(c/2 + d*x/2, 2)/d - (2*a**2 - 2*b**2)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_565():
    f = (a + b*cos(c + d*x))**2/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a*b*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 4*a*b*elliptic_e(c/2 + d*x/2, 2)/d + (2*a**2 + 6*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_566():
    f = (a + b*cos(c + d*x))**2/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 4*a*b*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a*b*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (6*a**2 + 10*b**2)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - (6*a**2 + 10*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_567():
    f = (a + b*cos(c + d*x))**3*cos(c + d*x)**(sympy.S(3)/2)
    F = 40*a*b**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d) + 2*a*(7*a**2 + 15*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*a*(7*a**2 + 15*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*b**2*(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(9*d) + 2*b*(27*a**2 + 7*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 2*b*(27*a**2 + 7*b**2)*elliptic_e(c/2 + d*x/2, 2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_568():
    f = (a + b*cos(c + d*x))**3*sqrt(cos(c + d*x))
    F = 32*a*b**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + 2*a*(5*a**2 + 9*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*b**2*(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*b*(21*a**2 + 5*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*b*(21*a**2 + 5*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_569():
    f = (a + b*cos(c + d*x))**3/sqrt(cos(c + d*x))
    F = 8*a*b**2*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 2*a*(a**2 + b**2)*elliptic_f(c/2 + d*x/2, 2)/d + 2*b**2*(a + b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 6*b*(5*a**2 + b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_570():
    f = (a + b*cos(c + d*x))**3/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*(a + b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*a*(a**2 - 3*b**2)*elliptic_e(c/2 + d*x/2, 2)/d - 2*b*(3*a**2 - b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*b*(9*a**2 + b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_571():
    f = (a + b*cos(c + d*x))**3/cos(c + d*x)**(sympy.S(5)/2)
    F = 16*a**2*b*sin(c + d*x)/(3*d*sqrt(cos(c + d*x))) + 2*a**2*(a + b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(a**2 + 9*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d) - 2*b*(3*a**2 - b**2)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_572():
    f = (a + b*cos(c + d*x))**3/cos(c + d*x)**(sympy.S(7)/2)
    F = 8*a**2*b*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(a + b*cos(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 6*a*(a**2 + 5*b**2)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 6*a*(a**2 + 5*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*b*(a**2 + b**2)*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_573():
    f = (a + b*cos(c + d*x))**3/cos(c + d*x)**(sympy.S(9)/2)
    F = 32*a**2*b*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(a + b*cos(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 2*a*(5*a**2 + 21*b**2)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(5*a**2 + 21*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*b*(9*a**2 + 5*b**2)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 2*b*(9*a**2 + 5*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_574():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_575():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))
    F = ((Integer(2) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_576():
    f = sqrt(cos(c + d*x))/(a + b*cos(c + d*x))
    F = ((Integer(2) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_577():
    f = 1/((a + b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = (Integer(2) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_578():
    f = 1/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(2) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_579():
    f = 1/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(2) * Symbol('b') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_580():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(4))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_581():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**2
    F = ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_582():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((Symbol('a') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_583():
    f = sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**2
    F = (sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_584():
    f = 1/((a + b*cos(c + d*x))**2*sqrt(cos(c + d*x)))
    F = (Integer(-1) * ((Symbol('b') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_585():
    f = 1/((a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_586():
    f = 1/((a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_587():
    f = cos(c + d*x)**(sympy.S(9)/2)/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((Symbol('a') * ((Integer(35) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(24) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(105) * (Symbol('a'))**(Integer(6))) + (Integer(-1) * (Integer(223) * (Symbol('a'))**(Integer(4)) * (Symbol('b'))**(Integer(2)))) + (Integer(128) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(4))) + (Integer(8) * (Symbol('b'))**(Integer(6)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(12) * (Symbol('b'))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(3)) * ((Integer(35) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(63) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(5)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(35) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(13) * (Symbol('b'))**(Integer(2))))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_588():
    f = cos(c + d*x)**(sympy.S(7)/2)/(a + b*cos(c + d*x))**3
    F = ((((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(11) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Integer(15) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(11) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_589():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_590():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4))))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_591():
    f = sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**3
    F = ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_592():
    f = 1/((a + b*cos(c + d*x))**3*sqrt(cos(c + d*x)))
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * ((Integer(5) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Symbol('b'))**(Integer(4))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_593():
    f = 1/((a + b*cos(c + d*x))**3*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(35) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_594():
    f = 1/((a + b*cos(c + d*x))**3*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Symbol('b') * ((Integer(24) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(63) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(24) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(13) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_595():
    f = sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_596():
    f = sqrt(a + b*cos(c + d*x))*sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_597():
    f = sqrt(a + b*cos(c + d*x))/sqrt(cos(c + d*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_598():
    f = sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_599():
    f = sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_600():
    f = sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(15*a*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a + 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a**2 - 2*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_601():
    f = sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 2*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(35*a*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(50*a**2 - 8*b**2)*sin(c + d*x)/(105*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(25*a**2 + 6*a*b + 8*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**3*d) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(19*a**2 + 8*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_602():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * ((Integer(3) * Symbol('a')) + (Integer(8) * Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(12) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_603():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = (Integer(-1) * ((Integer(5) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(5) * Symbol('a')) + (Integer(2) * Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_604():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a')) + Symbol('b')) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_605():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * (Integer(2) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_606():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(8*a - 8*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*(2*a - 6*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_607():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 4*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3*a - b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(5*a*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3*a**2 + b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_608():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 16*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(50*a**2 + 6*b**2)*sin(c + d*x)/(105*a*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(25*a**2 - 57*a*b - 6*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**2*d) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(4*a - 4*b)*(41*a**2 - 3*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_609():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2)) + 20*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(63*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*cos(c + d*x))*(98*a**2 + 6*b**2)*sin(c + d*x)/(315*a*d*cos(c + d*x)**(sympy.S(5)/2)) + 8*b*sqrt(a + b*cos(c + d*x))*(22*a**2 - b**2)*sin(c + d*x)/(315*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*a**3 - 39*a**2*b - 6*a*b**2 - 8*b**3)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**3*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*a**4 + 33*a**2*b**2 + 8*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_610():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(26) * Symbol('a') * Symbol('b')) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_611():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = (Integer(-1) * ((Integer(9) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(9) * Symbol('a') * Symbol('b')) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + ((Integer(9) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_612():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('a') * Symbol('b'))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_613():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(14) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('b'))) + (Integer(9) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_614():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 22*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(15*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a**2 + 23*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a**2 - 8*a*b + 15*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_615():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 6*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(10*a**2 + 18*b**2)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(5*a**2 - 24*a*b + 3*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(21*a*d) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(29*a**2 + 3*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(21*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_616():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = 2*a**2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2)) + 38*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(63*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*cos(c + d*x))*(98*a**2 + 150*b**2)*sin(c + d*x)/(315*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*b*sqrt(a + b*cos(c + d*x))*(163*a**2 + 5*b**2)*sin(c + d*x)/(315*a*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*a**3 - 114*a**2*b + 165*a*b**2 + 10*b**3)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*a**4 + 279*a**2*b**2 - 10*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_617():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(13)/2)
    F = 2*a**2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(11*d*cos(c + d*x)**(sympy.S(11)/2)) + 46*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(99*d*cos(c + d*x)**(sympy.S(9)/2)) + sqrt(a + b*cos(c + d*x))*(162*a**2 + 226*b**2)*sin(c + d*x)/(693*d*cos(c + d*x)**(sympy.S(7)/2)) + 2*b*sqrt(a + b*cos(c + d*x))*(229*a**2 + 3*b**2)*sin(c + d*x)/(693*a*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(270*a**4 + 410*a**2*b**2 - 8*b**4)*sin(c + d*x)/(693*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(135*a**4 - 606*a**3*b + 57*a**2*b**2 + 6*a*b**3 + 8*b**4)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(693*a**3*d) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(741*a**4 + 51*a**2*b**2 + 8*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(693*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_618():
    f = sqrt(cos(c + d*x))/sqrt(a + b*cos(c + d*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_619():
    f = 1/(sqrt(a + b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_620():
    f = 1/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = -2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_621():
    f = 1/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*a*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(a + 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d) - b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(4*a - 4*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_622():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * Symbol('a')) + Symbol('b')) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_623():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_624():
    f = sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b)) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_625():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = -2*b*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b)) + 2*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_626():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*b**2*sin(c + d*x)/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a + 4*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a**2 - 4*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**3*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_627():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*b**2*sin(c + d*x)/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*cos(c + d*x))*(2*a**2 - 8*b**2)*sin(c + d*x)/(3*a**2*d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(a + 4*b)*(2*a + 4*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(a + b)) - 2*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(5*a**2 - 8*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_628():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(7)/2))
    F = 2*b**2*sin(c + d*x)/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(2*a**2 - 12*b**2)*sin(c + d*x)/(5*a**2*d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(5)/2)) - 2*b*sqrt(a + b*cos(c + d*x))*(3*a**2 - 8*b**2)*sin(c + d*x)/(5*a**3*d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a + 8*b)*(a**2 + 4*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(5*a**4*d*sqrt(a + b)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**4 + 16*a**2*b**2 - 32*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(5*a**5*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_629():
    f = cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Symbol('a') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_630():
    f = cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -8*a*b*sin(c + d*x)/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + 2*a*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + 8*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a - 6*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_631():
    f = sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*b*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + (6*a**2 + 2*b**2)*sin(c + d*x)/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**2 + 2*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_632():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = 2*b**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - 4*b*(3*a**2 - b**2)*sin(c + d*x)/(3*a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**2 - 6*a*b - 4*b**2)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*(a - b)*(a + b)**(sympy.S(3)/2)) + 4*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(3*a**2 - b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_633():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*b**2*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(cos(c + d*x))) + 8*b**2*(2*a**2 - b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**3 + 18*a**2*b - 12*a*b**2 - 16*b**3)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*(a - b)*(a + b)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**4 - 30*a**2*b**2 + 16*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_634():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*b**2*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*cos(c + d*x)**(sympy.S(3)/2)) + 4*b**2*(5*a**2 - 3*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*cos(c + d*x))*(2*a**4 - 26*a**2*b**2 + 16*b**4)*sin(c + d*x)/(3*a**3*d*(a**2 - b**2)**2*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a**4 + 18*a**3*b + 32*a**2*b**2 - 24*a*b**3 - 32*b**4)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*(a - b)*(a + b)**(sympy.S(3)/2)) - 8*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a**4 - 7*a**2*b**2 + 4*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**5*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_635():
    f = 1/(sqrt(3*cos(c + d*x) + 2)*sqrt(cos(c + d*x)))
    F = 2*sqrt(5)*elliptic_f(asin(sin(c + d*x)/(cos(c + d*x) + 1)), sympy.S(1)/5)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_636():
    f = 1/(sqrt(3*cos(c + d*x) - 2)*sqrt(cos(c + d*x)))
    F = 2*elliptic_f(asin(sin(c + d*x)/(cos(c + d*x) + 1)), 5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_637():
    f = 1/(sqrt(2 - 3*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = -2*sqrt(5)*sqrt(-cos(c + d*x))*elliptic_f(asin(sin(c + d*x)/(1 - cos(c + d*x))), sympy.S(1)/5)/(5*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_638():
    f = 1/(sqrt(-3*cos(c + d*x) - 2)*sqrt(cos(c + d*x)))
    F = -2*sqrt(-cos(c + d*x))*elliptic_f(asin(sin(c + d*x)/(1 - cos(c + d*x))), 5)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_639():
    f = 1/(sqrt(2*cos(c + d*x) + 3)*sqrt(cos(c + d*x)))
    F = 2*sqrt(-tan(c + d*x)**2)*cot(c + d*x)*elliptic_f(asin(sqrt(5)*sqrt(2*cos(c + d*x) + 3)/(5*sqrt(cos(c + d*x)))), -5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_640():
    f = 1/(sqrt(3 - 2*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = 2*sqrt(5)*sqrt(-tan(c + d*x)**2)*cot(c + d*x)*elliptic_f(asin(sqrt(3 - 2*cos(c + d*x))/sqrt(cos(c + d*x))), sympy.S(-1)/5)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_641():
    f = 1/(sqrt(2*cos(c + d*x) - 3)*sqrt(cos(c + d*x)))
    F = -2*sqrt(5)*sqrt(-cos(c + d*x))*sqrt(-tan(c + d*x)**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(2*cos(c + d*x) - 3)/sqrt(-cos(c + d*x))), sympy.S(-1)/5)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_642():
    f = 1/(sqrt(-2*cos(c + d*x) - 3)*sqrt(cos(c + d*x)))
    F = -2*sqrt(-cos(c + d*x))*sqrt(-tan(c + d*x)**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(5)*sqrt(-2*cos(c + d*x) - 3)/(5*sqrt(-cos(c + d*x)))), -5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_643():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(3*cos(c + d*x) + 2))
    F = 2*sqrt(5)*sqrt(cos(c + d*x))*elliptic_f(asin(sin(c + d*x)/(cos(c + d*x) + 1)), sympy.S(1)/5)/(5*d*sqrt(-cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_644():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(3*cos(c + d*x) - 2))
    F = 2*sqrt(cos(c + d*x))*elliptic_f(asin(sin(c + d*x)/(cos(c + d*x) + 1)), 5)/(d*sqrt(-cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_645():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(2 - 3*cos(c + d*x)))
    F = -2*sqrt(5)*elliptic_f(asin(sin(c + d*x)/(1 - cos(c + d*x))), sympy.S(1)/5)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_646():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(-3*cos(c + d*x) - 2))
    F = -2*elliptic_f(asin(sin(c + d*x)/(1 - cos(c + d*x))), 5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_647():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(2*cos(c + d*x) + 3))
    F = 2*sqrt(-tan(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)*csc(c + d*x)*elliptic_f(asin(sqrt(5)*sqrt(2*cos(c + d*x) + 3)/(5*sqrt(cos(c + d*x)))), -5)/(d*sqrt(-cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_648():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(3 - 2*cos(c + d*x)))
    F = 2*sqrt(5)*sqrt(-tan(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)*csc(c + d*x)*elliptic_f(asin(sqrt(3 - 2*cos(c + d*x))/sqrt(cos(c + d*x))), sympy.S(-1)/5)/(5*d*sqrt(-cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_649():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(2*cos(c + d*x) - 3))
    F = -2*sqrt(5)*sqrt(-tan(c + d*x)**2)*cot(c + d*x)*elliptic_f(asin(sqrt(2*cos(c + d*x) - 3)/sqrt(-cos(c + d*x))), sympy.S(-1)/5)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_650():
    f = 1/(sqrt(-cos(c + d*x))*sqrt(-2*cos(c + d*x) - 3))
    F = -2*sqrt(-tan(c + d*x)**2)*cot(c + d*x)*elliptic_f(asin(sqrt(5)*sqrt(-2*cos(c + d*x) - 3)/(5*sqrt(-cos(c + d*x)))), -5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_651():
    f = sqrt(cos(c + d*x))/sqrt(3*cos(c + d*x) + 2)
    F = Integer(-1) * ((Integer(4) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(3))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(2) + (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(5)) * sympy.sqrt((Integer(-1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_652():
    f = sqrt(cos(c + d*x))/sqrt(3*cos(c + d*x) - 2)
    F = Integer(-1) * ((Integer(4) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(3))**(Integer(-1)), sympy.asin((sympy.sqrt((Integer(-2) + (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))), (Integer(5))**(Integer(-1))) * sympy.sqrt((Integer(-1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * sympy.sqrt(Integer(5)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_653():
    f = sqrt(cos(c + d*x))/sqrt(2 - 3*cos(c + d*x))
    F = Integer(-1) * ((Integer(4) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(3))**(Integer(-1)), sympy.asin((sympy.sqrt((Integer(2) + (Integer(-1) * (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * (sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(5))**(Integer(-1))) * sympy.sqrt((Integer(-1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * sympy.sqrt(Integer(5)) * Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_654():
    f = sqrt(cos(c + d*x))/sqrt(-3*cos(c + d*x) - 2)
    F = Integer(-1) * ((Integer(4) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(3))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-2) + (Integer(-1) * (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), Integer(5)) * sympy.sqrt((Integer(-1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_655():
    f = sqrt(cos(c + d*x))/sqrt(2*cos(c + d*x) + 3)
    F = Integer(-1) * ((Integer(3) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(3) + (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-5)) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_656():
    f = sqrt(cos(c + d*x))/sqrt(3 - 2*cos(c + d*x))
    F = (Integer(3) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * (sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))), (Integer(-1) * (Integer(5))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_657():
    f = sqrt(cos(c + d*x))/sqrt(2*cos(c + d*x) - 3)
    F = (Integer(3) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-3) + (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * (Integer(5))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_658():
    f = sqrt(cos(c + d*x))/sqrt(-2*cos(c + d*x) - 3)
    F = Integer(-1) * ((Integer(3) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-3) + (Integer(-1) * (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), Integer(-5)) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('d') * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_659():
    f = sqrt(-cos(c + d*x))/sqrt(3*cos(c + d*x) + 2)
    F = Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(3))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(2) + (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(5)) * sympy.sqrt((Integer(-1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_660():
    f = sqrt(-cos(c + d*x))/sqrt(3*cos(c + d*x) - 2)
    F = Integer(-1) * ((Integer(4) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(3))**(Integer(-1)), sympy.asin((sympy.sqrt((Integer(-2) + (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))), (Integer(5))**(Integer(-1))) * sympy.sqrt((Integer(-1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * sympy.sqrt(Integer(5)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_661():
    f = sqrt(-cos(c + d*x))/sqrt(2 - 3*cos(c + d*x))
    F = Integer(-1) * ((Integer(4) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(3))**(Integer(-1)), sympy.asin((sympy.sqrt((Integer(2) + (Integer(-1) * (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * (sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(5))**(Integer(-1))) * sympy.sqrt((Integer(-1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Integer(3) * sympy.sqrt(Integer(5)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_662():
    f = sqrt(-cos(c + d*x))/sqrt(-3*cos(c + d*x) - 2)
    F = Integer(-1) * ((Integer(4) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(3))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-2) + (Integer(-1) * (Integer(3) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), Integer(5)) * sympy.sqrt((Integer(-1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x))))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_663():
    f = sqrt(-cos(c + d*x))/sqrt(2*cos(c + d*x) + 3)
    F = Integer(-1) * ((Integer(3) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(3) + (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), Integer(-5)) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_664():
    f = sqrt(-cos(c + d*x))/sqrt(3 - 2*cos(c + d*x))
    F = (Integer(3) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(3) + (Integer(-1) * (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * (sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))))**(Integer(-1)))), (Integer(-1) * (Integer(5))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_665():
    f = sqrt(-cos(c + d*x))/sqrt(2*cos(c + d*x) - 3)
    F = (Integer(3) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-3) + (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * (Integer(5))**(Integer(-1)))) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt(Integer(5)) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_666():
    f = sqrt(-cos(c + d*x))/sqrt(-2*cos(c + d*x) - 3)
    F = Integer(-1) * ((Integer(3) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Integer(5) * (Integer(2))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(-3) + (Integer(-1) * (Integer(2) * sympy.cos((Symbol('c') + (Symbol('d') * x))))))) * ((sympy.sqrt(Integer(5)) * sympy.sqrt((Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), Integer(-5)) * sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt((Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_667():
    f = cos(c + d*x)**(sympy.S(2)/3)/(a + b*cos(c + d*x))
    F = a*(cos(c + d*x)**2)**(sympy.S(1)/6)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(1)/3)) - b*sin(c + d*x)*cos(c + d*x)**(sympy.S(2)/3)*appellf1(sympy.S.Half, sympy.S(-1)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*(cos(c + d*x)**2)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_668():
    f = cos(c + d*x)**(sympy.S(1)/3)/(a + b*cos(c + d*x))
    F = a*(cos(c + d*x)**2)**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(2)/3)) - b*sin(c + d*x)*cos(c + d*x)**(sympy.S(1)/3)*appellf1(sympy.S.Half, sympy.S(-1)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*(cos(c + d*x)**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_669():
    f = 1/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(1)/3))
    F = a*(cos(c + d*x)**2)**(sympy.S(2)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(2)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(4)/3)) - b*(cos(c + d*x)**2)**(sympy.S(1)/6)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_670():
    f = 1/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(2)/3))
    F = a*(cos(c + d*x)**2)**(sympy.S(5)/6)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(5)/6, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(5)/3)) - b*(cos(c + d*x)**2)**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, 1, sympy.S(3)/2, sin(c + d*x)**2, -b**2*sin(c + d*x)**2/(a**2 - b**2))/(d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(2)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_671():
    f = cos(c + d*x)**(sympy.S(7)/3)/sqrt(a + b*cos(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_672():
    f = cos(c + d*x)**(sympy.S(5)/3)/sqrt(a + b*cos(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_673():
    f = cos(c + d*x)**(sympy.S(4)/3)/sqrt(a + b*cos(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_674():
    f = cos(c + d*x)**(sympy.S(2)/3)/sqrt(a + b*cos(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_675():
    f = cos(c + d*x)**(sympy.S(1)/3)/sqrt(a + b*cos(c + d*x))
    F = sympy.Function('Unintegrable')(((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_676():
    f = 1/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(1)/3))
    F = sympy.Function('Unintegrable')((((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3))**(Integer(-1))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_677():
    f = 1/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(2)/3))
    F = sympy.Function('Unintegrable')((((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_678():
    f = 1/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(4)/3))
    F = sympy.Function('Unintegrable')((((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_679():
    f = 1/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/3))
    F = sympy.Function('Unintegrable')((((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_680():
    f = 1/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/3))
    F = sympy.Function('Unintegrable')((((sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(3))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_681():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 6*A*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 6*A*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*B*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_682():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*B*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_683():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*A*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_684():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))
    F = 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_685():
    f = (A + B*cos(c + d*x))/sqrt(sec(c + d*x))
    F = 2*A*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*B*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_686():
    f = (A + B*cos(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*B*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_687():
    f = (A + B*cos(c + d*x))/sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*A*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 10*B*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*B*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 10*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_688():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 4*a*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 12*a*b*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 12*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (10*a**2 + 14*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + (10*a**2 + 14*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_689():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 4*a*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 4*a*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (6*a**2 + 10*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (6*a**2 + 10*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_690():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 4*a*b*sin(c + d*x)*sqrt(sec(c + d*x))/d - 4*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (2*a**2 + 6*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_691():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/d + 4*a*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d - (2*a**2 - 2*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_692():
    f = (a + b*cos(c + d*x))**2*sqrt(sec(c + d*x))
    F = 4*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + (6*a**2 + 2*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_693():
    f = (a + b*cos(c + d*x))**2/sqrt(sec(c + d*x))
    F = 4*a*b*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 4*a*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*b**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + (10*a**2 + 6*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_694():
    f = (a + b*cos(c + d*x))**2/sec(c + d*x)**(sympy.S(3)/2)
    F = 4*a*b*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 12*a*b*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + (14*a**2 + 10*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (14*a**2 + 10*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_695():
    f = (a + b*cos(c + d*x))**2/sec(c + d*x)**(sympy.S(5)/2)
    F = 20*a*b*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 4*a*b*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 20*a*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*b**2*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + (18*a**2 + 14*b**2)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + (18*a**2 + 14*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_696():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(9)/2)
    F = 32*a**2*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + 2*a**2*(a*sec(c + d*x) + b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*(5*a**2 + 21*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + 2*a*(5*a**2 + 21*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*b*(9*a**2 + 5*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 2*b*(9*a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_697():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(7)/2)
    F = 8*a**2*b*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a**2*(a*sec(c + d*x) + b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*a*(a**2 + 5*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 6*a*(a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b*(a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_698():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(5)/2)
    F = 16*a**2*b*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*a**2*(a*sec(c + d*x) + b)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*a*(a**2 + 9*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - 2*b*(3*a**2 - b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_699():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2)
    F = -2*a*(a**2 - 3*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*(3*a**2 - b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*b**2*(a*sec(c + d*x) + b)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*b*(9*a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_700():
    f = (a + b*cos(c + d*x))**3*sqrt(sec(c + d*x))
    F = 8*a*b**2*sin(c + d*x)/(5*d*sqrt(sec(c + d*x))) + 2*a*(a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*b**2*(a*sec(c + d*x) + b)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*b*(5*a**2 + b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_701():
    f = (a + b*cos(c + d*x))**3/sqrt(sec(c + d*x))
    F = 32*a*b**2*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*(5*a**2 + 9*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*(a*sec(c + d*x) + b)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*b*(21*a**2 + 5*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*b*(21*a**2 + 5*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_702():
    f = (a + b*cos(c + d*x))**3/sec(c + d*x)**(sympy.S(3)/2)
    F = 40*a*b**2*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a*(7*a**2 + 15*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*(7*a**2 + 15*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*b**2*(a*sec(c + d*x) + b)*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*b*(27*a**2 + 7*b**2)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*b*(27*a**2 + 7*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_703():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))
    F = ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_704():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_705():
    f = sqrt(sec(c + d*x))/(a + b*cos(c + d*x))
    F = (Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_706():
    f = 1/((a + b*cos(c + d*x))*sqrt(sec(c + d*x)))
    F = ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_707():
    f = 1/((a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_708():
    f = 1/((a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_709():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**2
    F = ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_710():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_711():
    f = sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_712():
    f = 1/((a + b*cos(c + d*x))**2*sqrt(sec(c + d*x)))
    F = ((sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_713():
    f = 1/((a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Symbol('a') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_714():
    f = 1/((a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(5)/2))
    F = ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_715():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**3
    F = ((Symbol('b') * ((Integer(24) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(63) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(24) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(35) * (Symbol('b'))**(Integer(4)))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(13) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_716():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(35) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(15) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_717():
    f = sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((Integer(3) * Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(3) * ((Integer(5) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Symbol('b'))**(Integer(4))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_718():
    f = 1/((a + b*cos(c + d*x))**3*sqrt(sec(c + d*x)))
    F = ((((Integer(5) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_719():
    f = 1/((a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Symbol('a'))**(Integer(2)) + (Integer(5) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_720():
    f = 1/((a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(5)/2))
    F = (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(8) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(4)) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(5) * (Symbol('b'))**(Integer(4)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_721():
    f = sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*a*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a + 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**2*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a**2 - 2*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_722():
    f = sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*sqrt(sec(c + d*x))) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_723():
    f = sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(sec(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_724():
    f = sqrt(a + b*cos(c + d*x))*sqrt(sec(c + d*x))
    F = Integer(-1) * (((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)) * (Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')((Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_725():
    f = sqrt(a + b*cos(c + d*x))/sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_726():
    f = sqrt(a + b*cos(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_727():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 16*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + sqrt(a + b*cos(c + d*x))*(50*a**2 + 6*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*a*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(25*a**2 - 57*a*b - 6*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**2*d*sqrt(sec(c + d*x))) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(4*a - 4*b)*(41*a**2 - 3*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_728():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 4*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3*a - b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(5*a*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3*a**2 + b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(5*a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_729():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(8*a - 8*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*(2*a - 6*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_730():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * (Integer(2) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_731():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a')) + Symbol('b')) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_732():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = ((Integer(-5) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(5) * Symbol('a')) + (Integer(2) * Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_733():
    f = (a + b*cos(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * ((Integer(3) * Symbol('a')) + (Integer(8) * Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(12) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_734():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(11)/2)
    F = 2*a**2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d) + 38*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d) + sqrt(a + b*cos(c + d*x))*(98*a**2 + 150*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(315*d) + 2*b*sqrt(a + b*cos(c + d*x))*(163*a**2 + 5*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*a*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*a**3 - 114*a**2*b + 165*a*b**2 + 10*b**3)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**2*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*a**4 + 279*a**2*b**2 - 10*b**4)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_735():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*a**2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 6*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(7*d) + sqrt(a + b*cos(c + d*x))*(10*a**2 + 18*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(5*a**2 - 24*a*b + 3*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(21*a*d*sqrt(sec(c + d*x))) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(29*a**2 + 3*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(21*a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_736():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*a**2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 22*a*b*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a**2 + 23*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a*d*sqrt(sec(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*a**2 - 8*a*b + 15*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_737():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(14) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('b'))) + (Integer(9) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_738():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(6) * Symbol('a') * Symbol('b'))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_739():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = ((Integer(-9) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(9) * Symbol('a') * Symbol('b')) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(9) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_740():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = ((Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(24) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(26) * Symbol('a') * Symbol('b')) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(13) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(33) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_741():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(-1) * ((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(284) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(192) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * (Symbol('a'))**(Integer(3))) + (Integer(118) * (Symbol('a'))**(Integer(2)) * Symbol('b')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2))) + (Integer(72) * (Symbol('b'))**(Integer(3)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(5) * (Symbol('a'))**(Integer(4))) + (Integer(-1) * (Integer(120) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(48) * (Symbol('b'))**(Integer(4))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(17) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(59) * (Symbol('a'))**(Integer(2))) + (Integer(36) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(96) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(15) * (Symbol('a'))**(Integer(2))) + (Integer(284) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_742():
    f = sec(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*cos(c + d*x))
    F = 2*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(a + 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*sqrt(sec(c + d*x))) + b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(-4*a + 4*b)*sqrt(a + b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_743():
    f = sec(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*cos(c + d*x))
    F = -2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_744():
    f = sqrt(sec(c + d*x))/sqrt(a + b*cos(c + d*x))
    F = 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_745():
    f = 1/(sqrt(a + b*cos(c + d*x))*sqrt(sec(c + d*x)))
    F = (Integer(-2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_746():
    f = 1/(sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (sympy.sin((Symbol('c') + (Symbol('d') * x))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_747():
    f = 1/(sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * Symbol('a')) + (Integer(-1) * (Integer(2) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_748():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(2*a**2 - 8*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(a**2 - b**2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(a + 4*b)*(2*a + 4*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(a + b)*sqrt(sec(c + d*x))) - 2*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(5*a**2 - 8*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*sqrt(a + b)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_749():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b**2*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a + 4*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a**2 - 4*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**3*d*sqrt(a + b)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_750():
    f = sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = -2*b*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b)*sqrt(sec(c + d*x))) + 2*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_751():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = 2*a*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b)*sqrt(sec(c + d*x))) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_752():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_753():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(3) * Symbol('a')) + Symbol('b')) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_754():
    f = sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 4*b**2*(5*a**2 - 3*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt(a + b*cos(c + d*x))*(2*a**4 - 26*a**2*b**2 + 16*b**4)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**3*d*(a**2 - b**2)**2) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a**4 + 18*a**3*b + 32*a**2*b**2 - 24*a*b**3 - 32*b**4)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) - 8*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a**4 - 7*a**2*b**2 + 4*b**4)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**5*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_755():
    f = sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b**2*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 8*b**2*(2*a**2 - b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**3 + 18*a**2*b - 12*a*b**2 - 16*b**3)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**4 - 30*a**2*b**2 + 16*b**4)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_756():
    f = sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b**2*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(sec(c + d*x))) - 4*b*(3*a**2 - b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**2 - 6*a*b - 4*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + 4*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(3*a**2 - b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_757():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = -2*b*sin(c + d*x)/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*sqrt(sec(c + d*x))) + (6*a**2 + 2*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*a**2 + 2*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_758():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = -8*a*b*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + 2*a*sin(c + d*x)/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*sqrt(sec(c + d*x))) + 8*b*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a - 6*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_759():
    f = 1/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Symbol('a') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_760():
    f = (a + b*cos(c + d*x))**4*cos(c + d*x)**m
    F = 2*a*b**3*(m + 5)*sin(c + d*x)*cos(c + d*x)**(m + 2)/(d*(m + 3)*(m + 4)) - 4*a*b*(a**2*(m + 3) + b**2*(m + 2))*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*(m + 3)*sqrt(sin(c + d*x)**2)) + b**2*(a + b*cos(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 4)) + b**2*(a**2*(5*m + 22) + b**2*(m + 3))*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 2)*(m + 4)) - (a**4*(m**2 + 6*m + 8) + 6*a**2*b**2*(m**2 + 5*m + 4) + b**4*(m**2 + 4*m + 3))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*(m + 2)*(m + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_761():
    f = (a + b*cos(c + d*x))**3*cos(c + d*x)**m
    F = a*b**2*(2*m + 7)*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 2)*(m + 3)) - a*(a**2*(m + 2) + 3*b**2*(m + 1))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*(m + 2)*sqrt(sin(c + d*x)**2)) + b**2*(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 3)) - b*(3*a**2*(m + 3) + b**2*(m + 2))*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*(m + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_762():
    f = (a + b*cos(c + d*x))**2*cos(c + d*x)**m
    F = -2*a*b*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*sqrt(sin(c + d*x)**2)) + b**2*sin(c + d*x)*cos(c + d*x)**(m + 1)/(d*(m + 2)) - (a**2*(m + 2) + b**2*(m + 1))*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*(m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_763():
    f = (a + b*cos(c + d*x))*cos(c + d*x)**m
    F = -a*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + 1)*sqrt(sin(c + d*x)**2)) - b*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(c + d*x)**2)/(d*(m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_764():
    f = cos(c + d*x)**m/(a + b*cos(c + d*x))
    F = a*(cos(c + d*x)**2)**(sympy.S.Half - m/2)*sin(c + d*x)*cos(c + d*x)**(m - 1)*appellf1(sympy.S.Half, 1, sympy.S.Half - m/2, sympy.S(3)/2, -b**2*sin(c + d*x)**2/(a**2 - b**2), sin(c + d*x)**2)/(d*(a**2 - b**2)) - b*sin(c + d*x)*cos(c + d*x)**m*appellf1(sympy.S.Half, 1, -m/2, sympy.S(3)/2, -b**2*sin(c + d*x)**2/(a**2 - b**2), sin(c + d*x)**2)/(d*(a**2 - b**2)*(cos(c + d*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_765():
    f = cos(c + d*x)**m/(a + b*cos(c + d*x))**2
    F = a**2*(cos(c + d*x)**2)**(sympy.S.Half - m/2)*sin(c + d*x)*cos(c + d*x)**(m - 1)*appellf1(sympy.S.Half, 2, sympy.S.Half - m/2, sympy.S(3)/2, -b**2*sin(c + d*x)**2/(a**2 - b**2), sin(c + d*x)**2)/(d*(a**2 - b**2)**2) - 2*a*b*sin(c + d*x)*cos(c + d*x)**m*appellf1(sympy.S.Half, 2, -m/2, sympy.S(3)/2, -b**2*sin(c + d*x)**2/(a**2 - b**2), sin(c + d*x)**2)/(d*(a**2 - b**2)**2*(cos(c + d*x)**2)**(m/2)) + b**2*(cos(c + d*x)**2)**(-m/2 + sympy.S(-1)/2)*sin(c + d*x)*cos(c + d*x)**(m + 1)*appellf1(sympy.S.Half, 2, -m/2 + sympy.S(-1)/2, sympy.S(3)/2, -b**2*sin(c + d*x)**2/(a**2 - b**2), sin(c + d*x)**2)/(d*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_766():
    f = (a + b*cos(c + d*x))**3*sec(c + d*x)**m
    F = -a**2*b*(1 - 2*m)*sin(c + d*x)*sec(c + d*x)**(m - 2)/(d*(1 - m)*(2 - m)) - a**2*(a*sec(c + d*x) + b)*sin(c + d*x)*sec(c + d*x)**(m - 2)/(d*(1 - m)) - a*(a**2*(2 - m) + 3*b**2*(1 - m))*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/2 - m/2), (sympy.S(5)/2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 3)/(d*(1 - m)*(3 - m)*sqrt(sin(c + d*x)**2)) - b*(3*a**2*(3 - m) + b**2*(2 - m))*sin(c + d*x)*hyper((sympy.S.Half, 2 - m/2), (3 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 4)/(d*(2 - m)*(4 - m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_767():
    f = (a + b*cos(c + d*x))**2*sec(c + d*x)**m
    F = -a**2*sin(c + d*x)*sec(c + d*x)**(m - 1)/(d*(1 - m)) - 2*a*b*sin(c + d*x)*hyper((sympy.S.Half, 1 - m/2), (2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 2)/(d*(2 - m)*sqrt(sin(c + d*x)**2)) - (a**2*(2 - m) + b**2*(1 - m))*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(3)/2 - m/2), (sympy.S(5)/2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 3)/(d*(1 - m)*(3 - m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_768():
    f = (a + b*cos(c + d*x))*sec(c + d*x)**m
    F = -a*sin(c + d*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 1)/(d*(1 - m)*sqrt(sin(c + d*x)**2)) - b*sin(c + d*x)*hyper((sympy.S.Half, 1 - m/2), (2 - m/2,), cos(c + d*x)**2)*sec(c + d*x)**(m - 2)/(d*(2 - m)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_769():
    f = sqrt(1 - cos(x))/sqrt(a - cos(x))
    F = -2*atan(sin(x)/(sqrt(1 - cos(x))*sqrt(a - cos(x))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_770():
    f = sqrt((1 - cos(x))/(a - cos(x)))
    F = -2*sqrt((1 - cos(x))/(a - cos(x)))*sqrt(a - cos(x))*atan(sin(x)/(sqrt(1 - cos(x))*sqrt(a - cos(x))))/sqrt(1 - cos(x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_771():
    f = (B*cos(c + d*x) - B/2)*(a*cos(c + d*x) + a)
    F = B*a*sin(c + d*x)*cos(c + d*x)/(2*d) + B*a*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_772():
    f = (B*cos(c + d*x) - 4*B/5)*(a*cos(c + d*x) + a)**4
    F = B*(a*cos(c + d*x) + a)**4*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_773():
    f = (a*cos(c + d*x) + a)**n*(-B*n/(n + 1) + B*cos(c + d*x))
    F = B*(a*cos(c + d*x) + a)**n*sin(c + d*x)/(d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_774():
    f = (B*cos(c + d*x) - 3*B/2)/(a*cos(c + d*x) + a)**3
    F = -B*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_775():
    f = (B*cos(c + d*x) - 3*B/5)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*B*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_776():
    f = (B*cos(c + d*x) + B)/sqrt(a*cos(c + d*x) + a)
    F = 2*B*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_777():
    f = (B*cos(c + d*x) - 5*B/3)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -2*B*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_778():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(2)/3)
    F = 3*B*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)/(5*d) + 2*2**(sympy.S(1)/6)*(5*A + 2*B)*(a*cos(c + d*x) + a)**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(5*d*(cos(c + d*x) + 1)**(sympy.S(7)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_779():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*B*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)/(4*d) + 2**(sympy.S(5)/6)*(4*A + B)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(4*d*(cos(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_780():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(1)/3)
    F = 3*B*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(1)/6)*(2*A - B)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*(cos(c + d*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_781():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(2)/3)
    F = (3*A - 3*B)*sin(c + d*x)/(d*(a*cos(c + d*x) + a)**(sympy.S(2)/3)) - 2**(sympy.S(5)/6)*(A - 2*B)*(a*cos(c + d*x) + a)**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(3)/2,), sympy.S.Half - cos(c + d*x)/2)/(a*d*(cos(c + d*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_782():
    f = (B*cos(c + d*x) + B*b/a)/(a + b*cos(c + d*x))
    F = B*x/b - 2*B*sqrt(a - b)*sqrt(a + b)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_783():
    f = (a + b*cos(c + d*x))/(a*cos(c + d*x) + b)**2
    F = sin(c + d*x)/(d*(a*cos(c + d*x) + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_784():
    f = (cos(c + d*x) + 3)/(2 - cos(c + d*x))
    F = -x + 5*sqrt(3)*x/3 + 10*sqrt(3)*atan(sin(c + d*x)/(-cos(c + d*x) + sqrt(3) + 2))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_785():
    f = (B*a + B*b*cos(c + d*x))/sqrt(a + b*cos(c + d*x))
    F = 2*B*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_786():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(2)/3)
    F = sqrt(2)*B*(a + b)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-5)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(b*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(2)/3)*(A*b - B*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(b*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_787():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(1)/3)
    F = sqrt(2)*B*(a + b)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-4)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(b*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*(a + b*cos(c + d*x))**(sympy.S(1)/3)*(A*b - B*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(b*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_788():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(1)/3)
    F = sqrt(2)*B*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-2)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(b*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*(A*b - B*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(b*d*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_789():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(2)/3)
    F = sqrt(2)*B*(a + b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S(-1)/3, sympy.S.Half, sympy.S(3)/2, b*(1 - cos(c + d*x))/(a + b), sympy.S.Half - cos(c + d*x)/2)/(b*d*((a + b*cos(c + d*x))/(a + b))**(sympy.S(1)/3)*sqrt(cos(c + d*x) + 1)) + sqrt(2)*((a + b*cos(c + d*x))/(a + b))**(sympy.S(2)/3)*(A*b - B*a)*sin(c + d*x)*appellf1(sympy.S.Half, sympy.S.Half, sympy.S(2)/3, sympy.S(3)/2, sympy.S.Half - cos(c + d*x)/2, b*(1 - cos(c + d*x))/(a + b))/(b*d*(a + b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_790():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*cos(c + d*x)**2
    F = 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*A*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d) + 10*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_791():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*cos(c + d*x)
    F = 2*A*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_792():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))
    F = 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_793():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*sec(c + d*x)
    F = 2*A*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_794():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*sec(c + d*x)**2
    F = 2*A*b*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_795():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*sec(c + d*x)**3
    F = 2*A*b**2*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*A*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*b*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_796():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*A*b*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*b**2*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_797():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*cos(c + d*x)
    F = 6*A*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*A*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 10*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_798():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))
    F = 2*A*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 6*B*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_799():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*sec(c + d*x)
    F = 2*A*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_800():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**2
    F = 2*A*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) + 2*B*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_801():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**3
    F = 2*A*b**2*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*A*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_802():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*A*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*b**2*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_803():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*A*b**2*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*A*b*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*b**3*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b**2*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_804():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))
    F = 6*A*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*A*b*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 10*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_805():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))*sec(c + d*x)
    F = 2*A*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d) + 6*B*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*b*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_806():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**2
    F = 2*A*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_807():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**3
    F = 2*A*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) + 2*B*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_808():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**4
    F = 2*A*b**3*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*A*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x))) + 2*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_809():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**5
    F = 2*A*b**4*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*A*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*b**3*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_810():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))*sec(c + d*x)**6
    F = 2*A*b**5*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*A*b**3*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*A*b**2*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(cos(c + d*x))) + 2*B*b**4*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*b**3*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_811():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x))) + 2*A*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**2*d) + 10*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_812():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b*d) + 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_813():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_814():
    f = (A + B*cos(c + d*x))/sqrt(b*cos(c + d*x))
    F = 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_815():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/sqrt(b*cos(c + d*x))
    F = 2*A*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_816():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/sqrt(b*cos(c + d*x))
    F = 2*A*b*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x))) + 2*B*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_817():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/sqrt(b*cos(c + d*x))
    F = 2*A*b**2*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*A*sin(c + d*x)/(5*d*sqrt(b*cos(c + d*x))) - 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b*d*sqrt(cos(c + d*x))) + 2*B*b*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_818():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x))) + 2*A*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**3*d) + 10*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b**2*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_819():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x))) + 2*A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**2*d) + 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_820():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_821():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_822():
    f = (A + B*cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))) - 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_823():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x))) + 2*B*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_824():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*A*b*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*A*sin(c + d*x)/(5*b*d*sqrt(b*cos(c + d*x))) - 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**2*d*sqrt(cos(c + d*x))) + 2*B*sin(c + d*x)/(3*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_825():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**5/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x))) + 2*A*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**4*d) + 10*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*b**2*d*sqrt(b*cos(c + d*x))) + 10*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(21*b**3*d) + 2*B*(b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_826():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x))) + 2*A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**3*d) + 6*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x))) + 2*B*(b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_827():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_828():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b**2*d*sqrt(b*cos(c + d*x))) + 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_829():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))) - 2*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_830():
    f = (A + B*cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(3*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*A*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x))) + 2*B*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))) - 2*B*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(b**3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_831():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*A*sin(c + d*x)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*A*sin(c + d*x)/(5*b**2*d*sqrt(b*cos(c + d*x))) - 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**3*d*sqrt(cos(c + d*x))) + 2*B*sin(c + d*x)/(3*b*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_832():
    f = (A + B*cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(7)/2)
    F = 2*A*sin(c + d*x)/(5*b*d*(b*cos(c + d*x))**(sympy.S(5)/2)) + 6*A*sin(c + d*x)/(5*b**3*d*sqrt(b*cos(c + d*x))) - 6*A*sqrt(b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(5*b**4*d*sqrt(cos(c + d*x))) + 2*B*sin(c + d*x)/(3*b**2*d*(b*cos(c + d*x))**(sympy.S(3)/2)) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*b**3*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_833():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = -A*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 3*B*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*B*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_834():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = A*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + A*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d) - B*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_835():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))*sqrt(cos(c + d*x))
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + B*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_836():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))/sqrt(cos(c + d*x))
    F = A*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_837():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = A*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + B*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_838():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_839():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + A*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x))) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_840():
    f = sqrt(b*cos(c + d*x))*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(9)/2)
    F = A*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + A*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_841():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = -A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 3*B*b*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_842():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))*sqrt(cos(c + d*x))
    F = A*b*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d) - B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_843():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))/sqrt(cos(c + d*x))
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + B*b*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_844():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = A*b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_845():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = A*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + B*b*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_846():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_847():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(9)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + A*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x))) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_848():
    f = (b*cos(c + d*x))**(sympy.S(3)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(11)/2)
    F = A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + A*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*b*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*b*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_849():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))*sqrt(cos(c + d*x))
    F = -A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 3*B*b**2*x*sqrt(b*cos(c + d*x))/(8*sqrt(cos(c + d*x))) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + 3*B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_850():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))/sqrt(cos(c + d*x))
    F = A*b**2*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d) - B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*sqrt(cos(c + d*x))) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_851():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + B*b**2*x*sqrt(b*cos(c + d*x))/(2*sqrt(cos(c + d*x))) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_852():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = A*b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x)) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_853():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x))) + B*b**2*x*sqrt(b*cos(c + d*x))/sqrt(cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_854():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(9)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_855():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(11)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + A*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x))) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_856():
    f = (b*cos(c + d*x))**(sympy.S(5)/2)*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(13)/2)
    F = A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)**3/(3*d*cos(c + d*x)**(sympy.S(7)/2)) + A*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(d*cos(c + d*x)**(sympy.S(3)/2)) + B*b**2*sqrt(b*cos(c + d*x))*sin(c + d*x)/(2*d*cos(c + d*x)**(sympy.S(5)/2)) + B*b**2*sqrt(b*cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_857():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/sqrt(b*cos(c + d*x))
    F = A*x*sqrt(cos(c + d*x))/(2*sqrt(b*cos(c + d*x))) + A*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(b*cos(c + d*x))) - B*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_858():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/sqrt(b*cos(c + d*x))
    F = A*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/(2*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_859():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    F = A*x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x)) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_860():
    f = (A + B*cos(c + d*x))/(sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/sqrt(b*cos(c + d*x))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_861():
    f = (A + B*cos(c + d*x))/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_862():
    f = (A + B*cos(c + d*x))/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)/(2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_863():
    f = (A + B*cos(c + d*x))/(sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = A*sin(c + d*x)**3/(3*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + A*sin(c + d*x)/(d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sin(c + d*x)/(2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_864():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*x*sqrt(cos(c + d*x))/(2*b*sqrt(b*cos(c + d*x))) + A*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b*d*sqrt(b*cos(c + d*x))) - B*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_865():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/(2*b*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_866():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_867():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(3)/2)
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b*d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/(b*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_868():
    f = (A + B*cos(c + d*x))/((b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_869():
    f = (A + B*cos(c + d*x))/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)/(2*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_870():
    f = (A + B*cos(c + d*x))/((b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = A*sin(c + d*x)**3/(3*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + A*sin(c + d*x)/(b*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sin(c + d*x)/(2*b*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_871():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(9)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*x*sqrt(cos(c + d*x))/(2*b**2*sqrt(b*cos(c + d*x))) + A*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b**2*d*sqrt(b*cos(c + d*x))) - B*sin(c + d*x)**3*sqrt(cos(c + d*x))/(3*b**2*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_872():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/(2*b**2*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_873():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_874():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x))) + B*x*sqrt(cos(c + d*x))/(b**2*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_875():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(5)/2)
    F = A*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_876():
    f = (A + B*cos(c + d*x))/((b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = A*sin(c + d*x)/(2*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + A*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b**2*d*sqrt(b*cos(c + d*x))) + B*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_877():
    f = (A + B*cos(c + d*x))/((b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = A*sin(c + d*x)**3/(3*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)) + A*sin(c + d*x)/(b**2*d*sqrt(b*cos(c + d*x))*sqrt(cos(c + d*x))) + B*sin(c + d*x)/(2*b**2*d*sqrt(b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)) + B*sqrt(cos(c + d*x))*atanh(sin(c + d*x))/(2*b**2*d*sqrt(b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_878():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x))*cos(c + d*x)**2
    F = -3*A*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**3*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(13)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(13)/6), (sympy.S(19)/6,), cos(c + d*x)**2)/(13*b**4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_879():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x))*cos(c + d*x)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b**2*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_880():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x))
    F = -3*A*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*b*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_881():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x))*sec(c + d*x)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_882():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x))*sec(c + d*x)**2
    F = 3*A*b*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_883():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x))*sec(c + d*x)**3
    F = 3*A*b**2*sin(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(1)/6,), cos(c + d*x)**2)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2)) + 3*B*b*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_884():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x))*cos(c + d*x)**2
    F = -3*A*(b*cos(c + d*x))**(sympy.S(13)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(13)/6), (sympy.S(19)/6,), cos(c + d*x)**2)/(13*b**3*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(16)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(8)/3), (sympy.S(11)/3,), cos(c + d*x)**2)/(16*b**4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_885():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x))*cos(c + d*x)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**2*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(13)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(13)/6), (sympy.S(19)/6,), cos(c + d*x)**2)/(13*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_886():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x))
    F = -3*A*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_887():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x))*sec(c + d*x)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_888():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x))*sec(c + d*x)**2
    F = -3*A*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_889():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x))*sec(c + d*x)**3
    F = 3*A*b**2*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2)) - 3*B*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_890():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b**3*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(10)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/3), (sympy.S(8)/3,), cos(c + d*x)**2)/(10*b**4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_891():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*b**2*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(7)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(7)/6), (sympy.S(13)/6,), cos(c + d*x)**2)/(7*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_892():
    f = (A + B*cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(4)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(c + d*x)**2)/(4*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_893():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*A*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*hyper((sympy.S(1)/6, sympy.S.Half), (sympy.S(7)/6,), cos(c + d*x)**2)/(b*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_894():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*A*b*sin(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(1)/6,), cos(c + d*x)**2)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2)) + 3*B*sin(c + d*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(c + d*x)**2)/(2*d*(b*cos(c + d*x))**(sympy.S(2)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_895():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = 3*A*b**2*sin(c + d*x)*hyper((sympy.S(-4)/3, sympy.S.Half), (sympy.S(-1)/3,), cos(c + d*x)**2)/(8*d*(b*cos(c + d*x))**(sympy.S(8)/3)*sqrt(sin(c + d*x)**2)) + 3*B*b*sin(c + d*x)*hyper((sympy.S(-5)/6, sympy.S.Half), (sympy.S(1)/6,), cos(c + d*x)**2)/(5*d*(b*cos(c + d*x))**(sympy.S(5)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_896():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**3*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(8)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(c + d*x)**2)/(8*b**4*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_897():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = -3*A*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*b**2*d*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(5)/3)*sin(c + d*x)*hyper((sympy.S.Half, sympy.S(5)/6), (sympy.S(11)/6,), cos(c + d*x)**2)/(5*b**3*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_898():
    f = (A + B*cos(c + d*x))/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(c + d*x)**2)/(2*b**2*d*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_899():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2)) + 3*B*sin(c + d*x)*hyper((sympy.S(-1)/6, sympy.S.Half), (sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_900():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*b*sin(c + d*x)*hyper((sympy.S(-7)/6, sympy.S.Half), (sympy.S(-1)/6,), cos(c + d*x)**2)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2)) + 3*B*sin(c + d*x)*hyper((sympy.S(-2)/3, sympy.S.Half), (sympy.S(1)/3,), cos(c + d*x)**2)/(4*d*(b*cos(c + d*x))**(sympy.S(4)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_901():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*b**2*sin(c + d*x)*hyper((sympy.S(-5)/3, sympy.S.Half), (sympy.S(-2)/3,), cos(c + d*x)**2)/(10*d*(b*cos(c + d*x))**(sympy.S(10)/3)*sqrt(sin(c + d*x)**2)) + 3*B*b*sin(c + d*x)*hyper((sympy.S(-7)/6, sympy.S.Half), (sympy.S(-1)/6,), cos(c + d*x)**2)/(7*d*(b*cos(c + d*x))**(sympy.S(7)/3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_902():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*cos(c + d*x)**m
    F = -A*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + n/2 + sympy.S.Half), (m/2 + n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(d*(m + n + 1)*sqrt(sin(c + d*x)**2)) - B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + n/2 + 1), (m/2 + n/2 + 2,), cos(c + d*x)**2)/(d*(m + n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_903():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*cos(c + d*x)**2
    F = -A*(b*cos(c + d*x))**(n + 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*sqrt(sin(c + d*x)**2)) - B*(b*cos(c + d*x))**(n + 4)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 2), (n/2 + 3,), cos(c + d*x)**2)/(b**4*d*(n + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_904():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*cos(c + d*x)
    F = -A*(b*cos(c + d*x))**(n + 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*sqrt(sin(c + d*x)**2)) - B*(b*cos(c + d*x))**(n + 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(3)/2), (n/2 + sympy.S(5)/2,), cos(c + d*x)**2)/(b**3*d*(n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_905():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))
    F = -A*(b*cos(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2)) - B*(b*cos(c + d*x))**(n + 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + 1), (n/2 + 2,), cos(c + d*x)**2)/(b**2*d*(n + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_906():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*sec(c + d*x)
    F = -A*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2)) - B*(b*cos(c + d*x))**(n + 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S.Half), (n/2 + sympy.S(3)/2,), cos(c + d*x)**2)/(b*d*(n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_907():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*sec(c + d*x)**2
    F = A*b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2)) - B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2), (n/2 + 1,), cos(c + d*x)**2)/(d*n*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_908():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*sec(c + d*x)**3
    F = A*b**2*(b*cos(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 - 1), (n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2)) + B*b*(b*cos(c + d*x))**(n - 1)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/2), (n/2 + sympy.S.Half,), cos(c + d*x)**2)/(d*(1 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_909():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*sec(c + d*x)**4
    F = A*b**3*(b*cos(c + d*x))**(n - 3)*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/2), (n/2 + sympy.S(-1)/2,), cos(c + d*x)**2)/(d*(3 - n)*sqrt(sin(c + d*x)**2)) + B*b**2*(b*cos(c + d*x))**(n - 2)*sin(c + d*x)*hyper((sympy.S.Half, n/2 - 1), (n/2,), cos(c + d*x)**2)/(d*(2 - n)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_910():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = -2*A*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)*hyper((sympy.S.Half, n/2 + sympy.S(7)/4), (n/2 + sympy.S(11)/4,), cos(c + d*x)**2)/(d*(2*n + 7)*sqrt(sin(c + d*x)**2)) - 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)*hyper((sympy.S.Half, n/2 + sympy.S(9)/4), (n/2 + sympy.S(13)/4,), cos(c + d*x)**2)/(d*(2*n + 9)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_911():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = -2*A*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*sqrt(sin(c + d*x)**2)) - 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)*hyper((sympy.S.Half, n/2 + sympy.S(7)/4), (n/2 + sympy.S(11)/4,), cos(c + d*x)**2)/(d*(2*n + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_912():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))*sqrt(cos(c + d*x))
    F = -2*A*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2)) - 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)*hyper((sympy.S.Half, n/2 + sympy.S(5)/4), (n/2 + sympy.S(9)/4,), cos(c + d*x)**2)/(d*(2*n + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_913():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))/sqrt(cos(c + d*x))
    F = -2*A*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*sqrt(sin(c + d*x)**2)) - 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)*hyper((sympy.S.Half, n/2 + sympy.S(3)/4), (n/2 + sympy.S(7)/4,), cos(c + d*x)**2)/(d*(2*n + 3)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_914():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x))) - 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*sqrt(cos(c + d*x))*hyper((sympy.S.Half, n/2 + sympy.S(1)/4), (n/2 + sympy.S(5)/4,), cos(c + d*x)**2)/(d*(2*n + 1)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_915():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2)) + 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-1)/4), (n/2 + sympy.S(3)/4,), cos(c + d*x)**2)/(d*(1 - 2*n)*sqrt(sin(c + d*x)**2)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_916():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-5)/4), (n/2 + sympy.S(-1)/4,), cos(c + d*x)**2)/(d*(5 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2)) + 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-3)/4), (n/2 + sympy.S(1)/4,), cos(c + d*x)**2)/(d*(3 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_917():
    f = (b*cos(c + d*x))**n*(A + B*cos(c + d*x))/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-7)/4), (n/2 + sympy.S(-3)/4,), cos(c + d*x)**2)/(d*(7 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(7)/2)) + 2*B*(b*cos(c + d*x))**n*sin(c + d*x)*hyper((sympy.S.Half, n/2 + sympy.S(-5)/4), (n/2 + sympy.S(-1)/4,), cos(c + d*x)**2)/(d*(5 - 2*n)*sqrt(sin(c + d*x)**2)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_918():
    f = (b*cos(c + d*x))**(sympy.S(4)/3)*(A + B*cos(c + d*x))*cos(c + d*x)**m
    F = -3*A*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(3*m + 7)*sqrt(sin(c + d*x)**2)) - 3*B*b*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 3)*hyper((sympy.S.Half, m/2 + sympy.S(5)/3), (m/2 + sympy.S(8)/3,), cos(c + d*x)**2)/(d*(3*m + 10)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_919():
    f = (b*cos(c + d*x))**(sympy.S(2)/3)*(A + B*cos(c + d*x))*cos(c + d*x)**m
    F = -3*A*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(d*(3*m + 5)*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(2)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(4)/3), (m/2 + sympy.S(7)/3,), cos(c + d*x)**2)/(d*(3*m + 8)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_920():
    f = (b*cos(c + d*x))**(sympy.S(1)/3)*(A + B*cos(c + d*x))*cos(c + d*x)**m
    F = -3*A*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(2)/3), (m/2 + sympy.S(5)/3,), cos(c + d*x)**2)/(d*(3*m + 4)*sqrt(sin(c + d*x)**2)) - 3*B*(b*cos(c + d*x))**(sympy.S(1)/3)*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(7)/6), (m/2 + sympy.S(13)/6,), cos(c + d*x)**2)/(d*(3*m + 7)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_921():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(1)/3)
    F = -3*A*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/3), (m/2 + sympy.S(4)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(5)/6), (m/2 + sympy.S(11)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 5)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_922():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(2)/3)
    F = -3*A*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/6), (m/2 + sympy.S(7)/6,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 1)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*cos(c + d*x)**(m + 2)*hyper((sympy.S.Half, m/2 + sympy.S(2)/3), (m/2 + sympy.S(5)/3,), cos(c + d*x)**2)/(d*(b*cos(c + d*x))**(sympy.S(2)/3)*(3*m + 4)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_2_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_923():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**m/(b*cos(c + d*x))**(sympy.S(4)/3)
    F = 3*A*sin(c + d*x)*cos(c + d*x)**m*hyper((sympy.S.Half, m/2 + sympy.S(-1)/6), (m/2 + sympy.S(5)/6,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(1 - 3*m)*sqrt(sin(c + d*x)**2)) - 3*B*sin(c + d*x)*cos(c + d*x)**(m + 1)*hyper((sympy.S.Half, m/2 + sympy.S(1)/3), (m/2 + sympy.S(4)/3,), cos(c + d*x)**2)/(b*d*(b*cos(c + d*x))**(sympy.S(1)/3)*(3*m + 2)*sqrt(sin(c + d*x)**2))
    assert integrate(f, x) == F

