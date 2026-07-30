"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.2.3 (g sec)^p (a+b sec)^m (c+d sec)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, g, m, n, p = symbols('a b c d e f g m n p')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_1():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**4*sec(e + f*x)
    F = a*c**4*tan(e + f*x)**5/(5*f) + 4*a*c**4*tan(e + f*x)**3/(3*f) - 3*a*c**4*tan(e + f*x)*sec(e + f*x)**3/(4*f) - a*c**4*tan(e + f*x)*sec(e + f*x)/(8*f) + 7*a*c**4*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_2():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**3*sec(e + f*x)
    F = 2*a*c**3*tan(e + f*x)**3/(3*f) - a*c**3*tan(e + f*x)*sec(e + f*x)**3/(4*f) - 3*a*c**3*tan(e + f*x)*sec(e + f*x)/(8*f) + 5*a*c**3*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_3():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**2*sec(e + f*x)
    F = a*c**2*tan(e + f*x)**3/(3*f) - a*c**2*tan(e + f*x)*sec(e + f*x)/(2*f) + a*c**2*atanh(sin(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_4():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -a*c*tan(e + f*x)*sec(e + f*x)/(2*f) + a*c*atanh(sin(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_5():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)
    F = -2*a*tan(e + f*x)/(f*(-c*sec(e + f*x) + c)) - a*atanh(sin(e + f*x))/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_6():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**2
    F = -(a*sec(e + f*x) + a)*tan(e + f*x)/(3*f*(-c*sec(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_7():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**3
    F = -(a*sec(e + f*x) + a)*tan(e + f*x)/(5*f*(-c*sec(e + f*x) + c)**3) - (a*sec(e + f*x) + a)*tan(e + f*x)/(15*c*f*(-c*sec(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_8():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**4
    F = -(a*sec(e + f*x) + a)*tan(e + f*x)/(7*f*(-c*sec(e + f*x) + c)**4) - (2*a*sec(e + f*x) + 2*a)*tan(e + f*x)/(105*f*(-c**2*sec(e + f*x) + c**2)**2) - (2*a*sec(e + f*x) + 2*a)*tan(e + f*x)/(35*c*f*(-c*sec(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_9():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**5
    F = -(a*sec(e + f*x) + a)*tan(e + f*x)/(9*f*(-c*sec(e + f*x) + c)**5) - (a*sec(e + f*x) + a)*tan(e + f*x)/(21*c*f*(-c*sec(e + f*x) + c)**4) - (2*a*sec(e + f*x) + 2*a)*tan(e + f*x)/(315*c*f*(-c**2*sec(e + f*x) + c**2)**2) - (2*a*sec(e + f*x) + 2*a)*tan(e + f*x)/(105*c**2*f*(-c*sec(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_10():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**5*sec(e + f*x)
    F = -a**2*c**5*tan(e + f*x)**7/(7*f) - 4*a**2*c**5*tan(e + f*x)**5/(5*f) + a**2*c**5*tan(e + f*x)**3*sec(e + f*x)**3/(2*f) + a**2*c**5*tan(e + f*x)**3*sec(e + f*x)/(4*f) - 3*a**2*c**5*tan(e + f*x)*sec(e + f*x)**3/(8*f) - 3*a**2*c**5*tan(e + f*x)*sec(e + f*x)/(16*f) + 9*a**2*c**5*atanh(sin(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_11():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**4*sec(e + f*x)
    F = -2*a**2*c**4*tan(e + f*x)**5/(5*f) + a**2*c**4*tan(e + f*x)**3*sec(e + f*x)**3/(6*f) + a**2*c**4*tan(e + f*x)**3*sec(e + f*x)/(4*f) - a**2*c**4*tan(e + f*x)*sec(e + f*x)**3/(8*f) - 5*a**2*c**4*tan(e + f*x)*sec(e + f*x)/(16*f) + 7*a**2*c**4*atanh(sin(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_12():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**3*sec(e + f*x)
    F = -a**2*c**3*tan(e + f*x)**5/(5*f) + a**2*c**3*tan(e + f*x)**3*sec(e + f*x)/(4*f) - 3*a**2*c**3*tan(e + f*x)*sec(e + f*x)/(8*f) + 3*a**2*c**3*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_13():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**2*sec(e + f*x)
    F = a**2*c**2*tan(e + f*x)**3*sec(e + f*x)/(4*f) - 3*a**2*c**2*tan(e + f*x)*sec(e + f*x)/(8*f) + 3*a**2*c**2*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_14():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -a**2*c*tan(e + f*x)**3/(3*f) - a**2*c*tan(e + f*x)*sec(e + f*x)/(2*f) + a**2*c*atanh(sin(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_15():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)/(-c*sec(e + f*x) + c)
    F = -3*a**2*tan(e + f*x)/(c*f) - 3*a**2*atanh(sin(e + f*x))/(c*f) - (2*a**2*sec(e + f*x) + 2*a**2)*tan(e + f*x)/(f*(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_16():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)/(-c*sec(e + f*x) + c)**2
    F = 2*a**2*tan(e + f*x)/(f*(-c**2*sec(e + f*x) + c**2)) + a**2*atanh(sin(e + f*x))/(c**2*f) - (2*a**2*sec(e + f*x) + 2*a**2)*tan(e + f*x)/(3*f*(-c*sec(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_17():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)/(-c*sec(e + f*x) + c)**3
    F = -(a*sec(e + f*x) + a)**2*tan(e + f*x)/(5*f*(-c*sec(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_18():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)/(-c*sec(e + f*x) + c)**4
    F = -(a*sec(e + f*x) + a)**2*tan(e + f*x)/(7*f*(-c*sec(e + f*x) + c)**4) - (a*sec(e + f*x) + a)**2*tan(e + f*x)/(35*c*f*(-c*sec(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_19():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)/(-c*sec(e + f*x) + c)**5
    F = -(a*sec(e + f*x) + a)**2*tan(e + f*x)/(9*f*(-c*sec(e + f*x) + c)**5) - 2*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(63*c*f*(-c*sec(e + f*x) + c)**4) - 2*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(315*c**2*f*(-c*sec(e + f*x) + c)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_20():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)/(-c*sec(e + f*x) + c)**6
    F = -2*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(1155*f*(-c**2*sec(e + f*x) + c**2)**3) - (a*sec(e + f*x) + a)**2*tan(e + f*x)/(11*f*(-c*sec(e + f*x) + c)**6) - (a*sec(e + f*x) + a)**2*tan(e + f*x)/(33*c*f*(-c*sec(e + f*x) + c)**5) - 2*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(231*c**2*f*(-c*sec(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_21():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**6*sec(e + f*x)
    F = a**3*c**6*tan(e + f*x)**9/(9*f) + 4*a**3*c**6*tan(e + f*x)**7/(7*f) - 3*a**3*c**6*tan(e + f*x)**5*sec(e + f*x)**3/(8*f) - a**3*c**6*tan(e + f*x)**5*sec(e + f*x)/(6*f) + 5*a**3*c**6*tan(e + f*x)**3*sec(e + f*x)**3/(16*f) + 5*a**3*c**6*tan(e + f*x)**3*sec(e + f*x)/(24*f) - 15*a**3*c**6*tan(e + f*x)*sec(e + f*x)**3/(64*f) - 25*a**3*c**6*tan(e + f*x)*sec(e + f*x)/(128*f) + 55*a**3*c**6*atanh(sin(e + f*x))/(128*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_22():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**5*sec(e + f*x)
    F = 2*a**3*c**5*tan(e + f*x)**7/(7*f) - a**3*c**5*tan(e + f*x)**5*sec(e + f*x)**3/(8*f) - a**3*c**5*tan(e + f*x)**5*sec(e + f*x)/(6*f) + 5*a**3*c**5*tan(e + f*x)**3*sec(e + f*x)**3/(48*f) + 5*a**3*c**5*tan(e + f*x)**3*sec(e + f*x)/(24*f) - 5*a**3*c**5*tan(e + f*x)*sec(e + f*x)**3/(64*f) - 35*a**3*c**5*tan(e + f*x)*sec(e + f*x)/(128*f) + 45*a**3*c**5*atanh(sin(e + f*x))/(128*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_23():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**4*sec(e + f*x)
    F = a**3*c**4*tan(e + f*x)**7/(7*f) - a**3*c**4*tan(e + f*x)**5*sec(e + f*x)/(6*f) + 5*a**3*c**4*tan(e + f*x)**3*sec(e + f*x)/(24*f) - 5*a**3*c**4*tan(e + f*x)*sec(e + f*x)/(16*f) + 5*a**3*c**4*atanh(sin(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_24():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**3*sec(e + f*x)
    F = -a**3*c**3*tan(e + f*x)**5*sec(e + f*x)/(6*f) + 5*a**3*c**3*tan(e + f*x)**3*sec(e + f*x)/(24*f) - 5*a**3*c**3*tan(e + f*x)*sec(e + f*x)/(16*f) + 5*a**3*c**3*atanh(sin(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_25():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**2*sec(e + f*x)
    F = a**3*c**2*tan(e + f*x)**5/(5*f) + a**3*c**2*tan(e + f*x)**3*sec(e + f*x)/(4*f) - 3*a**3*c**2*tan(e + f*x)*sec(e + f*x)/(8*f) + 3*a**3*c**2*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_26():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -2*a**3*c*tan(e + f*x)**3/(3*f) - a**3*c*tan(e + f*x)*sec(e + f*x)**3/(4*f) - 3*a**3*c*tan(e + f*x)*sec(e + f*x)/(8*f) + 5*a**3*c*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_27():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)
    F = -5*a**3*tan(e + f*x)*sec(e + f*x)/(2*c*f) - 10*a**3*tan(e + f*x)/(c*f) - 15*a**3*atanh(sin(e + f*x))/(2*c*f) - 2*a*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(f*(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_28():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**2
    F = 5*a**3*tan(e + f*x)/(c**2*f) + 5*a**3*atanh(sin(e + f*x))/(c**2*f) - 2*a*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(3*f*(-c*sec(e + f*x) + c)**2) + (10*a**3*sec(e + f*x) + 10*a**3)*tan(e + f*x)/(3*f*(-c**2*sec(e + f*x) + c**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_29():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**3
    F = -2*a**3*tan(e + f*x)/(f*(-c**3*sec(e + f*x) + c**3)) - a**3*atanh(sin(e + f*x))/(c**3*f) - 2*a*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(5*f*(-c*sec(e + f*x) + c)**3) + (2*a**3*sec(e + f*x) + 2*a**3)*tan(e + f*x)/(3*c*f*(-c*sec(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_30():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**4
    F = -(a*sec(e + f*x) + a)**3*tan(e + f*x)/(7*f*(-c*sec(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_31():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**5
    F = -(a*sec(e + f*x) + a)**3*tan(e + f*x)/(9*f*(-c*sec(e + f*x) + c)**5) - (a*sec(e + f*x) + a)**3*tan(e + f*x)/(63*c*f*(-c*sec(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_32():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**6
    F = -(a*sec(e + f*x) + a)**3*tan(e + f*x)/(11*f*(-c*sec(e + f*x) + c)**6) - 2*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(99*c*f*(-c*sec(e + f*x) + c)**5) - 2*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(693*c**2*f*(-c*sec(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_33():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**7
    F = -(a*sec(e + f*x) + a)**3*tan(e + f*x)/(13*f*(-c*sec(e + f*x) + c)**7) - 3*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(143*c*f*(-c*sec(e + f*x) + c)**6) - 2*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(429*c**2*f*(-c*sec(e + f*x) + c)**5) - 2*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(3003*c**3*f*(-c*sec(e + f*x) + c)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_34():
    f = (-c*sec(e + f*x) + c)**4*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = 2*c*(-c*sec(e + f*x) + c)**3*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) + 7*c**4*tan(e + f*x)**3/(3*a*f) - 21*c**4*tan(e + f*x)*sec(e + f*x)/(2*a*f) + 28*c**4*tan(e + f*x)/(a*f) - 35*c**4*atanh(sin(e + f*x))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_35():
    f = (-c*sec(e + f*x) + c)**3*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = 2*c*(-c*sec(e + f*x) + c)**2*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) - 5*c**3*tan(e + f*x)*sec(e + f*x)/(2*a*f) + 10*c**3*tan(e + f*x)/(a*f) - 15*c**3*atanh(sin(e + f*x))/(2*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_36():
    f = (-c*sec(e + f*x) + c)**2*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = (-2*c**2*sec(e + f*x) + 2*c**2)*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) + 3*c**2*tan(e + f*x)/(a*f) - 3*c**2*atanh(sin(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_37():
    f = (-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = 2*c*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) - c*atanh(sin(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_38():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c))
    F = csc(e + f*x)/(a*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_39():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**2)
    F = -cot(e + f*x)**3/(3*a*c**2*f) - csc(e + f*x)**3/(3*a*c**2*f) + csc(e + f*x)/(a*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_40():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**3)
    F = 2*cot(e + f*x)**5/(5*a*c**3*f) + 2*csc(e + f*x)**5/(5*a*c**3*f) - csc(e + f*x)**3/(a*c**3*f) + csc(e + f*x)/(a*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_41():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**4)
    F = -4*cot(e + f*x)**7/(7*a*c**4*f) - cot(e + f*x)**5/(5*a*c**4*f) - 4*csc(e + f*x)**7/(7*a*c**4*f) + 9*csc(e + f*x)**5/(5*a*c**4*f) - 2*csc(e + f*x)**3/(a*c**4*f) + csc(e + f*x)/(a*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_42():
    f = (-c*sec(e + f*x) + c)**5*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = -6*c**2*(-c*sec(e + f*x) + c)**3*tan(e + f*x)/(f*(a**2*sec(e + f*x) + a**2)) + 2*c*(-c*sec(e + f*x) + c)**4*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - 7*c**5*tan(e + f*x)**3/(a**2*f) + 63*c**5*tan(e + f*x)*sec(e + f*x)/(2*a**2*f) - 84*c**5*tan(e + f*x)/(a**2*f) + 105*c**5*atanh(sin(e + f*x))/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_43():
    f = (-c*sec(e + f*x) + c)**4*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = 2*c*(-c*sec(e + f*x) + c)**3*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - 14*(-c**2*sec(e + f*x) + c**2)**2*tan(e + f*x)/(3*f*(a**2*sec(e + f*x) + a**2)) + 35*c**4*tan(e + f*x)*sec(e + f*x)/(6*a**2*f) - 70*c**4*tan(e + f*x)/(3*a**2*f) + 35*c**4*atanh(sin(e + f*x))/(2*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_44():
    f = (-c*sec(e + f*x) + c)**3*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = 2*c*(-c*sec(e + f*x) + c)**2*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - (-10*c**3*sec(e + f*x) + 10*c**3)*tan(e + f*x)/(3*f*(a**2*sec(e + f*x) + a**2)) - 5*c**3*tan(e + f*x)/(a**2*f) + 5*c**3*atanh(sin(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_45():
    f = (-c*sec(e + f*x) + c)**2*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = -2*c**2*tan(e + f*x)/(f*(a**2*sec(e + f*x) + a**2)) + (-2*c**2*sec(e + f*x) + 2*c**2)*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) + c**2*atanh(sin(e + f*x))/(a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_46():
    f = (-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = (-c*sec(e + f*x) + c)*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_47():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c))
    F = cot(e + f*x)**3/(3*a**2*c*f) - csc(e + f*x)**3/(3*a**2*c*f) + csc(e + f*x)/(a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_48():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**2)
    F = -csc(e + f*x)**3/(3*a**2*c**2*f) + csc(e + f*x)/(a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_49():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**3)
    F = cot(e + f*x)**5/(5*a**2*c**3*f) + csc(e + f*x)**5/(5*a**2*c**3*f) - 2*csc(e + f*x)**3/(3*a**2*c**3*f) + csc(e + f*x)/(a**2*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_50():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**4)
    F = -2*cot(e + f*x)**7/(7*a**2*c**4*f) - 2*csc(e + f*x)**7/(7*a**2*c**4*f) + csc(e + f*x)**5/(a**2*c**4*f) - 4*csc(e + f*x)**3/(3*a**2*c**4*f) + csc(e + f*x)/(a**2*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_51():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**5)
    F = 4*cot(e + f*x)**9/(9*a**2*c**5*f) + cot(e + f*x)**7/(7*a**2*c**5*f) + 4*csc(e + f*x)**9/(9*a**2*c**5*f) - 13*csc(e + f*x)**7/(7*a**2*c**5*f) + 3*csc(e + f*x)**5/(a**2*c**5*f) - 7*csc(e + f*x)**3/(3*a**2*c**5*f) + csc(e + f*x)/(a**2*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_52():
    f = (-c*sec(e + f*x) + c)**6*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 2*c*(-c*sec(e + f*x) + c)**5*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) + 66*(-c**2*sec(e + f*x) + c**2)**3*tan(e + f*x)/(5*f*(a**3*sec(e + f*x) + a**3)) - 22*c**2*(-c*sec(e + f*x) + c)**4*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2) + 77*c**6*tan(e + f*x)**3/(5*a**3*f) - 693*c**6*tan(e + f*x)*sec(e + f*x)/(10*a**3*f) + 924*c**6*tan(e + f*x)/(5*a**3*f) - 231*c**6*atanh(sin(e + f*x))/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_53():
    f = (-c*sec(e + f*x) + c)**5*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 42*c*(-c**2*sec(e + f*x) + c**2)**2*tan(e + f*x)/(5*f*(a**3*sec(e + f*x) + a**3)) + 2*c*(-c*sec(e + f*x) + c)**4*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) - 6*c**2*(-c*sec(e + f*x) + c)**3*tan(e + f*x)/(5*a*f*(a*sec(e + f*x) + a)**2) - 21*c**5*tan(e + f*x)*sec(e + f*x)/(2*a**3*f) + 42*c**5*tan(e + f*x)/(a**3*f) - 63*c**5*atanh(sin(e + f*x))/(2*a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_54():
    f = (-c*sec(e + f*x) + c)**4*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 2*c*(-c*sec(e + f*x) + c)**3*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) + (-14*c**4*sec(e + f*x) + 14*c**4)*tan(e + f*x)/(3*f*(a**3*sec(e + f*x) + a**3)) - 14*(-c**2*sec(e + f*x) + c**2)**2*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2) + 7*c**4*tan(e + f*x)/(a**3*f) - 7*c**4*atanh(sin(e + f*x))/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_55():
    f = (-c*sec(e + f*x) + c)**3*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 2*c**3*tan(e + f*x)/(f*(a**3*sec(e + f*x) + a**3)) + 2*c*(-c*sec(e + f*x) + c)**2*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) - (-2*c**3*sec(e + f*x) + 2*c**3)*tan(e + f*x)/(3*a*f*(a*sec(e + f*x) + a)**2) - c**3*atanh(sin(e + f*x))/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_56():
    f = (-c*sec(e + f*x) + c)**2*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = (-c*sec(e + f*x) + c)**2*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_57():
    f = (-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = (-c*sec(e + f*x) + c)*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) + (-c*sec(e + f*x) + c)*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_58():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c))
    F = -2*cot(e + f*x)**5/(5*a**3*c*f) + 2*csc(e + f*x)**5/(5*a**3*c*f) - csc(e + f*x)**3/(a**3*c*f) + csc(e + f*x)/(a**3*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_59():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**2)
    F = -cot(e + f*x)**5/(5*a**3*c**2*f) + csc(e + f*x)**5/(5*a**3*c**2*f) - 2*csc(e + f*x)**3/(3*a**3*c**2*f) + csc(e + f*x)/(a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_60():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**3)
    F = csc(e + f*x)**5/(5*a**3*c**3*f) - 2*csc(e + f*x)**3/(3*a**3*c**3*f) + csc(e + f*x)/(a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_61():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**4)
    F = -cot(e + f*x)**7/(7*a**3*c**4*f) - csc(e + f*x)**7/(7*a**3*c**4*f) + 3*csc(e + f*x)**5/(5*a**3*c**4*f) - csc(e + f*x)**3/(a**3*c**4*f) + csc(e + f*x)/(a**3*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_62():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**5)
    F = 2*cot(e + f*x)**9/(9*a**3*c**5*f) + 2*csc(e + f*x)**9/(9*a**3*c**5*f) - csc(e + f*x)**7/(a**3*c**5*f) + 9*csc(e + f*x)**5/(5*a**3*c**5*f) - 5*csc(e + f*x)**3/(3*a**3*c**5*f) + csc(e + f*x)/(a**3*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_63():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**6)
    F = -4*cot(e + f*x)**11/(11*a**3*c**6*f) - cot(e + f*x)**9/(9*a**3*c**6*f) - 4*csc(e + f*x)**11/(11*a**3*c**6*f) + 17*csc(e + f*x)**9/(9*a**3*c**6*f) - 4*csc(e + f*x)**7/(a**3*c**6*f) + 22*csc(e + f*x)**5/(5*a**3*c**6*f) - 8*csc(e + f*x)**3/(3*a**3*c**6*f) + csc(e + f*x)/(a**3*c**6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_64():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)
    F = -256*c**4*(a*sec(e + f*x) + a)*tan(e + f*x)/(315*f*sqrt(-c*sec(e + f*x) + c)) - 64*c**3*(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(105*f) - 8*c**2*(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(21*f) - 2*c*(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(9*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_65():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)
    F = -64*c**3*(a*sec(e + f*x) + a)*tan(e + f*x)/(105*f*sqrt(-c*sec(e + f*x) + c)) - 16*c**2*(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(35*f) - 2*c*(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_66():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)
    F = -8*c**2*(a*sec(e + f*x) + a)*tan(e + f*x)/(15*f*sqrt(-c*sec(e + f*x) + c)) - 2*c*(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_67():
    f = (a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -2*c*(a*sec(e + f*x) + a)*tan(e + f*x)/(3*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_68():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/sqrt(-c*sec(e + f*x) + c)
    F = 2*a*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c)) - 2*sqrt(2)*a*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_69():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -a*tan(e + f*x)/(f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + sqrt(2)*a*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_70():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = -a*tan(e + f*x)/(2*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) + a*tan(e + f*x)/(8*c*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + sqrt(2)*a*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(16*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_71():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)
    F = -256*c**4*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(1155*f*sqrt(-c*sec(e + f*x) + c)) - 64*c**3*(a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(231*f) - 8*c**2*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(33*f) - 2*c*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(11*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_72():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)
    F = -64*c**3*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(315*f*sqrt(-c*sec(e + f*x) + c)) - 16*c**2*(a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(63*f) - 2*c*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(9*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_73():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)
    F = -8*c**2*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(35*f*sqrt(-c*sec(e + f*x) + c)) - 2*c*(a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_74():
    f = (a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -2*c*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(5*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_75():
    f = (a*sec(e + f*x) + a)**2*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(7)/2)
    F = a**2*tan(e + f*x)/(4*c*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - a**2*tan(e + f*x)/(16*c**2*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - sqrt(2)*a**2*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(32*c**(sympy.S(7)/2)*f) - (a**2*sec(e + f*x) + a**2)*tan(e + f*x)/(3*f*(-c*sec(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_76():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)
    F = -256*c**4*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(3003*f*sqrt(-c*sec(e + f*x) + c)) - 64*c**3*(a*sec(e + f*x) + a)**3*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(429*f) - 24*c**2*(a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(143*f) - 2*c*(a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(13*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_77():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)
    F = -64*c**3*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(693*f*sqrt(-c*sec(e + f*x) + c)) - 16*c**2*(a*sec(e + f*x) + a)**3*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(99*f) - 2*c*(a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(11*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_78():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)
    F = -8*c**2*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(63*f*sqrt(-c*sec(e + f*x) + c)) - 2*c*(a*sec(e + f*x) + a)**3*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(9*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_79():
    f = (a*sec(e + f*x) + a)**3*sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -2*c*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(7*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_80():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/sqrt(-c*sec(e + f*x) + c)
    F = 8*a**3*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c)) - 8*sqrt(2)*a**3*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(sqrt(c)*f) + 2*a*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(5*f*sqrt(-c*sec(e + f*x) + c)) + (4*a**3*sec(e + f*x) + 4*a**3)*tan(e + f*x)/(3*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_81():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -10*a**3*tan(e + f*x)/(c*f*sqrt(-c*sec(e + f*x) + c)) + 10*sqrt(2)*a**3*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(c**(sympy.S(3)/2)*f) - a*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - (5*a**3*sec(e + f*x) + 5*a**3)*tan(e + f*x)/(3*c*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_82():
    f = (a*sec(e + f*x) + a)**3*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = 15*a**3*tan(e + f*x)/(4*c**2*f*sqrt(-c*sec(e + f*x) + c)) - 15*sqrt(2)*a**3*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(4*c**(sympy.S(5)/2)*f) - a*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(2*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) + (5*a**3*sec(e + f*x) + 5*a**3)*tan(e + f*x)/(4*c*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_83():
    f = (-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = 2*c*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) + 128*c**4*tan(e + f*x)/(5*a*f*sqrt(-c*sec(e + f*x) + c)) + 32*c**3*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(5*a*f) + 12*c**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(5*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_84():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = 2*c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) + 32*c**3*tan(e + f*x)/(3*a*f*sqrt(-c*sec(e + f*x) + c)) + 8*c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(3*a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_85():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = 2*c*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) + 4*c**2*tan(e + f*x)/(a*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_86():
    f = sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = 2*c*tan(e + f*x)/(f*(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_87():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    F = tan(e + f*x)/(f*(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(2*a*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_88():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = tan(e + f*x)/(f*(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 3*tan(e + f*x)/(4*a*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 3*sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(8*a*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_89():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = tan(e + f*x)/(f*(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - 5*tan(e + f*x)/(8*a*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - 15*tan(e + f*x)/(32*a*c*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 15*sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(64*a*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_90():
    f = (-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = -4*c**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(f*(a**2*sec(e + f*x) + a**2)) + 2*c*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - 64*c**4*tan(e + f*x)/(3*a**2*f*sqrt(-c*sec(e + f*x) + c)) - 16*c**3*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_91():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = -8*c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(3*f*(a**2*sec(e + f*x) + a**2)) + 2*c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - 16*c**3*tan(e + f*x)/(3*a**2*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_92():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = -4*c**2*tan(e + f*x)/(3*f*(a**2*sec(e + f*x) + a**2)*sqrt(-c*sec(e + f*x) + c)) + 2*c*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_93():
    f = sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = 2*c*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_94():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c))
    F = tan(e + f*x)/(2*f*(a**2*sec(e + f*x) + a**2)*sqrt(-c*sec(e + f*x) + c)) + tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c)) - sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(4*a**2*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_95():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = 5*tan(e + f*x)/(6*f*(a**2*sec(e + f*x) + a**2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 5*tan(e + f*x)/(8*a**2*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 5*sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(16*a**2*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_96():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = 7*tan(e + f*x)/(6*f*(a**2*sec(e + f*x) + a**2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) + tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - 35*tan(e + f*x)/(48*a**2*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - 35*tan(e + f*x)/(64*a**2*c*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 35*sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(128*a**2*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_97():
    f = (-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 16*c**3*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(5*f*(a**3*sec(e + f*x) + a**3)) + 2*c*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) - 4*c**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(5*a*f*(a*sec(e + f*x) + a)**2) + 32*c**4*tan(e + f*x)/(5*a**3*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_98():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 16*c**3*tan(e + f*x)/(15*f*(a**3*sec(e + f*x) + a**3)*sqrt(-c*sec(e + f*x) + c)) + 2*c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) - 8*c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_99():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 2*c*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) - 4*c**2*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_100():
    f = sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = 2*c*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_101():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*sqrt(-c*sec(e + f*x) + c))
    F = tan(e + f*x)/(4*f*(a**3*sec(e + f*x) + a**3)*sqrt(-c*sec(e + f*x) + c)) + tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3*sqrt(-c*sec(e + f*x) + c)) + tan(e + f*x)/(6*a*f*(a*sec(e + f*x) + a)**2*sqrt(-c*sec(e + f*x) + c)) - sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(8*a**3*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_102():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = 7*tan(e + f*x)/(12*f*(a**3*sec(e + f*x) + a**3)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + 7*tan(e + f*x)/(30*a*f*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 7*tan(e + f*x)/(16*a**3*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 7*sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(32*a**3*c**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_103():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = 21*tan(e + f*x)/(20*f*(a**3*sec(e + f*x) + a**3)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) + tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) + 3*tan(e + f*x)/(10*a*f*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - 21*tan(e + f*x)/(32*a**3*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - 63*tan(e + f*x)/(128*a**3*c*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 63*sqrt(2)*atan(sqrt(2)*sqrt(c)*tan(e + f*x)/(2*sqrt(-c*sec(e + f*x) + c)))/(256*a**3*c**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_104():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)
    F = a*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_105():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)
    F = a*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_106():
    f = sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -c*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_107():
    f = sqrt(a*sec(e + f*x) + a)*sec(e + f*x)/sqrt(-c*sec(e + f*x) + c)
    F = a*log(1 - sec(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_108():
    f = sqrt(a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(2*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_109():
    f = sqrt(a*sec(e + f*x) + a)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = -a*tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_110():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)
    F = a**2*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*tan(e + f*x)/(10*f*sqrt(a*sec(e + f*x) + a)) + a*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*tan(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_111():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)
    F = a**2*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(6*f*sqrt(a*sec(e + f*x) + a)) + a*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_112():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)
    F = -c**2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(3*f*sqrt(-c*sec(e + f*x) + c)) - c*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_113():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -c*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_114():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/sqrt(-c*sec(e + f*x) + c)
    F = 2*a**2*log(1 - sec(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + a*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_115():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -a**2*log(1 - sec(e + f*x))*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_116():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = -(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(4*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_117():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(7)/2)
    F = -(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(6*f*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)) - (a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(24*c*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_118():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(9)/2)
    F = a**2*tan(e + f*x)/(12*c*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)) - a*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(4*f*(-c*sec(e + f*x) + c)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_119():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(11)/2)
    F = a**2*tan(e + f*x)/(20*c*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(9)/2)) - a*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(5*f*(-c*sec(e + f*x) + c)**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_120():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*sec(e + f*x)
    F = a**3*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*tan(e + f*x)/(15*f*sqrt(a*sec(e + f*x) + a)) + 2*a**2*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*tan(e + f*x)/(15*f) + a*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)*tan(e + f*x)/(6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_121():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)
    F = -2*c**3*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(15*f*sqrt(-c*sec(e + f*x) + c)) - c**2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(5*f) - c*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_122():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)
    F = -c**2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(6*f*sqrt(-c*sec(e + f*x) + c)) - c*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_123():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -c*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(3*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_124():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/sqrt(-c*sec(e + f*x) + c)
    F = 4*a**3*log(1 - sec(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 2*a**2*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c)) + a*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_125():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -4*a**3*log(1 - sec(e + f*x))*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - 2*a**2*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(c*f*sqrt(-c*sec(e + f*x) + c)) - a*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_126():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = a**3*log(1 - sec(e + f*x))*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + a**2*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(c*f*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - a*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_127():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(7)/2)
    F = -(a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(6*f*(-c*sec(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_128():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(9)/2)
    F = -(a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(8*f*(-c*sec(e + f*x) + c)**(sympy.S(9)/2)) - (a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(48*c*f*(-c*sec(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_129():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(11)/2)
    F = -(a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(10*f*(-c*sec(e + f*x) + c)**(sympy.S(11)/2)) - (a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(40*c*f*(-c*sec(e + f*x) + c)**(sympy.S(9)/2)) - (a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(240*c**2*f*(-c*sec(e + f*x) + c)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_130():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/sqrt(a*sec(e + f*x) + a)
    F = -4*c**3*log(sec(e + f*x) + 1)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - 2*c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) - c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_131():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/sqrt(a*sec(e + f*x) + a)
    F = -2*c**2*log(sec(e + f*x) + 1)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - c*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_132():
    f = sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)/sqrt(a*sec(e + f*x) + a)
    F = -c*log(sec(e + f*x) + 1)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_133():
    f = sec(e + f*x)/(sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    F = -tan(e + f*x)*atanh(cos(e + f*x))/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_134():
    f = sec(e + f*x)/(sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = -tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - tan(e + f*x)*atanh(cos(e + f*x))/(2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_135():
    f = sec(e + f*x)/(sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = -tan(e + f*x)/(4*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - tan(e + f*x)/(4*c*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - tan(e + f*x)*atanh(cos(e + f*x))/(4*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_136():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + 4*c**3*log(sec(e + f*x) + 1)*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 2*c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_137():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = c*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + c**2*log(sec(e + f*x) + 1)*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_138():
    f = sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(2*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_139():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c))
    F = tan(e + f*x)/(2*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)) - tan(e + f*x)*atanh(cos(e + f*x))/(2*a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_140():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = -tan(e + f*x)*atanh(cos(e + f*x))/(2*a*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + csc(e + f*x)/(2*a*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_141():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = -tan(e + f*x)/(4*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - 3*tan(e + f*x)*atanh(cos(e + f*x))/(8*a*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 3*csc(e + f*x)/(8*a*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_142():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) - c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - c**3*log(sec(e + f*x) + 1)*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_143():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(4*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_144():
    f = sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = c*tan(e + f*x)/(2*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_145():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c))
    F = tan(e + f*x)/(4*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)) + tan(e + f*x)/(4*a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)) - tan(e + f*x)*atanh(cos(e + f*x))/(4*a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_146():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = tan(e + f*x)/(4*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) - 3*tan(e + f*x)*atanh(cos(e + f*x))/(8*a**2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 3*csc(e + f*x)/(8*a**2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_147():
    f = sec(e + f*x)/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = -3*tan(e + f*x)*atanh(cos(e + f*x))/(8*a**2*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - cot(e + f*x)**2*csc(e + f*x)/(4*a**2*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 3*csc(e + f*x)/(8*a**2*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_148():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**n*sec(e + f*x)
    F = -2**(n + sympy.S.Half)*c*(1 - sec(e + f*x))**(sympy.S.Half - n)*(a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(n - 1)*tan(e + f*x)*hyper((m + sympy.S.Half, sympy.S.Half - n), (m + sympy.S(3)/2,), sec(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_149():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**2*sec(e + f*x)
    F = 2**(m + sympy.S.Half)*a*(a*sec(e + f*x) + a)**(m - 1)*(-c*sec(e + f*x) + c)**2*(sec(e + f*x) + 1)**(sympy.S.Half - m)*tan(e + f*x)*hyper((sympy.S(5)/2, sympy.S.Half - m), (sympy.S(7)/2,), sympy.S.Half - sec(e + f*x)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_150():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = 2**(m + sympy.S.Half)*a*(a*sec(e + f*x) + a)**(m - 1)*(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)**(sympy.S.Half - m)*tan(e + f*x)*hyper((sympy.S(3)/2, sympy.S.Half - m), (sympy.S(5)/2,), sympy.S.Half - sec(e + f*x)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_151():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)/(-c*sec(e + f*x) + c)
    F = -2**(m + sympy.S.Half)*a*(a*sec(e + f*x) + a)**(m - 1)*(sec(e + f*x) + 1)**(sympy.S.Half - m)*tan(e + f*x)*hyper((sympy.S(-1)/2, sympy.S.Half - m), (sympy.S.Half,), sympy.S.Half - sec(e + f*x)/2)/(f*(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_152():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)/(-c*sec(e + f*x) + c)**2
    F = -2**(m + sympy.S.Half)*a*(a*sec(e + f*x) + a)**(m - 1)*(sec(e + f*x) + 1)**(sympy.S.Half - m)*tan(e + f*x)*hyper((sympy.S(-3)/2, sympy.S.Half - m), (sympy.S(-1)/2,), sympy.S.Half - sec(e + f*x)/2)/(3*f*(-c*sec(e + f*x) + c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_153():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*sec(e + f*x)
    F = -64*c**3*(a*sec(e + f*x) + a)**m*tan(e + f*x)/(f*(2*m + 5)*sqrt(-c*sec(e + f*x) + c)*(4*m**2 + 8*m + 3)) - 16*c**2*(a*sec(e + f*x) + a)**m*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*(4*m**2 + 16*m + 15)) - 2*c*(a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(f*(2*m + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_154():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*sec(e + f*x)
    F = -8*c**2*(a*sec(e + f*x) + a)**m*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c)*(4*m**2 + 8*m + 3)) - 2*c*(a*sec(e + f*x) + a)**m*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*(2*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_155():
    f = (a*sec(e + f*x) + a)**m*sqrt(-c*sec(e + f*x) + c)*sec(e + f*x)
    F = -2*c*(a*sec(e + f*x) + a)**m*tan(e + f*x)/(f*(2*m + 1)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_156():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)/sqrt(-c*sec(e + f*x) + c)
    F = -(a*sec(e + f*x) + a)**m*tan(e + f*x)*hyper((1, m + sympy.S.Half), (m + sympy.S(3)/2,), sec(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_157():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -(a*sec(e + f*x) + a)**m*tan(e + f*x)*hyper((2, m + sympy.S.Half), (m + sympy.S(3)/2,), sec(e + f*x)/2 + sympy.S.Half)/(2*c*f*(2*m + 1)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_158():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = -(a*sec(e + f*x) + a)**m*tan(e + f*x)*hyper((3, m + sympy.S.Half), (m + sympy.S(3)/2,), sec(e + f*x)/2 + sympy.S.Half)/(4*c**2*f*(2*m + 1)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_159():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(-m - 3)*sec(e + f*x)
    F = -(a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(-m - 3)*tan(e + f*x)/(f*(2*m + 1)) + 2*(a*sec(e + f*x) + a)**(m + 1)*(-c*sec(e + f*x) + c)**(-m - 3)*tan(e + f*x)/(a*f*(4*m**2 + 8*m + 3)) - 2*(a*sec(e + f*x) + a)**(m + 2)*(-c*sec(e + f*x) + c)**(-m - 3)*tan(e + f*x)/(a**2*f*(2*m + 1)*(4*m**2 + 16*m + 15))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_160():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(-m - 2)*sec(e + f*x)
    F = -(a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(-m - 2)*tan(e + f*x)/(f*(2*m + 1)) + (a*sec(e + f*x) + a)**(m + 1)*(-c*sec(e + f*x) + c)**(-m - 2)*tan(e + f*x)/(a*f*(4*m**2 + 8*m + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_161():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(-m - 1)*sec(e + f*x)
    F = -(a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(-m - 1)*tan(e + f*x)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_162():
    f = (a*sec(e + f*x) + a)**m*sec(e + f*x)/(-c*sec(e + f*x) + c)**m
    F = -2**(sympy.S.Half - m)*c*(1 - sec(e + f*x))**(m + sympy.S.Half)*(a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(-m - 1)*tan(e + f*x)*hyper((m + sympy.S.Half, m + sympy.S.Half), (m + sympy.S(3)/2,), sec(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_163():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(1 - m)*sec(e + f*x)
    F = -2**(sympy.S(3)/2 - m)*c*(1 - sec(e + f*x))**(m + sympy.S(-1)/2)*(a*sec(e + f*x) + a)**m*tan(e + f*x)*hyper((m + sympy.S(-1)/2, m + sympy.S.Half), (m + sympy.S(3)/2,), sec(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*(-c*sec(e + f*x) + c)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_164():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(2 - m)*sec(e + f*x)
    F = -2**(sympy.S(5)/2 - m)*c**2*(1 - sec(e + f*x))**(m + sympy.S(-1)/2)*(a*sec(e + f*x) + a)**m*tan(e + f*x)*hyper((m + sympy.S(-3)/2, m + sympy.S.Half), (m + sympy.S(3)/2,), sec(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1)*(-c*sec(e + f*x) + c)**m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_165():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)*sec(e + f*x)**2
    F = -a**3*c*tan(e + f*x)**5/(5*f) - 2*a**3*c*tan(e + f*x)**3/(3*f) - a**3*c*tan(e + f*x)*sec(e + f*x)**3/(2*f) + a**3*c*tan(e + f*x)*sec(e + f*x)/(4*f) + a**3*c*atanh(sin(e + f*x))/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_166():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)*sec(e + f*x)**2
    F = -a**2*c*tan(e + f*x)**3/(3*f) - a**2*c*tan(e + f*x)*sec(e + f*x)**3/(4*f) + a**2*c*tan(e + f*x)*sec(e + f*x)/(8*f) + a**2*c*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_167():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)*sec(e + f*x)**2
    F = -a*c*tan(e + f*x)**3/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_168():
    f = (-c*sec(e + f*x) + c)*sec(e + f*x)**2/(a*sec(e + f*x) + a)
    F = -2*c*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) - c*tan(e + f*x)/(a*f) + 2*c*atanh(sin(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_169():
    f = (-c*sec(e + f*x) + c)*sec(e + f*x)**2/(a*sec(e + f*x) + a)**2
    F = -2*c*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - c*atanh(sin(e + f*x))/(a**2*f) + 7*c*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_170():
    f = (-c*sec(e + f*x) + c)*sec(e + f*x)**2/(a*sec(e + f*x) + a)**3
    F = -4*c*tan(e + f*x)/(15*f*(a**3*sec(e + f*x) + a**3)) - 2*c*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) + 11*c*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_171():
    f = (g*sec(e + f*x))**p*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)
    F = -a**2*c*(g*sec(e + f*x))**p*(cos(e + f*x)**2)**(p/2 + sympy.S(3)/2)*tan(e + f*x)**3*hyper((sympy.S(3)/2, p/2 + sympy.S(3)/2), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f) - a**2*c*(g*sec(e + f*x))**(p + 1)*(cos(e + f*x)**2)**(p/2 + 2)*tan(e + f*x)**3*hyper((sympy.S(3)/2, p/2 + 2), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f*g)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_172():
    f = (g*sec(e + f*x))**p*(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)
    F = -a*c*(g*sec(e + f*x))**p*(cos(e + f*x)**2)**(p/2 + sympy.S(3)/2)*tan(e + f*x)**3*hyper((sympy.S(3)/2, p/2 + sympy.S(3)/2), (sympy.S(5)/2,), sin(e + f*x)**2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_173():
    f = (g*sec(e + f*x))**p*(-c*sec(e + f*x) + c)/(a*sec(e + f*x) + a)
    F = -2*c*(g*sec(e + f*x))**p*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) - c*g*(g*sec(e + f*x))**(p - 1)*(1 - 2*p)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - p/2), (sympy.S(3)/2 - p/2,), cos(e + f*x)**2)/(a*f*(1 - p)*sqrt(sin(e + f*x)**2)) + 2*c*(g*sec(e + f*x))**p*sin(e + f*x)*hyper((sympy.S.Half, -p/2), (1 - p/2,), cos(e + f*x)**2)/(a*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_174():
    f = (g*sec(e + f*x))**p*(-c*sec(e + f*x) + c)/(a*sec(e + f*x) + a)**2
    F = -2*c*(g*sec(e + f*x))**p*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - c*g*(g*sec(e + f*x))**(p - 1)*(3 - 4*p)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - p/2), (sympy.S(3)/2 - p/2,), cos(e + f*x)**2)/(3*a**2*f*sqrt(sin(e + f*x)**2)) + c*(g*sec(e + f*x))**p*(5 - 4*p)*sin(e + f*x)*hyper((sympy.S.Half, -p/2), (1 - p/2,), cos(e + f*x)**2)/(3*a**2*f*sqrt(sin(e + f*x)**2)) - c*(g*sec(e + f*x))**p*(5 - 4*p)*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_175():
    f = sec(e + f*x)**2/(sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    F = log(tan(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_176():
    f = sqrt(a*sec(e + f*x) + a)*sec(e + f*x)/(c - d*sec(e + f*x))
    F = 2*sqrt(a)*atanh(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c - d)*sqrt(a*sec(e + f*x) + a)))/(sqrt(d)*f*sqrt(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_177():
    f = (c + d*sec(e + f*x))**4*(a*sec(e + f*x) + a)*sec(e + f*x)
    F = a*d*(24*c**3 + 130*c**2*d + 116*c*d**2 + 45*d**3)*tan(e + f*x)*sec(e + f*x)/(120*f) + a*(c + d*sec(e + f*x))**4*tan(e + f*x)/(5*f) + a*(c + d*sec(e + f*x))**3*(4*c + 5*d)*tan(e + f*x)/(20*f) + a*(c + d*sec(e + f*x))**2*(12*c**2 + 35*c*d + 16*d**2)*tan(e + f*x)/(60*f) + a*(8*c**4 + 16*c**3*d + 24*c**2*d**2 + 12*c*d**3 + 3*d**4)*atanh(sin(e + f*x))/(8*f) + a*(12*c**4 + 95*c**3*d + 112*c**2*d**2 + 80*c*d**3 + 16*d**4)*tan(e + f*x)/(30*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_178():
    f = (c + d*sec(e + f*x))**3*(a*sec(e + f*x) + a)*sec(e + f*x)
    F = a*d*(6*c**2 + 20*c*d + 9*d**2)*tan(e + f*x)*sec(e + f*x)/(24*f) + a*(c + d*sec(e + f*x))**3*tan(e + f*x)/(4*f) + a*(c + d*sec(e + f*x))**2*(3*c + 4*d)*tan(e + f*x)/(12*f) + a*(3*c**3 + 16*c**2*d + 12*c*d**2 + 4*d**3)*tan(e + f*x)/(6*f) + a*(8*c**3 + 12*c**2*d + 12*c*d**2 + 3*d**3)*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_179():
    f = (c + d*sec(e + f*x))**2*(a*sec(e + f*x) + a)*sec(e + f*x)
    F = a*d*(2*c + 3*d)*tan(e + f*x)*sec(e + f*x)/(6*f) + a*(c + d*sec(e + f*x))**2*tan(e + f*x)/(3*f) + 2*a*(c**2 + 3*c*d + d**2)*tan(e + f*x)/(3*f) + a*(2*c**2 + 2*c*d + d**2)*atanh(sin(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_180():
    f = (c + d*sec(e + f*x))*(a*sec(e + f*x) + a)*sec(e + f*x)
    F = a*d*tan(e + f*x)*sec(e + f*x)/(2*f) + a*(c + d)*tan(e + f*x)/f + a*(2*c + d)*atanh(sin(e + f*x))/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_181():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(c + d*sec(e + f*x))
    F = -2*a*sqrt(c - d)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(d*f*sqrt(c + d)) + a*atanh(sin(e + f*x))/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_182():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(c + d*sec(e + f*x))**2
    F = a*tan(e + f*x)/(f*(c + d)*(c + d*sec(e + f*x))) + 2*a*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*sqrt(c - d)*(c + d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_183():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(c + d*sec(e + f*x))**3
    F = a*(c - 2*d)*tan(e + f*x)/(f*(c + d)**2*(c + d*sec(e + f*x))*(2*c - 2*d)) + a*tan(e + f*x)/(f*(c + d*sec(e + f*x))**2*(2*c + 2*d)) + a*(2*c - d)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*(c - d)**(sympy.S(3)/2)*(c + d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_184():
    f = (a*sec(e + f*x) + a)*sec(e + f*x)/(c + d*sec(e + f*x))**4
    F = a*(c - 4*d)*(2*c - d)*tan(e + f*x)/(6*f*(c - d)**2*(c + d)**3*(c + d*sec(e + f*x))) + a*tan(e + f*x)/(f*(c + d*sec(e + f*x))**3*(3*c + 3*d)) + a*(2*c - 3*d)*tan(e + f*x)/(f*(c + d)**2*(c + d*sec(e + f*x))**2*(6*c - 6*d)) + a*(2*c**2 - 2*c*d + d**2)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*(c - d)**(sympy.S(5)/2)*(c + d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_185():
    f = (c + d*sec(e + f*x))*(a*sec(e + f*x) + a)**2*sec(e + f*x)
    F = a**2*(3*c + 2*d)*tan(e + f*x)*sec(e + f*x)/(6*f) + 2*a**2*(3*c + 2*d)*tan(e + f*x)/(3*f) + a**2*(3*c + 2*d)*atanh(sin(e + f*x))/(2*f) + d*(a*sec(e + f*x) + a)**2*tan(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_186():
    f = (c + d*sec(e + f*x))*(a*sec(e + f*x) + a)**3*sec(e + f*x)
    F = a**3*(4*c + 3*d)*tan(e + f*x)**3/(12*f) + 3*a**3*(4*c + 3*d)*tan(e + f*x)*sec(e + f*x)/(8*f) + a**3*(4*c + 3*d)*tan(e + f*x)/f + 5*a**3*(4*c + 3*d)*atanh(sin(e + f*x))/(8*f) + d*(a*sec(e + f*x) + a)**3*tan(e + f*x)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_187():
    f = (c + d*sec(e + f*x))*sec(e + f*x)/(a*sec(e + f*x) + a)
    F = (c - d)*tan(e + f*x)/(f*(a*sec(e + f*x) + a)) + d*atanh(sin(e + f*x))/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_188():
    f = (c + d*sec(e + f*x))*sec(e + f*x)/(a*sec(e + f*x) + a)**2
    F = (c - d)*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) + (c + 2*d)*tan(e + f*x)/(3*f*(a**2*sec(e + f*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_189():
    f = (c + d*sec(e + f*x))**2*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = (c - d)**2*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) + (2*c**2 + 6*c*d + 7*d**2)*tan(e + f*x)/(15*f*(a**3*sec(e + f*x) + a**3)) + (c + 4*d)*(2*c - 2*d)*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_190():
    f = (c + d*sec(e + f*x))*sec(e + f*x)/(a*sec(e + f*x) + a)**3
    F = (c - d)*tan(e + f*x)/(5*f*(a*sec(e + f*x) + a)**3) + (2*c + 3*d)*tan(e + f*x)/(15*f*(a**3*sec(e + f*x) + a**3)) + (2*c + 3*d)*tan(e + f*x)/(15*a*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_191():
    f = sqrt(a*sec(e + f*x) + a)*sec(e + f*x)/sqrt(c + d*sec(e + f*x))
    F = 2*sqrt(a)*atanh(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_192():
    f = sqrt(c + d*sec(e + f*x))*sec(e + f*x)/sqrt(a*sec(e + f*x) + a)
    F = 2*sqrt(d)*atanh(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f) + sqrt(2)*sqrt(c - d)*atan(sqrt(2)*sqrt(a)*sqrt(c - d)*tan(e + f*x)/(2*sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_193():
    f = sec(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = sqrt(2)*atan(sqrt(2)*sqrt(a)*sqrt(c - d)*tan(e + f*x)/(2*sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*sqrt(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_194():
    f = sec(e + f*x)**2/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(a)*sqrt(c - d)*tan(e + f*x)/(2*sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*sqrt(c - d)) + 2*atanh(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*sqrt(d)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_195():
    f = sqrt(a*sec(e + f*x) + a)*sec(e + f*x)/(c + d*sec(e + f*x))
    F = 2*sqrt(a)*atan(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(sqrt(d)*f*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_196():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)*sqrt(a*sec(e + f*x) + a)/(c + d*sec(e + f*x))
    F = -2*sqrt(a)*sqrt(c)*g**(sympy.S(3)/2)*atanh(sqrt(a)*sqrt(c)*sqrt(g)*tan(e + f*x)/(sqrt(g*sec(e + f*x))*sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(d*f*sqrt(c + d)) + 2*sqrt(a)*g**(sympy.S(3)/2)*atanh(sqrt(a)*sqrt(g)*tan(e + f*x)/(sqrt(g*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_197():
    f = sec(e + f*x)/((c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = -2*sqrt(d)*atan(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*(c - d)*sqrt(c + d)) + sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_198():
    f = sec(e + f*x)**2/((c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = 2*c*atan(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*sqrt(d)*f*(c - d)*sqrt(c + d)) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_199():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)/((c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = 2*sqrt(c)*g**(sympy.S(3)/2)*atanh(sqrt(a)*sqrt(c)*sqrt(g)*tan(e + f*x)/(sqrt(g*sec(e + f*x))*sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*(c - d)*sqrt(c + d)) - sqrt(2)*g**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(g)*tan(e + f*x)/(2*sqrt(g*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*(c - d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_200():
    f = (g*sec(e + f*x))**(sympy.S(5)/2)/((c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = -2*c**(sympy.S(3)/2)*g**(sympy.S(5)/2)*atanh(sqrt(a)*sqrt(c)*sqrt(g)*tan(e + f*x)/(sqrt(g*sec(e + f*x))*sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*d*f*(c - d)*sqrt(c + d)) + sqrt(2)*g**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sqrt(g)*tan(e + f*x)/(2*sqrt(g*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*(c - d)) + 2*g**(sympy.S(5)/2)*atanh(sqrt(a)*sqrt(g)*tan(e + f*x)/(sqrt(g*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*d*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_201():
    f = (a + b*sec(e + f*x))*(c + d*sec(e + f*x))**4*sec(e + f*x)
    F = b*(c + d*sec(e + f*x))**4*tan(e + f*x)/(5*f) + d*(130*a*c**2*d + 45*a*d**3 + 24*b*c**3 + 116*b*c*d**2)*tan(e + f*x)*sec(e + f*x)/(120*f) + (c + d*sec(e + f*x))**3*(5*a*d + 4*b*c)*tan(e + f*x)/(20*f) + (c + d*sec(e + f*x))**2*(35*a*c*d + 12*b*c**2 + 16*b*d**2)*tan(e + f*x)/(60*f) + (8*a*c**4 + 24*a*c**2*d**2 + 3*a*d**4 + 16*b*c**3*d + 12*b*c*d**3)*atanh(sin(e + f*x))/(8*f) + (95*a*c**3*d + 80*a*c*d**3 + 12*b*c**4 + 112*b*c**2*d**2 + 16*b*d**4)*tan(e + f*x)/(30*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_202():
    f = (a + b*sec(e + f*x))*(c + d*sec(e + f*x))**3*sec(e + f*x)
    F = b*(c + d*sec(e + f*x))**3*tan(e + f*x)/(4*f) + d*(20*a*c*d + 6*b*c**2 + 9*b*d**2)*tan(e + f*x)*sec(e + f*x)/(24*f) + (c + d*sec(e + f*x))**2*(4*a*d + 3*b*c)*tan(e + f*x)/(12*f) + (4*a*d*(4*c**2 + d**2) + 3*b*(c**3 + 4*c*d**2))*tan(e + f*x)/(6*f) + (8*a*c**3 + 12*a*c*d**2 + 12*b*c**2*d + 3*b*d**3)*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_203():
    f = (a + b*sec(e + f*x))*(c + d*sec(e + f*x))**2*sec(e + f*x)
    F = b*(c + d*sec(e + f*x))**2*tan(e + f*x)/(3*f) + d*(3*a*d + 2*b*c)*tan(e + f*x)*sec(e + f*x)/(6*f) + (a*(2*c**2 + d**2) + 2*b*c*d)*atanh(sin(e + f*x))/(2*f) + (6*a*c*d + 2*b*(c**2 + d**2))*tan(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_204():
    f = (a + b*sec(e + f*x))*(c + d*sec(e + f*x))*sec(e + f*x)
    F = b*d*tan(e + f*x)*sec(e + f*x)/(2*f) + (2*a*c + b*d)*atanh(sin(e + f*x))/(2*f) + (a*d + b*c)*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_205():
    f = (a + b*sec(e + f*x))*sec(e + f*x)/(c + d*sec(e + f*x))
    F = b*atanh(sin(e + f*x))/(d*f) - (-2*a*d + 2*b*c)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(d*f*sqrt(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_206():
    f = (a + b*sec(e + f*x))*sec(e + f*x)/(c + d*sec(e + f*x))**2
    F = (-a*d + b*c)*tan(e + f*x)/(f*(c + d*sec(e + f*x))*(c**2 - d**2)) + (2*a*c - 2*b*d)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*(c - d)**(sympy.S(3)/2)*(c + d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_207():
    f = (a + b*sec(e + f*x))*sec(e + f*x)/(c + d*sec(e + f*x))**3
    F = -(3*a*c*d - b*(c**2 + 2*d**2))*tan(e + f*x)/(2*f*(c + d*sec(e + f*x))*(c**2 - d**2)**2) + (-a*d + b*c)*tan(e + f*x)/(f*(c + d*sec(e + f*x))**2*(2*c**2 - 2*d**2)) - (-a*(2*c**2 + d**2) + 3*b*c*d)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*(c - d)**(sympy.S(5)/2)*(c + d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_208():
    f = (a + b*sec(e + f*x))*sec(e + f*x)/(c + d*sec(e + f*x))**4
    F = (-11*a*c**2*d - 4*a*d**3 + 2*b*c**3 + 13*b*c*d**2)*tan(e + f*x)/(6*f*(c + d*sec(e + f*x))*(c**2 - d**2)**3) + (-5*a*c*d + 2*b*c**2 + 3*b*d**2)*tan(e + f*x)/(6*f*(c + d*sec(e + f*x))**2*(c**2 - d**2)**2) + (-a*d + b*c)*tan(e + f*x)/(f*(c + d*sec(e + f*x))**3*(3*c**2 - 3*d**2)) + (2*a*c**3 + 3*a*c*d**2 - 4*b*c**2*d - b*d**3)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*(c - d)**(sympy.S(7)/2)*(c + d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_209():
    f = (c + d*sec(e + f*x))**4*sec(e + f*x)/(a + b*sec(e + f*x))
    F = d**4*tan(e + f*x)**3/(3*b*f) + d**4*tan(e + f*x)/(b*f) + d**3*(-a*d + 4*b*c)*tan(e + f*x)*sec(e + f*x)/(2*b**2*f) + d**3*(-a*d + 4*b*c)*atanh(sin(e + f*x))/(2*b**2*f) + d**2*(a**2*d**2 - 4*a*b*c*d + 6*b**2*c**2)*tan(e + f*x)/(b**3*f) + d*(-a*d + 2*b*c)*(a**2*d**2 - 2*a*b*c*d + 2*b**2*c**2)*atanh(sin(e + f*x))/(b**4*f) + 2*(-a*d + b*c)**4*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(b**4*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_210():
    f = (c + d*sec(e + f*x))**3*sec(e + f*x)/(a + b*sec(e + f*x))
    F = d**3*tan(e + f*x)*sec(e + f*x)/(2*b*f) + d**3*atanh(sin(e + f*x))/(2*b*f) + d**2*(-a*d + 3*b*c)*tan(e + f*x)/(b**2*f) + d*(a**2*d**2 - 3*a*b*c*d + 3*b**2*c**2)*atanh(sin(e + f*x))/(b**3*f) + 2*(-a*d + b*c)**3*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(b**3*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_211():
    f = (c + d*sec(e + f*x))**2*sec(e + f*x)/(a + b*sec(e + f*x))
    F = d**2*tan(e + f*x)/(b*f) + d*(-a*d + 2*b*c)*atanh(sin(e + f*x))/(b**2*f) + 2*(-a*d + b*c)**2*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(b**2*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_212():
    f = (c + d*sec(e + f*x))*sec(e + f*x)/(a + b*sec(e + f*x))
    F = d*atanh(sin(e + f*x))/(b*f) + (-2*a*d + 2*b*c)*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(b*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_213():
    f = sec(e + f*x)/((a + b*sec(e + f*x))*(c + d*sec(e + f*x)))
    F = 2*b*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(f*sqrt(a - b)*sqrt(a + b)*(-a*d + b*c)) - 2*d*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*sqrt(c - d)*sqrt(c + d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_214():
    f = sec(e + f*x)/((a + b*sec(e + f*x))*(c + d*sec(e + f*x))**2)
    F = 2*b**2*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(f*sqrt(a - b)*sqrt(a + b)*(-a*d + b*c)**2) + d**2*sin(e + f*x)/(f*(c**2 - d**2)*(-a*d + b*c)*(c*cos(e + f*x) + d)) - 2*d*(-a*c*d + 2*b*c**2 - b*d**2)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*(c - d)**(sympy.S(3)/2)*(c + d)**(sympy.S(3)/2)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_215():
    f = (c + d*sec(e + f*x))**5*sec(e + f*x)/(a + b*sec(e + f*x))**2
    F = d**5*tan(e + f*x)**3/(3*b**2*f) + d**5*tan(e + f*x)/(b**2*f) + d**4*(-2*a*d + 5*b*c)*tan(e + f*x)*sec(e + f*x)/(2*b**3*f) + d**4*(-2*a*d + 5*b*c)*atanh(sin(e + f*x))/(2*b**3*f) + d**3*(3*a**2*d**2 - 10*a*b*c*d + 10*b**2*c**2)*tan(e + f*x)/(b**4*f) - (-a*d + b*c)**5*sin(e + f*x)/(b**4*f*(a**2 - b**2)*(a*cos(e + f*x) + b)) + d**2*(-4*a**3*d**3 + 15*a**2*b*c*d**2 - 20*a*b**2*c**2*d + 10*b**3*c**3)*atanh(sin(e + f*x))/(b**5*f) + 2*(-a*d + b*c)**5*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*b**3*f*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + 2*(-a*d + b*c)**4*(4*a*d + b*c)*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*b**5*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_216():
    f = (c + d*sec(e + f*x))**4*sec(e + f*x)/(a + b*sec(e + f*x))**2
    F = d**4*tan(e + f*x)*sec(e + f*x)/(2*b**2*f) + d**4*atanh(sin(e + f*x))/(2*b**2*f) + 2*d**3*(-a*d + 2*b*c)*tan(e + f*x)/(b**3*f) - (-a*d + b*c)**4*sin(e + f*x)/(b**3*f*(a**2 - b**2)*(a*cos(e + f*x) + b)) + d**2*(3*a**2*d**2 - 8*a*b*c*d + 6*b**2*c**2)*atanh(sin(e + f*x))/(b**4*f) + 2*(-a*d + b*c)**4*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*b**2*f*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + 2*(-a*d + b*c)**3*(3*a*d + b*c)*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*b**4*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_217():
    f = (c + d*sec(e + f*x))**3*sec(e + f*x)/(a + b*sec(e + f*x))**2
    F = d**3*tan(e + f*x)/(b**2*f) - (-a*d + b*c)**3*sin(e + f*x)/(b**2*f*(a**2 - b**2)*(a*cos(e + f*x) + b)) + d**2*(-2*a*d + 3*b*c)*atanh(sin(e + f*x))/(b**3*f) + 2*(-a*d + b*c)**3*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*b*f*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + 2*(-a*d + b*c)**2*(2*a*d + b*c)*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*b**3*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_218():
    f = (c + d*sec(e + f*x))**2*sec(e + f*x)/(a + b*sec(e + f*x))**2
    F = -(-a*d + b*c)**2*sin(e + f*x)/(b*f*(a**2 - b**2)*(a*cos(e + f*x) + b)) + d**2*atanh(sin(e + f*x))/(b**2*f) + 2*(-a*d + b*c)**2*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*f*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + (-2*a**2*d**2 + 2*b**2*c**2)*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(a*b**2*f*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_219():
    f = (c + d*sec(e + f*x))*sec(e + f*x)/(a + b*sec(e + f*x))**2
    F = -(-a*d + b*c)*tan(e + f*x)/(f*(a + b*sec(e + f*x))*(a**2 - b**2)) + (2*a*c - 2*b*d)*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(f*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_220():
    f = sec(e + f*x)/((a + b*sec(e + f*x))**2*(c + d*sec(e + f*x)))
    F = -b**2*sin(e + f*x)/(f*(a**2 - b**2)*(-a*d + b*c)*(a*cos(e + f*x) + b)) + 2*b*(-2*a**2*d + a*b*c + b**2*d)*atanh(sqrt(a - b)*tan(e/2 + f*x/2)/sqrt(a + b))/(f*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)*(-a*d + b*c)**2) + 2*d**2*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(f*sqrt(c - d)*sqrt(c + d)*(-a*d + b*c)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_221():
    f = sqrt(a + b*sec(e + f*x))*sec(e + f*x)/(c + d*sec(e + f*x))
    F = ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_222():
    f = sqrt(a + b*sec(e + f*x))*sec(e + f*x)/sqrt(c + d*sec(e + f*x))
    F = (Integer(2) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(((Symbol('a') + Symbol('b')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + Symbol('b')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_223():
    f = sec(e + f*x)/(sqrt(a + b*sec(e + f*x))*sqrt(c + d*sec(e + f*x)))
    F = 2*sqrt((1 - sec(e + f*x))*(-a*d + b*c)/((a + b)*(c + d*sec(e + f*x))))*sqrt(-(-a*d + b*c)*(sec(e + f*x) + 1)/((a - b)*(c + d*sec(e + f*x))))*sqrt(a + b)*(c + d*sec(e + f*x))*cot(e + f*x)*elliptic_f(asin(sqrt(a + b*sec(e + f*x))*sqrt(c + d)/(sqrt(a + b)*sqrt(c + d*sec(e + f*x)))), (a + b)*(c - d)/((a - b)*(c + d)))/(f*sqrt(c + d)*(-a*d + b*c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_224():
    f = sec(e + f*x)/(sqrt(3*sec(e + f*x) + 2)*sqrt(5*sec(e + f*x) - 4))
    F = 2*sqrt((1 - sec(e + f*x))/(4 - 5*sec(e + f*x)))*sqrt((sec(e + f*x) + 1)/(4 - 5*sec(e + f*x)))*(4 - 5*sec(e + f*x))*cot(e + f*x)*elliptic_f(asin(sqrt(5)*sqrt(3*sec(e + f*x) + 2)/(5*sqrt(5*sec(e + f*x) - 4))), 45)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_225():
    f = sec(e + f*x)/(sqrt(4 - 5*sec(e + f*x))*sqrt(3*sec(e + f*x) + 2))
    F = 2*sqrt(5)*I*sqrt((1 - sec(e + f*x))/(3*sec(e + f*x) + 2))*sqrt((sec(e + f*x) + 1)/(3*sec(e + f*x) + 2))*(3*sec(e + f*x) + 2)*cot(e + f*x)*elliptic_f(I*asinh(sqrt(5)*sqrt(4 - 5*sec(e + f*x))/sqrt(3*sec(e + f*x) + 2)), sympy.S(1)/45)/(15*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_226():
    f = sec(e + f*x)**2/(sqrt(a + b*sec(e + f*x))*sqrt(c + d*sec(e + f*x)))
    F = ((Integer(2) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(((Symbol('a') + Symbol('b')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(((Symbol('a') + Symbol('b')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_227():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)*sqrt(c + d*sec(e + f*x))/(a + b*sec(e + f*x))
    F = ((Integer(2) * Symbol('d') * Symbol('g') * sympy.sqrt(((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('g') * sympy.sqrt(((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_228():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)/((a + b*sec(e + f*x))*sqrt(c + d*sec(e + f*x)))
    F = (Integer(2) * Symbol('g') * sympy.sqrt(((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_229():
    f = sqrt(g*sec(e + f*x))*sqrt(c + d*sec(e + f*x))/(a + b*cos(e + f*x))
    F = ((Integer(2) * Symbol('d') * sympy.sqrt(((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('a') * Symbol('c')) + (Integer(-1) * (Symbol('b') * Symbol('d')))) * sympy.sqrt(((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_230():
    f = sqrt(a + b*sec(e + f*x))*sec(e + f*x)/(c*sec(e + f*x) + c)
    F = sqrt(a + b*sec(e + f*x))*sqrt(1/(sec(e + f*x) + 1))*elliptic_e(asin(tan(e + f*x)/(sec(e + f*x) + 1)), (a - b)/(a + b))/(c*f*sqrt((a + b*sec(e + f*x))/((a + b)*(sec(e + f*x) + 1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_231():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)*sqrt(a + b*sec(e + f*x))/(c*sec(e + f*x) + c)
    F = ((Symbol('g') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('g') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('g') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('g') * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('f') * (Symbol('c') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_232():
    f = sec(e + f*x)/(sqrt(a + b*sec(e + f*x))*(c*sec(e + f*x) + c))
    F = -2*sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*sqrt(a + b)*cot(e + f*x)*elliptic_f(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/(c*f*(a - b)) + sqrt(a + b*sec(e + f*x))*sqrt(1/(sec(e + f*x) + 1))*elliptic_e(asin(tan(e + f*x)/(sec(e + f*x) + 1)), (a - b)/(a + b))/(c*f*sqrt((a + b*sec(e + f*x))/((a + b)*(sec(e + f*x) + 1)))*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_233():
    f = sec(e + f*x)**2/(sqrt(a + b*sec(e + f*x))*(c*sec(e + f*x) + c))
    F = 2*a*sqrt(b*(1 - sec(e + f*x))/(a + b))*sqrt(-b*(sec(e + f*x) + 1)/(a - b))*sqrt(a + b)*cot(e + f*x)*elliptic_f(asin(sqrt(a + b*sec(e + f*x))/sqrt(a + b)), (a + b)/(a - b))/(b*c*f*(a - b)) - sqrt(a + b*sec(e + f*x))*sqrt(1/(sec(e + f*x) + 1))*elliptic_e(asin(tan(e + f*x)/(sec(e + f*x) + 1)), (a - b)/(a + b))/(c*f*sqrt((a + b*sec(e + f*x))/((a + b)*(sec(e + f*x) + 1)))*(a - b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_234():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)/(sqrt(a + b*sec(e + f*x))*(c*sec(e + f*x) + c))
    F = -g*sqrt(g*sec(e + f*x))*(a*cos(e + f*x) + b)*sin(e + f*x)/(f*(a - b)*sqrt(a + b*sec(e + f*x))*(c*cos(e + f*x) + c)) + g*sqrt(g*sec(e + f*x))*sqrt((a*cos(e + f*x) + b)/(a + b))*elliptic_f(e/2 + f*x/2, 2*a/(a + b))/(c*f*sqrt(a + b*sec(e + f*x))) + g*sqrt(g*sec(e + f*x))*(a*cos(e + f*x) + b)*elliptic_e(e/2 + f*x/2, 2*a/(a + b))/(c*f*sqrt((a*cos(e + f*x) + b)/(a + b))*(a - b)*sqrt(a + b*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_235():
    f = (g*sec(e + f*x))**(sympy.S(5)/2)/(sqrt(a + b*sec(e + f*x))*(c*sec(e + f*x) + c))
    F = (Integer(-1) * (((Symbol('g'))**(Integer(2)) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('c') * Symbol('f') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('g'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('g'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (((Symbol('g'))**(Integer(2)) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('f') * (Symbol('c') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_236():
    f = sqrt(a + b*sec(e + f*x))*sec(e + f*x)/(c + d*sec(e + f*x))
    F = ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_237():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)*sqrt(a + b*sec(e + f*x))/(c + d*sec(e + f*x))
    F = ((Integer(2) * Symbol('b') * Symbol('g') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('d') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('g') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('d') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_238():
    f = sec(e + f*x)/(sqrt(a + b*sec(e + f*x))*(c + d*sec(e + f*x)))
    F = (Integer(2) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_239():
    f = sec(e + f*x)**2/(sqrt(a + b*sec(e + f*x))*(c + d*sec(e + f*x)))
    F = ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('c') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('d') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_240():
    f = (g*sec(e + f*x))**(sympy.S(3)/2)/(sqrt(a + b*sec(e + f*x))*(c + d*sec(e + f*x)))
    F = (Integer(2) * Symbol('g') * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_241():
    f = (g*sec(e + f*x))**(sympy.S(5)/2)/(sqrt(a + b*sec(e + f*x))*(c + d*sec(e + f*x)))
    F = ((Integer(2) * (Symbol('g'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('d') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('g'))**(Integer(2)) * sympy.sqrt(((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('c')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('e') + (Symbol('f') * x))), ((Integer(2) * Symbol('a')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Symbol('g') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('d') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_242():
    f = tan(e + f*x)**4*sec(e + f*x)/(-c*sec(e + f*x) + c)**7
    F = cot(e/2 + f*x/2)**9/(36*c**7*f) - cot(e/2 + f*x/2)**7/(14*c**7*f) + cot(e/2 + f*x/2)**5/(20*c**7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_3_g_sec_pow_p_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_243():
    f = tan(e + f*x)**4*sec(e + f*x)/(-c*sec(e + f*x) + c)**8
    F = -cot(e/2 + f*x/2)**11/(88*c**8*f) + cot(e/2 + f*x/2)**9/(24*c**8*f) - 3*cot(e/2 + f*x/2)**7/(56*c**8*f) + cot(e/2 + f*x/2)**5/(40*c**8*f)
    assert integrate(f, x) == F

