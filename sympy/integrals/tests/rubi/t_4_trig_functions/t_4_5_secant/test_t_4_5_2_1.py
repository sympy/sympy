"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.5 Secant/4.5.2.1 (a+b sec)^m (c+d sec)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL, slow

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n, p = symbols('a b c d e f m n p')

def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_1():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**5
    F = a**2*c**5*x + 3*a**2*c**5*tan(e + f*x)**5/(5*f) - a**2*c**5*tan(e + f*x)**3*sec(e + f*x)**3/(6*f) - 3*a**2*c**5*tan(e + f*x)**3*sec(e + f*x)/(4*f) + a**2*c**5*tan(e + f*x)**3/(3*f) + a**2*c**5*tan(e + f*x)*sec(e + f*x)**3/(8*f) + 17*a**2*c**5*tan(e + f*x)*sec(e + f*x)/(16*f) - a**2*c**5*tan(e + f*x)/f - 19*a**2*c**5*atanh(sin(e + f*x))/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_2():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**4
    F = a**2*c**4*x + a**2*c**4*tan(e + f*x)**5/(5*f) - a**2*c**4*tan(e + f*x)**3*sec(e + f*x)/(2*f) + a**2*c**4*tan(e + f*x)**3/(3*f) + 3*a**2*c**4*tan(e + f*x)*sec(e + f*x)/(4*f) - a**2*c**4*tan(e + f*x)/f - 3*a**2*c**4*atanh(sin(e + f*x))/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_3():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**3
    F = a**2*c**3*x - 3*a**2*c**3*atanh(sin(e + f*x))/(8*f) + a**2*(-3*c**3*sec(e + f*x) + 4*c**3)*tan(e + f*x)**3/(12*f) - a**2*(-3*c**3*sec(e + f*x) + 8*c**3)*tan(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_4():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**2
    F = a**2*c**2*x + a**2*c**2*tan(e + f*x)**3/(3*f) - a**2*c**2*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_5():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)
    F = a**2*c*x + a**2*c*atanh(sin(e + f*x))/(2*f) - c*(a**2*sec(e + f*x) + 2*a**2)*tan(e + f*x)/(2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_6():
    f = (a*sec(e + f*x) + a)**2/(-c*sec(e + f*x) + c)
    F = a**2*x/c - a**2*atanh(sin(e + f*x))/(c*f) - 4*a**2*tan(e + f*x)/(c*f*(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_7():
    f = (a*sec(e + f*x) + a)**2/(-c*sec(e + f*x) + c)**2
    F = a**2*x/c**2 - 4*a**2*tan(e + f*x)/(3*c**2*f*(1 - sec(e + f*x))) - 4*a**2*tan(e + f*x)/(3*c**2*f*(1 - sec(e + f*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_8():
    f = (a*sec(e + f*x) + a)**2/(-c*sec(e + f*x) + c)**3
    F = a**2*x/c**3 - 23*a**2*tan(e + f*x)/(15*c**3*f*(1 - sec(e + f*x))) - 8*a**2*tan(e + f*x)/(15*c**3*f*(1 - sec(e + f*x))**2) - 4*a**2*tan(e + f*x)/(5*c**3*f*(1 - sec(e + f*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_9():
    f = (a*sec(e + f*x) + a)**2/(-c*sec(e + f*x) + c)**4
    F = a**2*x/c**4 - 164*a**2*tan(e + f*x)/(105*c**4*f*(1 - sec(e + f*x))) - 59*a**2*tan(e + f*x)/(105*c**4*f*(1 - sec(e + f*x))**2) - 12*a**2*tan(e + f*x)/(35*c**4*f*(1 - sec(e + f*x))**3) - 4*a**2*tan(e + f*x)/(7*c**4*f*(1 - sec(e + f*x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_10():
    f = (a*sec(e + f*x) + a)**2/(-c*sec(e + f*x) + c)**5
    F = a**2*x/c**5 - 494*a**2*tan(e + f*x)/(315*c**5*f*(1 - sec(e + f*x))) - 179*a**2*tan(e + f*x)/(315*c**5*f*(1 - sec(e + f*x))**2) - 37*a**2*tan(e + f*x)/(105*c**5*f*(1 - sec(e + f*x))**3) - 16*a**2*tan(e + f*x)/(63*c**5*f*(1 - sec(e + f*x))**4) - 4*a**2*tan(e + f*x)/(9*c**5*f*(1 - sec(e + f*x))**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_11():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**5
    F = a**3*c**5*x - a**3*c**5*tan(e + f*x)**7/(7*f) + a**3*c**5*tan(e + f*x)**5*sec(e + f*x)/(3*f) - a**3*c**5*tan(e + f*x)**5/(5*f) - 5*a**3*c**5*tan(e + f*x)**3*sec(e + f*x)/(12*f) + a**3*c**5*tan(e + f*x)**3/(3*f) + 5*a**3*c**5*tan(e + f*x)*sec(e + f*x)/(8*f) - a**3*c**5*tan(e + f*x)/f - 5*a**3*c**5*atanh(sin(e + f*x))/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_12():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**4
    F = a**3*c**4*x - 5*a**3*c**4*atanh(sin(e + f*x))/(16*f) - a**3*(-5*c**4*sec(e + f*x) + 6*c**4)*tan(e + f*x)**5/(30*f) + a**3*(-5*c**4*sec(e + f*x) + 8*c**4)*tan(e + f*x)**3/(24*f) - a**3*(-5*c**4*sec(e + f*x) + 16*c**4)*tan(e + f*x)/(16*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_13():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**3
    F = a**3*c**3*x - a**3*c**3*tan(e + f*x)**5/(5*f) + a**3*c**3*tan(e + f*x)**3/(3*f) - a**3*c**3*tan(e + f*x)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_14():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**2
    F = a**3*c**2*x + 3*a**3*c**2*atanh(sin(e + f*x))/(8*f) + c**2*(3*a**3*sec(e + f*x) + 4*a**3)*tan(e + f*x)**3/(12*f) - c**2*(3*a**3*sec(e + f*x) + 8*a**3)*tan(e + f*x)/(8*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_15():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)
    F = a**3*c*x - a**3*c*tan(e + f*x)**3/(3*f) - a**3*c*tan(e + f*x)*sec(e + f*x)/f - a**3*c*tan(e + f*x)/f + a**3*c*atanh(sin(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_16():
    f = (a*sec(e + f*x) + a)**3/(-c*sec(e + f*x) + c)
    F = a**3*x/c - a**3*tan(e + f*x)/(c*f) + 8*a**3*cot(e + f*x)/(c*f) - 4*a**3*atanh(sin(e + f*x))/(c*f) + 8*a**3*csc(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_17():
    f = (a*sec(e + f*x) + a)**3/(-c*sec(e + f*x) + c)**2
    F = a**3*x/c**2 + a**3*atanh(sin(e + f*x))/(c**2*f) + 4*a**3*tan(e + f*x)/(3*c**2*f*(1 - sec(e + f*x))) - 8*a**3*tan(e + f*x)/(3*c**2*f*(1 - sec(e + f*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_18():
    f = (a*sec(e + f*x) + a)**3/(-c*sec(e + f*x) + c)**3
    F = a**3*x/c**3 - 26*a**3*tan(e + f*x)/(15*c**3*f*(1 - sec(e + f*x))) + 4*a**3*tan(e + f*x)/(15*c**3*f*(1 - sec(e + f*x))**2) - 8*a**3*tan(e + f*x)/(5*c**3*f*(1 - sec(e + f*x))**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_19():
    f = (a*sec(e + f*x) + a)**3/(-c*sec(e + f*x) + c)**4
    F = a**3*x/c**4 - 167*a**3*tan(e + f*x)/(105*c**4*f*(1 - sec(e + f*x))) - 62*a**3*tan(e + f*x)/(105*c**4*f*(1 - sec(e + f*x))**2) + 4*a**3*tan(e + f*x)/(35*c**4*f*(1 - sec(e + f*x))**3) - 8*a**3*tan(e + f*x)/(7*c**4*f*(1 - sec(e + f*x))**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_20():
    f = (a*sec(e + f*x) + a)**3/(-c*sec(e + f*x) + c)**5
    F = a**3*x/c**5 - 496*a**3*tan(e + f*x)/(315*c**5*f*(1 - sec(e + f*x))) - 181*a**3*tan(e + f*x)/(315*c**5*f*(1 - sec(e + f*x))**2) - 38*a**3*tan(e + f*x)/(105*c**5*f*(1 - sec(e + f*x))**3) + 4*a**3*tan(e + f*x)/(63*c**5*f*(1 - sec(e + f*x))**4) - 8*a**3*tan(e + f*x)/(9*c**5*f*(1 - sec(e + f*x))**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_21():
    f = (-c*sec(e + f*x) + c)**4/(a*sec(e + f*x) + a)**2
    F = c**4*x/a**2 + c**4*tan(e + f*x)/(a**2*f) - 32*c**4*cot(e + f*x)**3/(3*a**2*f) - 16*c**4*cot(e + f*x)/(a**2*f) - 6*c**4*atanh(sin(e + f*x))/(a**2*f) + 32*c**4*csc(e + f*x)**3/(3*a**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_22():
    f = (-c*sec(e + f*x) + c)**3/(a*sec(e + f*x) + a)**2
    F = c**3*x/a**2 - c**3*atanh(sin(e + f*x))/(a**2*f) + 4*c**3*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1)) - 8*c**3*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_23():
    f = (-c*sec(e + f*x) + c)**2/(a*sec(e + f*x) + a)**2
    F = c**2*x/a**2 - 4*c**2*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1)) - 4*c**2*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_24():
    f = (-c*sec(e + f*x) + c)/(a*sec(e + f*x) + a)**2
    F = c*x/a**2 - 5*c*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1)) - 2*c*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_25():
    f = 1/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c))
    F = x/(a**2*c) - (1 - sec(e + f*x))*cot(e + f*x)**3/(3*a**2*c*f) + (3 - 2*sec(e + f*x))*cot(e + f*x)/(3*a**2*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_26():
    f = 1/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**2)
    F = x/(a**2*c**2) - cot(e + f*x)**3/(3*a**2*c**2*f) + cot(e + f*x)/(a**2*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_27():
    f = 1/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**3)
    F = x/(a**2*c**3) + (sec(e + f*x) + 1)*cot(e + f*x)**5/(5*a**2*c**3*f) - (4*sec(e + f*x) + 5)*cot(e + f*x)**3/(15*a**2*c**3*f) + (8*sec(e + f*x) + 15)*cot(e + f*x)/(15*a**2*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_28():
    f = 1/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**4)
    F = x/(a**2*c**4) - 2*cot(e + f*x)**7/(7*a**2*c**4*f) + cot(e + f*x)**5/(5*a**2*c**4*f) - cot(e + f*x)**3/(3*a**2*c**4*f) + cot(e + f*x)/(a**2*c**4*f) - 2*csc(e + f*x)**7/(7*a**2*c**4*f) + 6*csc(e + f*x)**5/(5*a**2*c**4*f) - 2*csc(e + f*x)**3/(a**2*c**4*f) + 2*csc(e + f*x)/(a**2*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_29():
    f = 1/((a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**5)
    F = x/(a**2*c**5) + 4*cot(e + f*x)**9/(9*a**2*c**5*f) - cot(e + f*x)**7/(7*a**2*c**5*f) + cot(e + f*x)**5/(5*a**2*c**5*f) - cot(e + f*x)**3/(3*a**2*c**5*f) + cot(e + f*x)/(a**2*c**5*f) + 4*csc(e + f*x)**9/(9*a**2*c**5*f) - 15*csc(e + f*x)**7/(7*a**2*c**5*f) + 21*csc(e + f*x)**5/(5*a**2*c**5*f) - 13*csc(e + f*x)**3/(3*a**2*c**5*f) + 3*csc(e + f*x)/(a**2*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_30():
    f = (-c*sec(e + f*x) + c)**5/(a*sec(e + f*x) + a)**3
    F = c**5*x/a**3 - c**5*tan(e + f*x)/(a**3*f) + 128*c**5*cot(e + f*x)**5/(5*a**3*f) + 128*c**5*cot(e + f*x)**3/(3*a**3*f) + 32*c**5*cot(e + f*x)/(a**3*f) + 8*c**5*atanh(sin(e + f*x))/(a**3*f) - 128*c**5*csc(e + f*x)**5/(5*a**3*f) + 64*c**5*csc(e + f*x)**3/(3*a**3*f) - 16*c**5*csc(e + f*x)/(a**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_31():
    f = (-c*sec(e + f*x) + c)**4/(a*sec(e + f*x) + a)**3
    F = c**4*x/a**3 + c**4*atanh(sin(e + f*x))/(a**3*f) - 23*c**4*tan(e + f*x)/(5*a**3*f*(sec(e + f*x) + 1)) + 14*c**4*tan(e + f*x)/(5*a**3*f*(sec(e + f*x) + 1)**2) - c**4*tan(e + f*x)*sec(e + f*x)**2/(5*a**3*f*(sec(e + f*x) + 1)**3) - 3*c**4*tan(e + f*x)/(a**3*f*(sec(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_32():
    f = (-c*sec(e + f*x) + c)**3/(a*sec(e + f*x) + a)**3
    F = c**3*x/a**3 - 26*c**3*tan(e + f*x)/(15*a**3*f*(sec(e + f*x) + 1)) + 4*c**3*tan(e + f*x)/(15*a**3*f*(sec(e + f*x) + 1)**2) - 8*c**3*tan(e + f*x)/(5*a**3*f*(sec(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_33():
    f = (-c*sec(e + f*x) + c)**2/(a*sec(e + f*x) + a)**3
    F = c**2*x/a**3 - 23*c**2*tan(e + f*x)/(15*a**3*f*(sec(e + f*x) + 1)) - 8*c**2*tan(e + f*x)/(15*a**3*f*(sec(e + f*x) + 1)**2) - 4*c**2*tan(e + f*x)/(5*a**3*f*(sec(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_34():
    f = (-c*sec(e + f*x) + c)/(a*sec(e + f*x) + a)**3
    F = c*x/a**3 - 8*c*tan(e + f*x)/(5*a**3*f*(sec(e + f*x) + 1)) - 3*c*tan(e + f*x)/(5*a**3*f*(sec(e + f*x) + 1)**2) - 2*c*tan(e + f*x)/(5*a**3*f*(sec(e + f*x) + 1)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_35():
    f = 1/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c))
    F = x/(a**3*c) + 2*cot(e + f*x)**5/(5*a**3*c*f) - cot(e + f*x)**3/(3*a**3*c*f) + cot(e + f*x)/(a**3*c*f) - 2*csc(e + f*x)**5/(5*a**3*c*f) + 4*csc(e + f*x)**3/(3*a**3*c*f) - 2*csc(e + f*x)/(a**3*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_36():
    f = 1/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**2)
    F = x/(a**3*c**2) + (1 - sec(e + f*x))*cot(e + f*x)**5/(5*a**3*c**2*f) - (5 - 4*sec(e + f*x))*cot(e + f*x)**3/(15*a**3*c**2*f) + (15 - 8*sec(e + f*x))*cot(e + f*x)/(15*a**3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_37():
    f = 1/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**3)
    F = x/(a**3*c**3) + cot(e + f*x)**5/(5*a**3*c**3*f) - cot(e + f*x)**3/(3*a**3*c**3*f) + cot(e + f*x)/(a**3*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_38():
    f = 1/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**4)
    F = x/(a**3*c**4) - (sec(e + f*x) + 1)*cot(e + f*x)**7/(7*a**3*c**4*f) + (6*sec(e + f*x) + 7)*cot(e + f*x)**5/(35*a**3*c**4*f) + (16*sec(e + f*x) + 35)*cot(e + f*x)/(35*a**3*c**4*f) - (24*sec(e + f*x) + 35)*cot(e + f*x)**3/(105*a**3*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_39():
    f = 1/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**5)
    F = x/(a**3*c**5) + 2*cot(e + f*x)**9/(9*a**3*c**5*f) - cot(e + f*x)**7/(7*a**3*c**5*f) + cot(e + f*x)**5/(5*a**3*c**5*f) - cot(e + f*x)**3/(3*a**3*c**5*f) + cot(e + f*x)/(a**3*c**5*f) + 2*csc(e + f*x)**9/(9*a**3*c**5*f) - 8*csc(e + f*x)**7/(7*a**3*c**5*f) + 12*csc(e + f*x)**5/(5*a**3*c**5*f) - 8*csc(e + f*x)**3/(3*a**3*c**5*f) + 2*csc(e + f*x)/(a**3*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_40():
    f = 1/((a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**6)
    F = x/(a**3*c**6) - 4*cot(e + f*x)**11/(11*a**3*c**6*f) + cot(e + f*x)**9/(9*a**3*c**6*f) - cot(e + f*x)**7/(7*a**3*c**6*f) + cot(e + f*x)**5/(5*a**3*c**6*f) - cot(e + f*x)**3/(3*a**3*c**6*f) + cot(e + f*x)/(a**3*c**6*f) - 4*csc(e + f*x)**11/(11*a**3*c**6*f) + 19*csc(e + f*x)**9/(9*a**3*c**6*f) - 36*csc(e + f*x)**7/(7*a**3*c**6*f) + 34*csc(e + f*x)**5/(5*a**3*c**6*f) - 16*csc(e + f*x)**3/(3*a**3*c**6*f) + 3*csc(e + f*x)/(a**3*c**6*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_41():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**4
    F = 2*sqrt(a)*c**4*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f + 2*a**4*c**4*tan(e + f*x)**7/(7*f*(a*sec(e + f*x) + a)**(sympy.S(7)/2)) - 2*a**3*c**4*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 2*a**2*c**4*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a*c**4*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_42():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**3
    F = 2*sqrt(a)*c**3*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f - 2*a**3*c**3*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 2*a**2*c**3*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a*c**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_43():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**2
    F = 2*sqrt(a)*c**2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f + 2*a**2*c**2*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a*c**2*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_44():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)
    F = 2*sqrt(a)*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f - 2*a*c*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_45():
    f = sqrt(a*sec(e + f*x) + a)/(-c*sec(e + f*x) + c)
    F = 2*sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c*f) + 2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_46():
    f = sqrt(a*sec(e + f*x) + a)/(-c*sec(e + f*x) + c)**2
    F = 2*sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**2*f) + 2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**2*f) - 2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*a*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_47():
    f = sqrt(a*sec(e + f*x) + a)/(-c*sec(e + f*x) + c)**3
    F = 2*sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**3*f) + 2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**3*f) - 2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*a*c**3*f) + 2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*a**2*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_48():
    f = sqrt(a*sec(e + f*x) + a)/(-c*sec(e + f*x) + c)**4
    F = 2*sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**4*f) + 2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**4*f) - 2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*a*c**4*f) + 2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*a**2*c**4*f) - 2*(a*sec(e + f*x) + a)**(sympy.S(7)/2)*cot(e + f*x)**7/(7*a**3*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_49():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**3
    F = 2*a**(sympy.S(3)/2)*c**3*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f - 2*a**5*c**3*tan(e + f*x)**7/(7*f*(a*sec(e + f*x) + a)**(sympy.S(7)/2)) - 2*a**4*c**3*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 2*a**3*c**3*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a**2*c**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_50():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**2
    F = 2*a**(sympy.S(3)/2)*c**2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f + 2*a**4*c**2*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 2*a**3*c**2*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a**2*c**2*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_51():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)
    F = 2*a**(sympy.S(3)/2)*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f - 2*a**3*c*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a**2*c*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_52():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(-c*sec(e + f*x) + c)
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c*f) + 4*a*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_53():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(-c*sec(e + f*x) + c)**2
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**2*f) + 2*a*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**2*f) - 4*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_54():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(-c*sec(e + f*x) + c)**3
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**3*f) + 2*a*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**3*f) - 2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*c**3*f) + 4*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*a*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_55():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(-c*sec(e + f*x) + c)**4
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**4*f) + 2*a*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**4*f) - 2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*c**4*f) + 2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*a*c**4*f) - 4*(a*sec(e + f*x) + a)**(sympy.S(7)/2)*cot(e + f*x)**7/(7*a**2*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_56():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**3
    F = 2*a**(sympy.S(5)/2)*c**3*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f - 2*a**7*c**3*tan(e + f*x)**9/(9*f*(a*sec(e + f*x) + a)**(sympy.S(9)/2)) - 6*a**6*c**3*tan(e + f*x)**7/(7*f*(a*sec(e + f*x) + a)**(sympy.S(7)/2)) - 2*a**5*c**3*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 2*a**4*c**3*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a**3*c**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_57():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**2
    F = 2*a**(sympy.S(5)/2)*c**2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f + 2*a**6*c**2*tan(e + f*x)**7/(7*f*(a*sec(e + f*x) + a)**(sympy.S(7)/2)) + 6*a**5*c**2*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 2*a**4*c**2*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a**3*c**2*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_58():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)
    F = 2*a**(sympy.S(5)/2)*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f - 2*a**5*c*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) - 2*a**4*c*tan(e + f*x)**3/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 2*a**3*c*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_59():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c*f) - 2*a**3*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)) + 8*a**2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_60():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**2
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**2*f) - 8*a*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_61():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**3
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**3*f) + 2*a**2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**3*f) + 8*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_62():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**4
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**4*f) + 2*a**2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**4*f) - 2*a*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*c**4*f) - 8*(a*sec(e + f*x) + a)**(sympy.S(7)/2)*cot(e + f*x)**7/(7*a*c**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_63():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**5
    F = 2*a**(sympy.S(5)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c**5*f) + 2*a**2*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(c**5*f) - 2*a*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*c**5*f) + 2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*c**5*f) + 8*(a*sec(e + f*x) + a)**(sympy.S(9)/2)*cot(e + f*x)**9/(9*a**2*c**5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_64():
    f = (-c*sec(e + f*x) + c)**4/sqrt(a*sec(e + f*x) + a)
    F = 2*a**2*c**4*tan(e + f*x)**5/(5*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) - 2*a*c**4*tan(e + f*x)**3/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + 14*c**4*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) + 2*c**4*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*f) - 16*sqrt(2)*c**4*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_65():
    f = (-c*sec(e + f*x) + c)**3/sqrt(a*sec(e + f*x) + a)
    F = -2*a*c**3*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + 6*c**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) + 2*c**3*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*f) - 8*sqrt(2)*c**3*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_66():
    f = (-c*sec(e + f*x) + c)**2/sqrt(a*sec(e + f*x) + a)
    F = 2*c**2*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) + 2*c**2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*f) - 4*sqrt(2)*c**2*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_67():
    f = (-c*sec(e + f*x) + c)/sqrt(a*sec(e + f*x) + a)
    F = 2*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*f) - 2*sqrt(2)*c*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_68():
    f = 1/(sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c))
    F = sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(a*c*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*c*f) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(2*sqrt(a)*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_69():
    f = 1/(sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**2)
    F = 3*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(2*a*c**2*f) - (a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(3*a**2*c**2*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*c**2*f) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(4*sqrt(a)*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_70():
    f = 1/(sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**3)
    F = 7*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(4*a*c**3*f) - (a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(2*a**2*c**3*f) + (a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(5*a**3*c**3*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*c**3*f) - sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(8*sqrt(a)*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_71():
    f = (-c*sec(e + f*x) + c)**4/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -a*c**4*sin(e + f*x)*tan(e + f*x)**4*sec(e/2 + f*x/2)**2/(f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) + 8*c**4*tan(e + f*x)**3/(3*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 14*c**4*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)) + 2*c**4*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(3)/2)*f) + 12*sqrt(2)*c**4*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_72():
    f = (-c*sec(e + f*x) + c)**3/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = c**3*sin(e + f*x)*tan(e + f*x)**2*sec(e/2 + f*x/2)**2/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - 4*c**3*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)) + 2*c**3*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(3)/2)*f) + 2*sqrt(2)*c**3*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_73():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c))
    F = sqrt(a*sec(e + f*x) + a)*cos(e + f*x)*cot(e + f*x)*sec(e/2 + f*x/2)**2/(4*a**2*c*f) + sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(4*a**2*c*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(3)/2)*c*f) - 7*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(8*a**(sympy.S(3)/2)*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_74():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**2)
    F = 7*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(8*a**2*c**2*f) - (a*sec(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)*cot(e + f*x)**3*sec(e/2 + f*x/2)**2/(4*a**3*c**2*f) + (a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(12*a**3*c**2*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(3)/2)*c**2*f) - 9*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(16*a**(sympy.S(3)/2)*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_75():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**3)
    F = 21*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(16*a**2*c**3*f) - 5*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(24*a**3*c**3*f) + (a*sec(e + f*x) + a)**(sympy.S(5)/2)*cos(e + f*x)*cot(e + f*x)**5*sec(e/2 + f*x/2)**2/(4*a**4*c**3*f) - 3*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*cot(e + f*x)**5/(20*a**4*c**3*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(3)/2)*c**3*f) - 11*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(32*a**(sympy.S(3)/2)*c**3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_76():
    f = (-c*sec(e + f*x) + c)**5/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = a*c**5*sin(e + f*x)**2*tan(e + f*x)**5*sec(e/2 + f*x/2)**4/(4*f*(a*sec(e + f*x) + a)**(sympy.S(7)/2)) + 3*c**5*sin(e + f*x)*tan(e + f*x)**4*sec(e/2 + f*x/2)**2/(4*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) - 19*c**5*tan(e + f*x)**3/(6*a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + 21*c**5*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)) + 2*c**5*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(5)/2)*f) - 23*sqrt(2)*c**5*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_77():
    f = (-c*sec(e + f*x) + c)**4/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -c**4*sin(e + f*x)**2*tan(e + f*x)**3*sec(e/2 + f*x/2)**4/(4*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) - c**4*sin(e + f*x)*tan(e + f*x)**2*sec(e/2 + f*x/2)**2/(4*a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + 7*c**4*tan(e + f*x)/(2*a**2*f*sqrt(a*sec(e + f*x) + a)) + 2*c**4*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(5)/2)*f) - 11*sqrt(2)*c**4*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(2*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_78():
    f = (-c*sec(e + f*x) + c)**3/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = c**3*sin(e + f*x)**2*tan(e + f*x)*sec(e/2 + f*x/2)**4/(4*a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) - c**3*sin(e + f*x)*sec(e/2 + f*x/2)**2/(4*a**2*f*sqrt(a*sec(e + f*x) + a)) + 2*c**3*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(5)/2)*f) - 7*sqrt(2)*c**3*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(4*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_79():
    f = (-c*sec(e + f*x) + c)**2/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -c**2*sin(e + f*x)*cos(e + f*x)*sec(e/2 + f*x/2)**4/(4*a**2*f*sqrt(a*sec(e + f*x) + a)) - 3*c**2*sin(e + f*x)*sec(e/2 + f*x/2)**2/(8*a**2*f*sqrt(a*sec(e + f*x) + a)) + 2*c**2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(5)/2)*f) - 11*sqrt(2)*c**2*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(8*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_80():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c))
    F = sqrt(a*sec(e + f*x) + a)*cos(e + f*x)**2*cot(e + f*x)*sec(e/2 + f*x/2)**4/(16*a**3*c*f) + 13*sqrt(a*sec(e + f*x) + a)*cos(e + f*x)*cot(e + f*x)*sec(e/2 + f*x/2)**2/(32*a**3*c*f) - 7*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(32*a**3*c*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(5)/2)*c*f) - 71*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(64*a**(sympy.S(5)/2)*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_81():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**2)
    F = 21*sqrt(a*sec(e + f*x) + a)*cot(e + f*x)/(64*a**3*c**2*f) - (a*sec(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)**2*cot(e + f*x)**3*sec(e/2 + f*x/2)**4/(16*a**4*c**2*f) - 15*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cos(e + f*x)*cot(e + f*x)**3*sec(e/2 + f*x/2)**2/(32*a**4*c**2*f) + 43*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*cot(e + f*x)**3/(96*a**4*c**2*f) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(5)/2)*c**2*f) - 107*sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(128*a**(sympy.S(5)/2)*c**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_82():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)
    F = a*c**4*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a*c**3*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) - a*c**2*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a)) - a*c*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_83():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = a*c**3*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a*c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) - a*c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_84():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = a*c**2*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a*c*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_85():
    f = sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)
    F = a*c*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_86():
    f = sqrt(a*sec(e + f*x) + a)/sqrt(-c*sec(e + f*x) + c)
    F = a*log(1 - cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_87():
    f = sqrt(a*sec(e + f*x) + a)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -a*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a*log(1 - cos(e + f*x))*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_88():
    f = sqrt(a*sec(e + f*x) + a)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = -a*tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - a*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a*log(1 - cos(e + f*x))*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_89():
    f = sqrt(a*sec(e + f*x) + a)/(-c*sec(e + f*x) + c)**(sympy.S(7)/2)
    F = -a*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)) - a*tan(e + f*x)/(2*c*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - a*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a*log(1 - cos(e + f*x))*tan(e + f*x)/(c**3*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_90():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = a**2*c**3*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a**2*c**2*sqrt(-c*sec(e + f*x) + c)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) - a**2*c*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(a*sec(e + f*x) + a)) + a**2*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_91():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = a**2*c**2*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + a**2*c**2*tan(e + f*x)**3/(2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_92():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)
    F = a**2*c*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a*c*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_93():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/sqrt(-c*sec(e + f*x) + c)
    F = 2*a**2*log(1 - sec(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + a**2*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_94():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -2*a**2*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a**2*log(1 - cos(e + f*x))*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_95():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = -a**2*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - a**2*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a**2*log(1 - cos(e + f*x))*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_96():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(-c*sec(e + f*x) + c)**(sympy.S(7)/2)
    F = -2*a**2*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)) - a**2*tan(e + f*x)/(2*c*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - a**2*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a**2*log(1 - cos(e + f*x))*tan(e + f*x)/(c**3*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_97():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = a**3*c**3*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a**3*c**3*tan(e + f*x)**5/(4*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + a**3*c**3*tan(e + f*x)**3/(2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_98():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = a**3*c**2*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a**2*c**2*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c)) - a*c**2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(-c*sec(e + f*x) + c)) + c**2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*tan(e + f*x)/(3*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_99():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)
    F = a**3*c*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - a**2*c*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(f*sqrt(-c*sec(e + f*x) + c)) - a*c*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(2*f*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_100():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/sqrt(-c*sec(e + f*x) + c)
    F = 4*a**3*log(1 - sec(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + a**3*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + a**3*tan(e + f*x)*sec(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_101():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**(sympy.S(3)/2)
    F = -4*a**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a**3*log(cos(e + f*x))*tan(e + f*x)/(c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_102():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**(sympy.S(5)/2)
    F = -2*a**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) + a**3*log(1 - cos(e + f*x))*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_103():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**(sympy.S(7)/2)
    F = -4*a**3*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)) - a**3*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a**3*log(1 - cos(e + f*x))*tan(e + f*x)/(c**3*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_104():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**(sympy.S(9)/2)
    F = -a**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(9)/2)) - a**3*tan(e + f*x)/(2*c**2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - a**3*tan(e + f*x)/(c**3*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a**3*log(1 - cos(e + f*x))*tan(e + f*x)/(c**4*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_105():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(-c*sec(e + f*x) + c)**(sympy.S(11)/2)
    F = -4*a**3*tan(e + f*x)/(5*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(11)/2)) - a**3*tan(e + f*x)/(3*c**2*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(7)/2)) - a**3*tan(e + f*x)/(2*c**3*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2)) - a**3*tan(e + f*x)/(c**4*f*sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2)) + a**3*log(1 - cos(e + f*x))*tan(e + f*x)/(c**5*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_106():
    f = (-c*sec(e + f*x) + c)**(sympy.S(7)/2)/sqrt(a*sec(e + f*x) + a)
    F = 8*c**4*log(sec(e + f*x) + 1)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + c**4*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + c**4*tan(e + f*x)*sec(e + f*x)**2/(2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - 4*c**4*tan(e + f*x)*sec(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_107():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)/sqrt(a*sec(e + f*x) + a)
    F = 4*c**3*log(sec(e + f*x) + 1)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + c**3*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - c**3*tan(e + f*x)*sec(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_108():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)/sqrt(a*sec(e + f*x) + a)
    F = 2*c**2*log(sec(e + f*x) + 1)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + c**2*log(cos(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_109():
    f = sqrt(-c*sec(e + f*x) + c)/sqrt(a*sec(e + f*x) + a)
    F = c*log(cos(e + f*x) + 1)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_110():
    f = 1/(sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    F = log(sin(e + f*x))*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_111():
    f = 1/(sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = 7*log(1 - sec(e + f*x))*tan(e + f*x)/(8*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + log(sec(e + f*x) + 1)*tan(e + f*x)/(8*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + log(cos(e + f*x))*tan(e + f*x)/(c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - 3*tan(e + f*x)/(4*c**2*f*(1 - sec(e + f*x))*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - tan(e + f*x)/(4*c**2*f*(1 - sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_112():
    f = (-c*sec(e + f*x) + c)**(sympy.S(7)/2)/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -4*c**4*log(sec(e + f*x) + 1)*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + c**4*log(cos(e + f*x))*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + c**4*tan(e + f*x)*sec(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - 8*c**4*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_113():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -4*c**3*tan(e + f*x)/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)) + c**3*log(cos(e + f*x))*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_114():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -2*c**2*tan(e + f*x)/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)) + c**2*log(cos(e + f*x) + 1)*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_115():
    f = sqrt(-c*sec(e + f*x) + c)/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -c*tan(e + f*x)/(f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)) + c*log(cos(e + f*x) + 1)*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_116():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c))
    F = log(1 - sec(e + f*x))*tan(e + f*x)/(4*a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 3*log(sec(e + f*x) + 1)*tan(e + f*x)/(4*a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + log(cos(e + f*x))*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - tan(e + f*x)/(2*a*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_117():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = log(sin(e + f*x))*tan(e + f*x)/(a*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + cot(e + f*x)/(2*a*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_118():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = 11*log(1 - sec(e + f*x))*tan(e + f*x)/(16*a*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 5*log(sec(e + f*x) + 1)*tan(e + f*x)/(16*a*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + log(cos(e + f*x))*tan(e + f*x)/(a*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - tan(e + f*x)/(8*a*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)) - tan(e + f*x)/(2*a*c**2*f*(1 - sec(e + f*x))*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - tan(e + f*x)/(8*a*c**2*f*(1 - sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_119():
    f = (-c*sec(e + f*x) + c)**(sympy.S(7)/2)/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = 2*c**4*log(sec(e + f*x) + 1)*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + c**4*log(cos(e + f*x))*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 4*c**4*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)) - 4*c**4*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_120():
    f = (-c*sec(e + f*x) + c)**(sympy.S(5)/2)/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -2*c**3*tan(e + f*x)/(f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)) + c**3*log(cos(e + f*x) + 1)*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_121():
    f = (-c*sec(e + f*x) + c)**(sympy.S(3)/2)/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -c**2*tan(e + f*x)/(f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)) - c**2*tan(e + f*x)/(a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)) + c**2*log(cos(e + f*x) + 1)*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_122():
    f = sqrt(-c*sec(e + f*x) + c)/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -c*tan(e + f*x)/(2*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c)) - c*tan(e + f*x)/(a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*sqrt(-c*sec(e + f*x) + c)) + c*log(cos(e + f*x) + 1)*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_123():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*sqrt(-c*sec(e + f*x) + c))
    F = log(1 - sec(e + f*x))*tan(e + f*x)/(8*a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 7*log(sec(e + f*x) + 1)*tan(e + f*x)/(8*a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + log(cos(e + f*x))*tan(e + f*x)/(a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - 3*tan(e + f*x)/(4*a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)) - tan(e + f*x)/(4*a**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_124():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(3)/2))
    F = 5*log(1 - sec(e + f*x))*tan(e + f*x)/(16*a**2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + 11*log(sec(e + f*x) + 1)*tan(e + f*x)/(16*a**2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + log(cos(e + f*x))*tan(e + f*x)/(a**2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - tan(e + f*x)/(2*a**2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)) - tan(e + f*x)/(8*a**2*c*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)*(sec(e + f*x) + 1)**2) - tan(e + f*x)/(8*a**2*c*f*(1 - sec(e + f*x))*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_125():
    f = 1/((a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**(sympy.S(5)/2))
    F = log(sin(e + f*x))*tan(e + f*x)/(a**2*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) - cot(e + f*x)**3/(4*a**2*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c)) + cot(e + f*x)/(2*a**2*c**2*f*sqrt(a*sec(e + f*x) + a)*sqrt(-c*sec(e + f*x) + c))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_126():
    f = (-c*sec(e + f*x) + c)**n*(sec(e + f*x) + 1)**m
    F = 2**(m + sympy.S.Half)*(-c*sec(e + f*x) + c)**n*tan(e + f*x)*appellf1(n + sympy.S.Half, 1, sympy.S.Half - m, n + sympy.S(3)/2, 1 - sec(e + f*x), sympy.S.Half - sec(e + f*x)/2)/(f*(2*n + 1)*sqrt(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_127():
    f = (a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**n
    F = 2**(n + sympy.S.Half)*c*(1 - sec(e + f*x))**(sympy.S.Half - n)*(a*sec(e + f*x) + a)**m*(-c*sec(e + f*x) + c)**(n - 1)*tan(e + f*x)*appellf1(m + sympy.S.Half, 1, sympy.S.Half - n, m + sympy.S(3)/2, sec(e + f*x) + 1, sec(e + f*x)/2 + sympy.S.Half)/(f*(2*m + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_128():
    f = (a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**n
    F = 2**(n + sympy.S.Half)*c*(1 - sec(e + f*x))**(sympy.S.Half - n)*(a*sec(e + f*x) + a)**3*(-c*sec(e + f*x) + c)**(n - 1)*tan(e + f*x)*appellf1(sympy.S(7)/2, 1, sympy.S.Half - n, sympy.S(9)/2, sec(e + f*x) + 1, sec(e + f*x)/2 + sympy.S.Half)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_129():
    f = (a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**n
    F = 2**(n + sympy.S.Half)*c*(1 - sec(e + f*x))**(sympy.S.Half - n)*(a*sec(e + f*x) + a)**2*(-c*sec(e + f*x) + c)**(n - 1)*tan(e + f*x)*appellf1(sympy.S(5)/2, 1, sympy.S.Half - n, sympy.S(7)/2, sec(e + f*x) + 1, sec(e + f*x)/2 + sympy.S.Half)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_130():
    f = (a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**n
    F = 2**(n + sympy.S.Half)*c*(1 - sec(e + f*x))**(sympy.S.Half - n)*(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**(n - 1)*tan(e + f*x)*appellf1(sympy.S(3)/2, 1, sympy.S.Half - n, sympy.S(5)/2, sec(e + f*x) + 1, sec(e + f*x)/2 + sympy.S.Half)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_131():
    f = (-c*sec(e + f*x) + c)**n/(a*sec(e + f*x) + a)
    F = -2**(n + sympy.S.Half)*c*(1 - sec(e + f*x))**(sympy.S.Half - n)*(-c*sec(e + f*x) + c)**(n - 1)*tan(e + f*x)*appellf1(sympy.S(-1)/2, 1, sympy.S.Half - n, sympy.S.Half, sec(e + f*x) + 1, sec(e + f*x)/2 + sympy.S.Half)/(f*(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_132():
    f = (-c*sec(e + f*x) + c)**n/(a*sec(e + f*x) + a)**2
    F = -2**(n + sympy.S.Half)*c*(1 - sec(e + f*x))**(sympy.S.Half - n)*(-c*sec(e + f*x) + c)**(n - 1)*tan(e + f*x)*appellf1(sympy.S(-3)/2, 1, sympy.S.Half - n, sympy.S(-1)/2, sec(e + f*x) + 1, sec(e + f*x)/2 + sympy.S.Half)/(3*f*(a*sec(e + f*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_133():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)*(-c*sec(e + f*x) + c)**n
    F = 2*a**3*(-c*sec(e + f*x) + c)**n*tan(e + f*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) + 6*a**3*(-c*sec(e + f*x) + c)**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) - 2*a**3*(-c*sec(e + f*x) + c)**(n + 1)*tan(e + f*x)/(c*f*(2*n + 3)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_134():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)*(-c*sec(e + f*x) + c)**n
    F = 2*a**2*(-c*sec(e + f*x) + c)**n*tan(e + f*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(-c*sec(e + f*x) + c)**n*tan(e + f*x)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_135():
    f = sqrt(a*sec(e + f*x) + a)*(-c*sec(e + f*x) + c)**n
    F = 2*a*(-c*sec(e + f*x) + c)**n*tan(e + f*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_136():
    f = (-c*sec(e + f*x) + c)**n/sqrt(a*sec(e + f*x) + a)
    F = -(-c*sec(e + f*x) + c)**n*tan(e + f*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), sympy.S.Half - sec(e + f*x)/2)/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) + 2*(-c*sec(e + f*x) + c)**n*tan(e + f*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), 1 - sec(e + f*x))/(f*(2*n + 1)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_137():
    f = (-c*sec(e + f*x) + c)**n/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -(5 - 2*n)*(-c*sec(e + f*x) + c)**n*tan(e + f*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), sympy.S.Half - sec(e + f*x)/2)/(4*a*f*(2*n + 1)*sqrt(a*sec(e + f*x) + a)) - (-c*sec(e + f*x) + c)**n*tan(e + f*x)/(2*a*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) + 2*(-c*sec(e + f*x) + c)**n*tan(e + f*x)*hyper((1, n + sympy.S.Half), (n + sympy.S(3)/2,), 1 - sec(e + f*x))/(a*f*(2*n + 1)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_138():
    f = sqrt(a*sec(e + f*x) + a)/(c*sec(e + f*x) + c)
    F = 2*sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c*f) - sqrt(2)*sqrt(a)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_139():
    f = (c + d*sec(e + f*x))**(sympy.S(3)/2)/(a*sec(e + f*x) + a)
    F = (Integer(-1) * ((Integer(2) * Symbol('c') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('c') * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('c') + Symbol('d'))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), ((Symbol('c') + (Integer(-1) * Symbol('d'))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('d') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('d') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.elliptic_e(sympy.asin((sympy.tan((Symbol('e') + (Symbol('f') * x))) * ((Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1)))), ((Symbol('c') + (Integer(-1) * Symbol('d'))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + Symbol('d')) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_140():
    f = sqrt(c + d*sec(e + f*x))/(a*sec(e + f*x) + a)
    F = (Integer(-1) * ((Integer(2) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('c') * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('c') + Symbol('d'))) * (sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), ((Symbol('c') + (Integer(-1) * Symbol('d'))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('d') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))))) * sympy.sqrt(((Symbol('d') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1)))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((sympy.elliptic_e(sympy.asin((sympy.tan((Symbol('e') + (Symbol('f') * x))) * ((Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1)))), ((Symbol('c') + (Integer(-1) * Symbol('d'))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + Symbol('d')) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_141():
    f = 1/(sqrt(c + d*sec(e + f*x))*(a*sec(e + f*x) + a))
    F = ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('c') + Symbol('d'))))**(Integer(-1)))), ((Symbol('c') + Symbol('d')) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('d') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('d') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Symbol('a') * (Symbol('c') + (Integer(-1) * Symbol('d'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('c') + Symbol('d')) * (Symbol('c'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('c') + Symbol('d'))))**(Integer(-1)))), ((Symbol('c') + Symbol('d')) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('d') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('d') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('c') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((sympy.elliptic_e(sympy.asin((sympy.tan((Symbol('e') + (Symbol('f') * x))) * ((Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1)))), ((Symbol('c') + (Integer(-1) * Symbol('d'))) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt(((Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Integer(-1))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * (Symbol('c') + (Integer(-1) * Symbol('d'))) * Symbol('f') * sympy.sqrt(((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + Symbol('d')) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(-1))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_142():
    f = (c + d*sec(e + f*x))**4*sqrt(a*sec(e + f*x) + a)
    F = 2*a**(sympy.S(3)/2)*c**4*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a*d*(2*c + d)*(2*c**2 + 2*c*d + d**2)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) - 2*d**2*(-a*sec(e + f*x) + a)*(6*c**2 + 8*c*d + 3*d**2)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)) + 2*d**3*(4*c + 3*d)*(-a*sec(e + f*x) + a)**2*tan(e + f*x)/(5*a*f*sqrt(a*sec(e + f*x) + a)) - 2*d**4*(-a*sec(e + f*x) + a)**3*tan(e + f*x)/(7*a**2*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_143():
    f = (c + d*sec(e + f*x))**3*sqrt(a*sec(e + f*x) + a)
    F = 2*a**(sympy.S(3)/2)*c**3*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a*d*(3*c**2 + 3*c*d + d**2)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) - 2*d**2*(3*c + 2*d)*(-a*sec(e + f*x) + a)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)) + 2*d**3*(-a*sec(e + f*x) + a)**2*tan(e + f*x)/(5*a*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_144():
    f = (c + d*sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)
    F = 2*a**(sympy.S(3)/2)*c**2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a*d*(2*c + d)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) - 2*d**2*(-a*sec(e + f*x) + a)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_145():
    f = (c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)
    F = 2*sqrt(a)*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f + 2*a*d*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_146():
    f = sqrt(a*sec(e + f*x) + a)/(c + d*sec(e + f*x))
    F = -2*sqrt(a)*sqrt(d)*atan(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(c*f*sqrt(c + d)) + 2*sqrt(a)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_147():
    f = sqrt(a*sec(e + f*x) + a)/(c + d*sec(e + f*x))**2
    F = -a**(sympy.S(3)/2)*sqrt(d)*(3*c + 2*d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**2*f*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**(sympy.S(3)/2)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**2*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - a*d*tan(e + f*x)/(c*f*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_148():
    f = sqrt(a*sec(e + f*x) + a)/(c + d*sec(e + f*x))**3
    F = -a**(sympy.S(3)/2)*sqrt(d)*(15*c**2 + 20*c*d + 8*d**2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(4*c**3*f*(c + d)**(sympy.S(5)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**(sympy.S(3)/2)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**3*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - a*d*tan(e + f*x)/(2*c*f*(c + d)*(c + d*sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)) - a*d*(7*c + 4*d)*tan(e + f*x)/(4*c**2*f*(c + d)**2*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_149():
    f = (c + d*sec(e + f*x))**3*(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(5)/2)*c**3*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(c + d*sec(e + f*x))**3*tan(e + f*x)/(7*f*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(c + d*sec(e + f*x))**2*(6*c + 13*d)*tan(e + f*x)/(35*f*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(72*c**3 + 486*c**2*d + 378*c*d**2 + 104*d**3 + d*(24*c**2 + 111*c*d + 52*d**2)*sec(e + f*x))*tan(e + f*x)/(105*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_150():
    f = (c + d*sec(e + f*x))**2*(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(5)/2)*c**2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(c + d*sec(e + f*x))**2*tan(e + f*x)/(5*f*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(12*c**2 + 50*c*d + 18*d**2 + d*(4*c + 9*d)*sec(e + f*x))*tan(e + f*x)/(15*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_151():
    f = (c + d*sec(e + f*x))*(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*a**(sympy.S(3)/2)*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f + 2*a**2*(3*c + 4*d)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)) + 2*a*d*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_152():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sec(e + f*x))
    F = 2*a**(sympy.S(3)/2)*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(c*f) + 2*a**(sympy.S(3)/2)*(c - d)*atan(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(c*sqrt(d)*f*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_153():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sec(e + f*x))**2
    F = 2*a**(sympy.S(5)/2)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**2*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + a**(sympy.S(5)/2)*(c**2 - 3*c*d - 2*d**2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**2*sqrt(d)*f*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + a**2*(c - d)*tan(e + f*x)/(c*f*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_154():
    f = (a*sec(e + f*x) + a)**(sympy.S(3)/2)/(c + d*sec(e + f*x))**3
    F = 2*a**(sympy.S(5)/2)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**3*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + a**(sympy.S(5)/2)*(3*c**3 - 15*c**2*d - 20*c*d**2 - 8*d**3)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(4*c**3*sqrt(d)*f*(c + d)**(sympy.S(5)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + a**2*(c - d)*tan(e + f*x)/(2*c*f*(c + d)*(c + d*sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)) + a**2*(3*c**2 - 7*c*d - 4*d**2)*tan(e + f*x)/(4*c**2*f*(c + d)**2*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_155():
    f = (c + d*sec(e + f*x))**3*(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(7)/2)*c**3*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**3*(3*c**3 + 12*c**2*d + 12*c*d**2 + 4*d**3)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) + 2*a*d*(-a*sec(e + f*x) + a)**2*(3*c**2 + 15*c*d + 13*d**2)*tan(e + f*x)/(5*f*sqrt(a*sec(e + f*x) + a)) - 6*d**2*(c + 2*d)*(-a*sec(e + f*x) + a)**3*tan(e + f*x)/(7*f*sqrt(a*sec(e + f*x) + a)) - (-a**3*sec(e + f*x) + a**3)*(2*c**3 + 24*c**2*d + 48*c*d**2 + 24*d**3)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)) + 2*d**3*(-a*sec(e + f*x) + a)**4*tan(e + f*x)/(9*a*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_156():
    f = (c + d*sec(e + f*x))**2*(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(7)/2)*c**2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**3*(c + 2*d)*(3*c + 2*d)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) + 2*a*d*(2*c + 5*d)*(-a*sec(e + f*x) + a)**2*tan(e + f*x)/(5*f*sqrt(a*sec(e + f*x) + a)) - 2*d**2*(-a*sec(e + f*x) + a)**3*tan(e + f*x)/(7*f*sqrt(a*sec(e + f*x) + a)) - (-a**3*sec(e + f*x) + a**3)*(2*c**2 + 16*c*d + 16*d**2)*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_157():
    f = (c + d*sec(e + f*x))*(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = 2*a**(sympy.S(5)/2)*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/f + 2*a**3*(35*c + 32*d)*tan(e + f*x)/(15*f*sqrt(a*sec(e + f*x) + a)) + 2*a**2*(5*c + 8*d)*sqrt(a*sec(e + f*x) + a)*tan(e + f*x)/(15*f) + 2*a*d*(a*sec(e + f*x) + a)**(sympy.S(3)/2)*tan(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_158():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sec(e + f*x))
    F = 2*a**(sympy.S(7)/2)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 2*a**(sympy.S(7)/2)*(c - d)**2*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c*d**(sympy.S(3)/2)*f*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**3*tan(e + f*x)/(d*f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_159():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sec(e + f*x))**2
    F = -a**(sympy.S(7)/2)*(c - d)**2*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c*d**(sympy.S(3)/2)*f*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**(sympy.S(7)/2)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**2*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**(sympy.S(7)/2)*(c - d)*sqrt(c + d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**2*d**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - a**3*(c - d)**2*tan(e + f*x)/(c*d*f*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_160():
    f = (a*sec(e + f*x) + a)**(sympy.S(5)/2)/(c + d*sec(e + f*x))**3
    F = -3*a**(sympy.S(7)/2)*(c - d)**2*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(4*c*d**(sympy.S(3)/2)*f*(c + d)**(sympy.S(5)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + a**(sympy.S(7)/2)*(c - d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**2*d**(sympy.S(3)/2)*f*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 2*a**(sympy.S(7)/2)*sqrt(d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**3*f*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*a**(sympy.S(7)/2)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**3*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - a**3*(c - d)**2*tan(e + f*x)/(2*c*d*f*(c + d)*(c + d*sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)) - 3*a**3*(c - d)**2*tan(e + f*x)/(4*c*d*f*(c + d)**2*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)) + a**3*(c - d)*tan(e + f*x)/(c**2*d*f*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_161():
    f = (c + d*sec(e + f*x))**3/sqrt(a*sec(e + f*x) + a)
    F = 2*sqrt(a)*c**3*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*sqrt(a)*(c - d)**3*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 2*d**3*(1 - sec(e + f*x))*tan(e + f*x)/(3*f*sqrt(a*sec(e + f*x) + a)) + 2*d**3*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a)) + d**2*(6*c - 2*d)*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_162():
    f = (c + d*sec(e + f*x))**2/sqrt(a*sec(e + f*x) + a)
    F = 2*sqrt(a)*c**2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*sqrt(a)*(c - d)**2*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*d**2*tan(e + f*x)/(f*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_163():
    f = (c + d*sec(e + f*x))/sqrt(a*sec(e + f*x) + a)
    F = 2*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*f) - sqrt(2)*(c - d)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_164():
    f = 1/((c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*(c - d)) + 2*d**(sympy.S(3)/2)*atan(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d)*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*c*f*(c - d)*sqrt(c + d)) + 2*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(sqrt(a)*c*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_165():
    f = 1/((c + d*sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a))
    F = -sqrt(2)*sqrt(a)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(f*(c - d)**2*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + sqrt(a)*d**(sympy.S(3)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c*f*(c - d)*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*sqrt(a)*d**(sympy.S(3)/2)*(2*c - d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**2*f*(c - d)**2*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*sqrt(a)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**2*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + d**2*tan(e + f*x)/(c*f*(c + d*sec(e + f*x))*(c**2 - d**2)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_166():
    f = 1/((c + d*sec(e + f*x))**3*sqrt(a*sec(e + f*x) + a))
    F = -sqrt(2)*sqrt(a)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(f*(c - d)**3*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 3*sqrt(a)*d**(sympy.S(3)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(4*c*f*(c - d)*(c + d)**(sympy.S(5)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + sqrt(a)*d**(sympy.S(3)/2)*(2*c - d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**2*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*sqrt(a)*d**(sympy.S(3)/2)*(3*c**2 - 3*c*d + d**2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(c**3*f*(c - d)**3*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*sqrt(a)*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(c**3*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + d**2*tan(e + f*x)/(2*c*f*(c + d*sec(e + f*x))**2*(c**2 - d**2)*sqrt(a*sec(e + f*x) + a)) + 3*d**2*tan(e + f*x)/(4*c*f*(c - d)*(c + d)**2*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)) + d**2*(2*c - d)*tan(e + f*x)/(c**2*f*(c - d)**2*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_167():
    f = (c + d*sec(e + f*x))**3/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = 2*d**3*tan(e + f*x)/(a*f*sqrt(a*sec(e + f*x) + a)) - (c - d)**3*tan(e + f*x)/(2*a*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) + 2*c**3*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(sqrt(a)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - d)**3*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*sqrt(a)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - d)**2*(c + 2*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(sqrt(a)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_168():
    f = (c + d*sec(e + f*x))**2/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c - d)**2*tan(e + f*x)/(2*a*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) + 2*c**2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(sqrt(a)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - d)**2*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*sqrt(a)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c**2 - d**2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(sqrt(a)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_169():
    f = (c + d*sec(e + f*x))/(a*sec(e + f*x) + a)**(sympy.S(3)/2)
    F = -(c - d)*tan(e + f*x)/(2*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + 2*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(3)/2)*f) - sqrt(2)*(5*c - d)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(4*a**(sympy.S(3)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_170():
    f = 1/((c + d*sec(e + f*x))*(a*sec(e + f*x) + a)**(sympy.S(3)/2))
    F = -tan(e + f*x)/(2*a*f*(c - d)*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - sqrt(2)*(c - 2*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(sqrt(a)*f*(c - d)**2*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*sqrt(a)*f*(c - d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 2*d**(sympy.S(5)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(sqrt(a)*c*f*(c - d)**2*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(sqrt(a)*c*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_171():
    f = 1/((c + d*sec(e + f*x))**2*(a*sec(e + f*x) + a)**(sympy.S(3)/2))
    F = -tan(e + f*x)/(2*a*f*(c - d)**2*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - d**3*tan(e + f*x)/(a*c*f*(c - d)**2*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - 3*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(sqrt(a)*f*(c - d)**3*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*sqrt(a)*f*(c - d)**2*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - d**(sympy.S(5)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(sqrt(a)*c*f*(c - d)**2*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - d**(sympy.S(5)/2)*(6*c - 2*d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(sqrt(a)*c**2*f*(c - d)**3*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(sqrt(a)*c**2*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_172():
    f = 1/((c + d*sec(e + f*x))**3*(a*sec(e + f*x) + a)**(sympy.S(3)/2))
    F = -tan(e + f*x)/(2*a*f*(c - d)**3*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - 3*d**3*tan(e + f*x)/(4*a*c*f*(c + d*sec(e + f*x))*(c**2 - d**2)**2*sqrt(a*sec(e + f*x) + a)) - d**3*tan(e + f*x)/(2*a*c*f*(c - d)**2*(c + d)*(c + d*sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)) - d**3*(3*c - d)*tan(e + f*x)/(a*c**2*f*(c - d)**3*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - 4*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(sqrt(a)*f*(c - d)**4*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*sqrt(a)*f*(c - d)**3*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 3*d**(sympy.S(5)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(4*sqrt(a)*c*f*(c - d)**2*(c + d)**(sympy.S(5)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - d**(sympy.S(5)/2)*(3*c - d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(sqrt(a)*c**2*f*(c - d)**3*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 2*d**(sympy.S(5)/2)*(6*c**2 - 4*c*d + d**2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(sqrt(a)*c**3*f*(c - d)**4*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(sqrt(a)*c**3*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_173():
    f = (c + d*sec(e + f*x))**3/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -3*(c - d)**3*tan(e + f*x)/(16*a**2*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - (c - d)**3*tan(e + f*x)/(4*a**2*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)**2) - (c - d)**2*(c + 2*d)*tan(e + f*x)/(2*a**2*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) + 2*c**3*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 3*sqrt(2)*(c - d)**3*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - d)**2*(c + 2*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c**3 - d**3)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_174():
    f = (c + d*sec(e + f*x))**2/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -3*(c - d)**2*tan(e + f*x)/(16*a**2*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - (c - d)**2*tan(e + f*x)/(4*a**2*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)**2) - (c**2 - d**2)*tan(e + f*x)/(2*a**2*f*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) + 2*c**2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*c**2*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 3*sqrt(2)*(c - d)**2*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c**2 - d**2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_175():
    f = (c + d*sec(e + f*x))/(a*sec(e + f*x) + a)**(sympy.S(5)/2)
    F = -(c - d)*tan(e + f*x)/(4*f*(a*sec(e + f*x) + a)**(sympy.S(5)/2)) - (11*c - 3*d)*tan(e + f*x)/(16*a*f*(a*sec(e + f*x) + a)**(sympy.S(3)/2)) + 2*c*atan(sqrt(a)*tan(e + f*x)/sqrt(a*sec(e + f*x) + a))/(a**(sympy.S(5)/2)*f) - sqrt(2)*(43*c - 3*d)*atan(sqrt(2)*sqrt(a)*tan(e + f*x)/(2*sqrt(a*sec(e + f*x) + a)))/(32*a**(sympy.S(5)/2)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_176():
    f = 1/((c + d*sec(e + f*x))*(a*sec(e + f*x) + a)**(sympy.S(5)/2))
    F = -(c - 2*d)*tan(e + f*x)/(2*a**2*f*(c - d)**2*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - 3*tan(e + f*x)/(16*a**2*f*(c - d)*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - tan(e + f*x)/(4*a**2*f*(c - d)*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)**2) - sqrt(2)*(c - 2*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**2*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 3*sqrt(2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*f*(c - d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c**2 - 3*c*d + 3*d**2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(a**(sympy.S(3)/2)*f*(c - d)**3*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*d**(sympy.S(7)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(a**(sympy.S(3)/2)*c*f*(c - d)**3*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*c*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_177():
    f = 1/((c + d*sec(e + f*x))**2*(a*sec(e + f*x) + a)**(sympy.S(5)/2))
    F = -(c - 3*d)*tan(e + f*x)/(2*a**2*f*(c - d)**3*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - 3*tan(e + f*x)/(16*a**2*f*(c - d)**2*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - tan(e + f*x)/(4*a**2*f*(c - d)**2*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)**2) + d**4*tan(e + f*x)/(a**2*c*f*(c - d)**3*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - 3*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**3*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 3*sqrt(2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*f*(c - d)**2*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c**2 - 4*c*d + 6*d**2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(a**(sympy.S(3)/2)*f*(c - d)**4*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + d**(sympy.S(7)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(a**(sympy.S(3)/2)*c*f*(c - d)**3*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + d**(sympy.S(7)/2)*(8*c - 2*d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(a**(sympy.S(3)/2)*c**2*f*(c - d)**4*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*c**2*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_178():
    f = 1/((c + d*sec(e + f*x))**3*(a*sec(e + f*x) + a)**(sympy.S(5)/2))
    F = -(c - 4*d)*tan(e + f*x)/(2*a**2*f*(c - d)**4*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - 3*tan(e + f*x)/(16*a**2*f*(c - d)**3*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)) - tan(e + f*x)/(4*a**2*f*(c - d)**3*sqrt(a*sec(e + f*x) + a)*(sec(e + f*x) + 1)**2) + d**4*tan(e + f*x)/(2*a**2*c*f*(c - d)**3*(c + d)*(c + d*sec(e + f*x))**2*sqrt(a*sec(e + f*x) + a)) + 3*d**4*tan(e + f*x)/(4*a**2*c*f*(c - d)**3*(c + d)**2*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)) + d**4*(4*c - d)*tan(e + f*x)/(a**2*c**2*f*(c - d)**4*(c + d)*(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c - 4*d)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(4*a**(sympy.S(3)/2)*f*(c - d)**4*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - 3*sqrt(2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*f*(c - d)**3*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) - sqrt(2)*(c**2 - 5*c*d + 10*d**2)*tan(e + f*x)*atanh(sqrt(2)*sqrt(-a*sec(e + f*x) + a)/(2*sqrt(a)))/(a**(sympy.S(3)/2)*f*(c - d)**5*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 3*d**(sympy.S(7)/2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(4*a**(sympy.S(3)/2)*c*f*(c - d)**3*(c + d)**(sympy.S(5)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + d**(sympy.S(7)/2)*(4*c - d)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(a**(sympy.S(3)/2)*c**2*f*(c - d)**4*(c + d)**(sympy.S(3)/2)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*d**(sympy.S(7)/2)*(10*c**2 - 5*c*d + d**2)*tan(e + f*x)*atanh(sqrt(d)*sqrt(-a*sec(e + f*x) + a)/(sqrt(a)*sqrt(c + d)))/(a**(sympy.S(3)/2)*c**3*f*(c - d)**5*sqrt(c + d)*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a)) + 2*tan(e + f*x)*atanh(sqrt(-a*sec(e + f*x) + a)/sqrt(a))/(a**(sympy.S(3)/2)*c**3*f*sqrt(-a*sec(e + f*x) + a)*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_179():
    f = sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)
    F = 2*sqrt(a)*sqrt(c)*atan(sqrt(a)*sqrt(c)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/f + 2*sqrt(a)*sqrt(d)*atanh(sqrt(a)*sqrt(d)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_180():
    f = sqrt(a*sec(e + f*x) + a)/sqrt(c + d*sec(e + f*x))
    F = 2*sqrt(a)*atan(sqrt(a)*sqrt(c)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_181():
    f = sqrt(a*sec(e + f*x) + a)/(c + d*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*sqrt(a)*atan(sqrt(a)*sqrt(c)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(c**(sympy.S(3)/2)*f) - 2*a*d*tan(e + f*x)/(c*f*(c + d)*sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_182():
    f = sqrt(c + d*sec(e + f*x))/sqrt(a*sec(e + f*x) + a)
    F = 2*sqrt(c)*atan(sqrt(a)*sqrt(c)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f) - sqrt(2)*sqrt(c - d)*atan(sqrt(2)*sqrt(a)*sqrt(c - d)*tan(e + f*x)/(2*sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_183():
    f = 1/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a))
    F = -sqrt(2)*atan(sqrt(2)*sqrt(a)*sqrt(c - d)*tan(e + f*x)/(2*sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*f*sqrt(c - d)) + 2*atan(sqrt(a)*sqrt(c)*tan(e + f*x)/(sqrt(c + d*sec(e + f*x))*sqrt(a*sec(e + f*x) + a)))/(sqrt(a)*sqrt(c)*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_184():
    f = (a + b*sec(e + f*x))/(c + d*sec(e + f*x))
    F = a*x/c + (-2*a*d + 2*b*c)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c*f*sqrt(c - d)*sqrt(c + d))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_185():
    f = (a + b*sec(e + f*x))/(c + d*sec(e + f*x))**2
    F = a*x/c**2 - d*(-a*d + b*c)*tan(e + f*x)/(c*f*(c + d*sec(e + f*x))*(c**2 - d**2)) + (-4*a*c**2*d + 2*a*d**3 + 2*b*c**3)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c**2*f*(c - d)**(sympy.S(3)/2)*(c + d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_186():
    f = (a + b*sec(e + f*x))/(c + d*sec(e + f*x))**3
    F = a*x/c**3 - d*(-a*d + b*c)*tan(e + f*x)/(2*c*f*(c + d*sec(e + f*x))**2*(c**2 - d**2)) - d*(-5*a*c**2*d + 2*a*d**3 + 3*b*c**3)*tan(e + f*x)/(2*c**2*f*(c + d*sec(e + f*x))*(c**2 - d**2)**2) + (-a*d*(6*c**4 - 5*c**2*d**2 + 2*d**4) + b*c**3*(2*c**2 + d**2))*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c**3*f*(c - d)**(sympy.S(5)/2)*(c + d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_187():
    f = (a + b*sec(e + f*x))**2/(c + d*sec(e + f*x))**2
    F = a**2*x/c**2 + (-a*d + b*c)**2*sin(e + f*x)/(c*f*(c**2 - d**2)*(c*cos(e + f*x) + d)) + (-2*a*d + 2*b*c)*(2*a*c**2 - a*d**2 - b*c*d)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c**2*f*(c - d)**(sympy.S(3)/2)*(c + d)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_188():
    f = (a + b*sec(e + f*x))**2/(c + d*sec(e + f*x))**3
    F = a**2*x/c**3 - d*(-a*d + b*c)**2*sin(e + f*x)/(2*c**2*f*(c**2 - d**2)*(c*cos(e + f*x) + d)**2) - (-a*d + b*c)*(3*a*d*(2*c**2 - d**2) - b*c*(2*c**2 + d**2))*sin(e + f*x)/(2*c**2*f*(c**2 - d**2)**2*(c*cos(e + f*x) + d)) - (a**2*(6*c**4*d - 5*c**2*d**3 + 2*d**5) - 2*a*b*c**3*(2*c**2 + d**2) + 3*b**2*c**4*d)*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c**3*f*(c - d)**(sympy.S(5)/2)*(c + d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_189():
    f = (a + b*sec(e + f*x))**2/(c + d*sec(e + f*x))**4
    F = a**2*x/c**4 + d**2*(a*cos(e + f*x) + b)**2*sin(e + f*x)/(3*c*f*(c**2 - d**2)*(c*cos(e + f*x) + d)**3) - d*(-a*d + b*c)*(-8*a*c**2*d + 3*a*d**3 + 6*b*c**3 - b*c*d**2)*sin(e + f*x)/(6*c**3*f*(c**2 - d**2)**2*(c*cos(e + f*x) + d)**2) - (-a**2*d**2*(34*c**4 - 28*c**2*d**2 + 9*d**4) + 2*a*b*c*d*(18*c**4 - 5*c**2*d**2 + 2*d**4) - b**2*(6*c**6 + 10*c**4*d**2 - c**2*d**4))*sin(e + f*x)/(6*c**3*f*(c**2 - d**2)**3*(c*cos(e + f*x) + d)) - (a**2*(8*c**6*d - 8*c**4*d**3 + 7*c**2*d**5 - 2*d**7) - a*b*(4*c**7 + 6*c**5*d**2) + b**2*c**4*d*(4*c**2 + d**2))*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c**4*f*(c - d)**(sympy.S(7)/2)*(c + d)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_190():
    f = (a + b*sec(e + f*x))**3/(c + d*sec(e + f*x))**3
    F = a**3*x/c**3 + (-a*d + b*c)**2*(a*cos(e + f*x) + b)*sin(e + f*x)/(2*c*f*(c**2 - d**2)*(c*cos(e + f*x) + d)**2) + (-a*d + b*c)**2*(5*a*c**2 - 2*a*d**2 - 3*b*c*d)*sin(e + f*x)/(2*c**2*f*(c**2 - d**2)**2*(c*cos(e + f*x) + d)) - (-a*d + b*c)*(-a**2*(6*c**4 - 5*c**2*d**2 + 2*d**4) + 2*a*b*c*d*(4*c**2 - d**2) - b**2*c**2*(c**2 + 2*d**2))*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c**3*f*(c - d)**(sympy.S(5)/2)*(c + d)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_191():
    f = (a + b*sec(e + f*x))**3/(c + d*sec(e + f*x))**4
    F = a**3*x/c**4 - d*(-a*d + b*c)*(a*cos(e + f*x) + b)**2*sin(e + f*x)/(3*c*f*(c**2 - d**2)*(c*cos(e + f*x) + d)**3) + (-a*d + b*c)**2*(-8*a*c**2*d + 3*a*d**3 + 3*b*c**3 + 2*b*c*d**2)*sin(e + f*x)/(6*c**3*f*(c**2 - d**2)**2*(c*cos(e + f*x) + d)**2) - (-a*d + b*c)*(a**2*(34*c**4*d - 28*c**2*d**3 + 9*d**5) - a*b*c*(18*c**4 + 17*c**2*d**2 - 5*d**4) + b**2*c**2*d*(13*c**2 + 2*d**2))*sin(e + f*x)/(6*c**3*f*(c**2 - d**2)**3*(c*cos(e + f*x) + d)) - (a**3*(8*c**6*d - 8*c**4*d**3 + 7*c**2*d**5 - 2*d**7) - a**2*b*(6*c**7 + 9*c**5*d**2) + 3*a*b**2*c**4*d*(4*c**2 + d**2) - b**3*c**5*(c**2 + 4*d**2))*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(c**4*f*sqrt(c - d)*sqrt(c + d)*(c**2 - d**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_192():
    f = (a + b*sec(e + f*x))**3/(c + d*sec(e + f*x))**5
    F = a**3*x/c**5 + d**2*(a*cos(e + f*x) + b)**3*sin(e + f*x)/(4*c*f*(c**2 - d**2)*(c*cos(e + f*x) + d)**4) - d*(a*cos(e + f*x) + b)**2*(-11*a*c**2*d + 4*a*d**3 + 8*b*c**3 - b*c*d**2)*sin(e + f*x)/(12*c**2*f*(c**2 - d**2)**2*(c*cos(e + f*x) + d)**3) - (-a*d + b*c)*(-a**2*d**2*(58*c**4 - 35*c**2*d**2 + 12*d**4) + 2*a*b*c*d*(32*c**4 + c**2*d**2 + 2*d**4) - b**2*(12*c**6 + 25*c**4*d**2 - 2*c**2*d**4))*sin(e + f*x)/(24*c**4*f*(c**2 - d**2)**3*(c*cos(e + f*x) + d)**2) - (-a**3*(212*c**6*d**2 - 210*c**4*d**4 + 139*c**2*d**6 - 36*d**8) + a**2*b*c*d*(272*c**6 + 10*c**4*d**2 + 49*c**2*d**4 - 16*d**6) - 3*a*b**2*c**2*(24*c**6 + 84*c**4*d**2 - 5*c**2*d**4 + 2*d**6) + b**3*c**3*d*(68*c**4 + 39*c**2*d**2 - 2*d**4))*sin(e + f*x)/(24*c**4*f*(c**2 - d**2)**4*(c*cos(e + f*x) + d)) - (a**3*(40*c**8*d - 40*c**6*d**3 + 63*c**4*d**5 - 36*c**2*d**7 + 8*d**9) - 3*a**2*b*c**5*(8*c**4 + 24*c**2*d**2 + 3*d**4) + 15*a*b**2*c**6*d*(4*c**2 + 3*d**2) - b**3*c**5*(4*c**4 + 27*c**2*d**2 + 4*d**4))*atanh(sqrt(c - d)*tan(e/2 + f*x/2)/sqrt(c + d))/(4*c**5*f*sqrt(c - d)*sqrt(c + d)*(c**2 - d**2)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_193():
    f = sqrt(a + b*sec(e + f*x))*(c + d*sec(e + f*x))
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * (Symbol('c') + (Integer(-1) * Symbol('d')))) + (Symbol('a') * Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_194():
    f = sqrt(a + b*sec(e + f*x))/(c + d*sec(e + f*x))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('c') * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_195():
    f = (a + b*sec(e + f*x))**(sympy.S(3)/2)*(c + d*sec(e + f*x))
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('b') * Symbol('c')) + (Integer(4) * Symbol('a') * Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a') * Symbol('b') * ((Integer(6) * Symbol('c')) + (Integer(-1) * (Integer(4) * Symbol('d'))))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(3) * Symbol('c')) + (Integer(-1) * Symbol('d'))))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * Symbol('b') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_196():
    f = (a + b*sec(e + f*x))**(sympy.S(3)/2)/(c + d*sec(e + f*x))
    F = ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('d') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('c') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * Symbol('d') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_197():
    f = (a + b*sec(e + f*x))**(sympy.S(5)/2)*(c + d*sec(e + f*x))
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(35) * Symbol('a') * Symbol('b') * Symbol('c')) + (Integer(23) * (Symbol('a'))**(Integer(2)) * Symbol('d')) + (Integer(9) * (Symbol('b'))**(Integer(2)) * Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * Symbol('b') * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(45) * Symbol('c')) + (Integer(-1) * (Integer(23) * Symbol('d'))))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(35) * Symbol('c')) + (Integer(-1) * (Integer(17) * Symbol('d')))))) + ((Symbol('b'))**(Integer(3)) * ((Integer(5) * Symbol('c')) + (Integer(-1) * (Integer(9) * Symbol('d'))))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(15) * Symbol('b') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (Symbol('f'))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Integer(5) * Symbol('b') * Symbol('c')) + (Integer(8) * Symbol('a') * Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(15) * Symbol('f')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(5) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_198():
    f = (c + d*sec(e + f*x))/sqrt(a + b*sec(e + f*x))
    F = ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('b') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_199():
    f = 1/(sqrt(a + b*sec(e + f*x))*(c + d*sec(e + f*x)))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('c') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('d')) * ((Symbol('c') + Symbol('d')))**(Integer(-1))), sympy.asin((sympy.sqrt((Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt(Integer(2)))**(Integer(-1)))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * (Symbol('c') + Symbol('d')) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_200():
    f = (c + d*sec(e + f*x))/(a + b*sec(e + f*x))**(sympy.S(3)/2)
    F = ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_201():
    f = (c + d*sec(e + f*x))/(a + b*sec(e + f*x))**(sympy.S(5)/2)
    F = ((Integer(2) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('c'))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('d')))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('c'))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('c'))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('d'))) + ((Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('a'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + Symbol('b'))))**(Integer(-1)))), ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))) * sympy.sqrt(((Symbol('b') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((Symbol('b') * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) * (((Symbol('a'))**(Integer(3)) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * Symbol('b') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('c'))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('d')))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_202():
    f = sqrt(a + b*sec(e + f*x))*sqrt(c + d*sec(e + f*x))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('f')))**(Integer(-1)))) + ((Integer(2) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(((sympy.sqrt(((Symbol('a') + Symbol('b')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt(((Symbol('a') + Symbol('b')) * ((Symbol('c') + Symbol('d')))**(Integer(-1)))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_203():
    f = sqrt(a + b*sec(e + f*x))/sqrt(c + d*sec(e + f*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_204():
    f = sqrt(a + b*sec(e + f*x))/(c + d*sec(e + f*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('c'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('c') * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_205():
    f = sqrt(a + b*sec(e + f*x))/(c + d*sec(e + f*x))**(sympy.S(5)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * ((Integer(6) * Symbol('b') * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(3)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('c'))**(Integer(2)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(3) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(2)))))) + (Integer(-1) * (Symbol('a') * Symbol('d') * ((Integer(9) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(6) * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(3) * (Symbol('d'))**(Integer(3))))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('c'))**(Integer(3)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_206():
    f = (a + b*sec(e + f*x))**(sympy.S(3)/2)/(c + d*sec(e + f*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('c') * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * ((Integer(2) * Symbol('c')) + (Integer(-1) * Symbol('d')))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(2)) * (Symbol('c') + (Integer(-1) * Symbol('d'))) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(2)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_207():
    f = (a + b*sec(e + f*x))**(sympy.S(3)/2)/(c + d*sec(e + f*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('b') * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(3)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('c'))**(Integer(2)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * ((Integer(3) * Symbol('c')) + Symbol('d'))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Integer(3) * (Symbol('c'))**(Integer(2))) + (Integer(2) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))) + ((Symbol('a'))**(Integer(2)) * Symbol('d') * ((Integer(9) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(6) * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(3) * (Symbol('d'))**(Integer(3)))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('c'))**(Integer(3)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_208():
    f = (a + b*sec(e + f*x))**(sympy.S(3)/2)/(c + d*sec(e + f*x))**(sympy.S(7)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d') * ((Integer(35) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(8) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(5) * (Symbol('d'))**(Integer(4))))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(58) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(41) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(15) * (Symbol('d'))**(Integer(4)))))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(15) * (Symbol('c'))**(Integer(6))) + (Integer(19) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4)))))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(15) * (Symbol('c'))**(Integer(3)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(3)) * ((Symbol('c') + Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * ((Integer(15) * (Symbol('c'))**(Integer(3))) + (Integer(10) * (Symbol('c'))**(Integer(2)) * Symbol('d')) + (Integer(9) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(3)))))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Integer(15) * (Symbol('c'))**(Integer(4))) + (Integer(20) * (Symbol('c'))**(Integer(3)) * Symbol('d')) + (Integer(-1) * (Integer(4) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(-1) * (Integer(4) * Symbol('c') * (Symbol('d'))**(Integer(3)))) + (Integer(5) * (Symbol('d'))**(Integer(4)))))) + ((Symbol('a'))**(Integer(2)) * Symbol('d') * ((Integer(60) * (Symbol('c'))**(Integer(5))) + (Integer(-1) * (Integer(2) * (Symbol('c'))**(Integer(4)) * Symbol('d'))) + (Integer(-1) * (Integer(66) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)))) + (Integer(25) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(3))) + (Integer(30) * Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(15) * (Symbol('d'))**(Integer(5))))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(15) * (Symbol('c'))**(Integer(4)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(3)) * ((Symbol('c') + Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(4)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(5) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Integer(10) * Symbol('b') * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(13) * Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(5) * Symbol('a') * (Symbol('d'))**(Integer(3)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(15) * (Symbol('c'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_209():
    f = (a + b*sec(e + f*x))**(sympy.S(5)/2)/(c + d*sec(e + f*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(7) * Symbol('a') * (Symbol('c'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('b') * Symbol('c') * Symbol('d'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(2))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('c'))**(Integer(2)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * (Symbol('c') + (Integer(3) * Symbol('d')))) + (Integer(-1) * (Symbol('a') * Symbol('b') * Symbol('c') * ((Integer(7) * (Symbol('c'))**(Integer(2))) + (Integer(4) * Symbol('c') * Symbol('d')) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(2))))))) + ((Symbol('a'))**(Integer(2)) * ((Integer(9) * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(2) * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(6) * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(3) * (Symbol('d'))**(Integer(3)))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(3) * (Symbol('c'))**(Integer(3)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(2)) * ((Symbol('c') + Symbol('d')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(3)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(3) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_210():
    f = (a + b*sec(e + f*x))**(sympy.S(5)/2)/(c + d*sec(e + f*x))**(sympy.S(7)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * Symbol('d') * ((Integer(29) * (Symbol('c'))**(Integer(2))) + (Integer(3) * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * (Symbol('a') * Symbol('b') * Symbol('c') * ((Integer(35) * (Symbol('c'))**(Integer(4))) + (Integer(34) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('d'))**(Integer(4))))))) + ((Symbol('a'))**(Integer(2)) * ((Integer(58) * (Symbol('c'))**(Integer(4)) * Symbol('d')) + (Integer(-1) * (Integer(41) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(3)))) + (Integer(15) * (Symbol('d'))**(Integer(5)))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(15) * (Symbol('c'))**(Integer(3)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(3)) * ((Symbol('c') + Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(4)) * ((Integer(5) * (Symbol('c'))**(Integer(2))) + (Integer(24) * Symbol('c') * Symbol('d')) + (Integer(3) * (Symbol('d'))**(Integer(2))))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * ((Integer(35) * (Symbol('c'))**(Integer(3))) + (Integer(42) * (Symbol('c'))**(Integer(2)) * Symbol('d')) + (Integer(21) * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('d'))**(Integer(3))))))) + ((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Integer(45) * (Symbol('c'))**(Integer(4))) + (Integer(48) * (Symbol('c'))**(Integer(3)) * Symbol('d')) + ((Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(8) * Symbol('c') * (Symbol('d'))**(Integer(3)))) + (Integer(10) * (Symbol('d'))**(Integer(4))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('d') * ((Integer(60) * (Symbol('c'))**(Integer(5))) + (Integer(-1) * (Integer(2) * (Symbol('c'))**(Integer(4)) * Symbol('d'))) + (Integer(-1) * (Integer(66) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2)))) + (Integer(25) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(3))) + (Integer(30) * Symbol('c') * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(15) * (Symbol('d'))**(Integer(5)))))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(15) * (Symbol('c'))**(Integer(4)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(3)) * ((Symbol('c') + Symbol('d')))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(4)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(5) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * ((Integer(5) * Symbol('b') * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(13) * Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(3) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2))) + (Integer(5) * Symbol('a') * (Symbol('d'))**(Integer(3)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(15) * (Symbol('c'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_211():
    f = (a + b*sec(e + f*x))**(sympy.S(5)/2)/(c + d*sec(e + f*x))**(sympy.S(9)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(3)) * Symbol('d') * ((Integer(133) * (Symbol('c'))**(Integer(4))) + (Integer(62) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('d'))**(Integer(4)))))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('c') * Symbol('d') * ((Integer(406) * (Symbol('c'))**(Integer(6))) + (Integer(73) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(132) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(35) * (Symbol('d'))**(Integer(6)))))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(2)) * ((Integer(245) * (Symbol('c'))**(Integer(6))) + (Integer(852) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(41) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(14) * (Symbol('d'))**(Integer(6)))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * ((Integer(582) * (Symbol('c'))**(Integer(6)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(485) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(4)))) + (Integer(392) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(6))) + (Integer(-1) * (Integer(105) * (Symbol('d'))**(Integer(8)))))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_e(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(105) * (Symbol('c'))**(Integer(4)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(4)) * ((Symbol('c') + Symbol('d')))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * (((Symbol('b'))**(Integer(3)) * (Symbol('c'))**(Integer(4)) * ((Integer(35) * (Symbol('c'))**(Integer(4))) + (Integer(231) * (Symbol('c'))**(Integer(3)) * Symbol('d')) + (Integer(67) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2))) + (Integer(57) * Symbol('c') * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(6) * (Symbol('d'))**(Integer(4)))))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('c'))**(Integer(3)) * ((Integer(245) * (Symbol('c'))**(Integer(5))) + (Integer(413) * (Symbol('c'))**(Integer(4)) * Symbol('d')) + (Integer(439) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(2))) + (Integer(53) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(3))) + (Integer(-1) * (Integer(12) * Symbol('c') * (Symbol('d'))**(Integer(4)))) + (Integer(14) * (Symbol('d'))**(Integer(5)))))) + ((Symbol('a'))**(Integer(2)) * Symbol('b') * (Symbol('c'))**(Integer(2)) * ((Integer(315) * (Symbol('c'))**(Integer(6))) + (Integer(497) * (Symbol('c'))**(Integer(5)) * Symbol('d')) + (Integer(219) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(73) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(3)))) + (Integer(208) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4))) + (Integer(56) * Symbol('c') * (Symbol('d'))**(Integer(5))) + (Integer(-1) * (Integer(70) * (Symbol('d'))**(Integer(6)))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('d') * ((Integer(525) * (Symbol('c'))**(Integer(7))) + (Integer(57) * (Symbol('c'))**(Integer(6)) * Symbol('d')) + (Integer(-1) * (Integer(699) * (Symbol('c'))**(Integer(5)) * (Symbol('d'))**(Integer(2)))) + (Integer(214) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(3))) + (Integer(672) * (Symbol('c'))**(Integer(3)) * (Symbol('d'))**(Integer(4))) + (Integer(-1) * (Integer(280) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(5)))) + (Integer(-1) * (Integer(210) * Symbol('c') * (Symbol('d'))**(Integer(6)))) + (Integer(105) * (Symbol('d'))**(Integer(7))))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(105) * (Symbol('c'))**(Integer(5)) * ((Symbol('c') + (Integer(-1) * Symbol('d'))))**(Integer(4)) * ((Symbol('c') + Symbol('d')))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.cos((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.csc((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * (((Symbol('c'))**(Integer(5)) * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f') * sympy.sqrt((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(7) * Symbol('c') * ((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))) * Symbol('f') * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Integer(3)) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('d') * ((Integer(14) * Symbol('b') * (Symbol('c'))**(Integer(3))) + (Integer(-1) * (Integer(19) * Symbol('a') * (Symbol('c'))**(Integer(2)) * Symbol('d'))) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('c') * (Symbol('d'))**(Integer(2)))) + (Integer(7) * Symbol('a') * (Symbol('d'))**(Integer(3)))) * (Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(35) * (Symbol('c'))**(Integer(2)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(2)) * Symbol('f') * ((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Integer(2)) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Integer(2) * Symbol('a') * Symbol('b') * Symbol('c') * Symbol('d') * ((Integer(91) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(2) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(7) * (Symbol('d'))**(Integer(4))))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * ((Integer(162) * (Symbol('c'))**(Integer(4))) + (Integer(-1) * (Integer(101) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(2)))) + (Integer(35) * (Symbol('d'))**(Integer(4)))))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * ((Integer(35) * (Symbol('c'))**(Integer(6))) + (Integer(67) * (Symbol('c'))**(Integer(4)) * (Symbol('d'))**(Integer(2))) + (Integer(-1) * (Integer(6) * (Symbol('c'))**(Integer(2)) * (Symbol('d'))**(Integer(4)))))))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(105) * (Symbol('c'))**(Integer(3)) * (((Symbol('c'))**(Integer(2)) + (Integer(-1) * (Symbol('d'))**(Integer(2)))))**(Integer(3)) * Symbol('f') * (Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_212():
    f = (c + d*sec(e + f*x))**(sympy.S(3)/2)/sqrt(a + b*sec(e + f*x))
    F = (Integer(-1) * ((Integer(2) * Symbol('c') * (Symbol('c') + Symbol('d')) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(sympy.sqrt((((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((((Symbol('a') + Symbol('b')) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(-1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('c') + Symbol('d')))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))) + ((Integer(2) * Symbol('d') * (Symbol('c') + Symbol('d')) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('b') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))), sympy.asin(sympy.sqrt((((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (((Symbol('a') + Symbol('b')) * (((Integer(-1) * Symbol('b')) * Symbol('c')) + (Symbol('a') * Symbol('d'))) * (Integer(-1) + sympy.sec((Symbol('e') + (Symbol('f') * x)))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('c') + Symbol('d')))**(Integer(2)) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)))))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('f') * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(sympy.sqrt((((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(-1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((Symbol('a') * Symbol('b') * Symbol('f') * sympy.sqrt((((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_213():
    f = sqrt(c + d*sec(e + f*x))/sqrt(a + b*sec(e + f*x))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((((Symbol('a') + Symbol('b')) * Symbol('c')) * ((Symbol('a') * (Symbol('c') + Symbol('d'))))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('c') + Symbol('d'))) * Symbol('f')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_214():
    f = 1/(sqrt(a + b*sec(e + f*x))*sqrt(c + d*sec(e + f*x)))
    F = (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')(((Symbol('a') * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * Symbol('c')))**(Integer(-1))), sympy.asin(((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') + Symbol('d')) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('c') + (Integer(-1) * Symbol('d'))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('c') * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.elliptic_f(sympy.asin(((sympy.sqrt((Symbol('c') + Symbol('d'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1)))), (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Integer(-1) * Symbol('d')))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + Symbol('d'))))**(Integer(-1)))) * sympy.sqrt(((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('a') + Symbol('b')) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * ((((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * (Integer(1) + sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))) * (Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('a') * sympy.sqrt((Symbol('c') + Symbol('d'))) * ((Symbol('b') * Symbol('c')) + (Integer(-1) * (Symbol('a') * Symbol('d')))) * Symbol('f')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_215():
    f = (a + b*sec(e + f*x))**(sympy.S(1)/3)/(c + d*sec(e + f*x))**(sympy.S(1)/3)
    F = (((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1))) * sympy.Function('Unintegrable')((((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1))) * (((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1))))**(Integer(-1))), x)) * ((((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1))) * ((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1)))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_216():
    f = (a + b*sec(e + f*x))**(sympy.S(1)/3)/(c + d*sec(e + f*x))**(sympy.S(4)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1))) * (((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_217():
    f = (a + b*sec(e + f*x))**(sympy.S(1)/3)/(c + d*sec(e + f*x))**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(3))**(Integer(-1))) * (((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(7) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_218():
    f = (a + b*sec(e + f*x))**(sympy.S(2)/3)/(c + d*sec(e + f*x))**(sympy.S(2)/3)
    F = (((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * sympy.Function('Unintegrable')((((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)) * ((((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_219():
    f = (a + b*sec(e + f*x))**(sympy.S(2)/3)/(c + d*sec(e + f*x))**(sympy.S(5)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(5) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_220():
    f = (a + b*sec(e + f*x))**(sympy.S(2)/3)/(c + d*sec(e + f*x))**(sympy.S(8)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(2) * (Integer(3))**(Integer(-1)))) * (((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(8) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_221():
    f = (a + b*sec(e + f*x))**(sympy.S(4)/3)/(c + d*sec(e + f*x))**(sympy.S(4)/3)
    F = (((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * sympy.Function('Unintegrable')((((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * (((Symbol('d') + (Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)) * ((((Symbol('b') + (Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * ((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_222():
    f = (a + b*sec(e + f*x))**(sympy.S(4)/3)/(c + d*sec(e + f*x))**(sympy.S(7)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * (((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(7) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_223():
    f = (a + b*sec(e + f*x))**(sympy.S(4)/3)/(c + d*sec(e + f*x))**(sympy.S(10)/3)
    F = sympy.Function('Unintegrable')((((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(4) * (Integer(3))**(Integer(-1)))) * (((Symbol('c') + (Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**((Integer(10) * (Integer(3))**(Integer(-1)))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_224():
    f = (c*(d*sec(e + f*x))**p)**n*(a*sec(e + f*x) + a)**m
    F = -(c*(d*sec(e + f*x))**p)**n*(a*sec(e + f*x) + a)**m*(sec(e + f*x) + 1)**(-m + sympy.S(-1)/2)*tan(e + f*x)*appellf1(n*p, sympy.S.Half, sympy.S.Half - m, n*p + 1, sec(e + f*x), -sec(e + f*x))/(f*n*p*sqrt(1 - sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_225():
    f = (c*(d*sec(e + f*x))**p)**n*(a*sec(e + f*x) + a)**3
    F = -a**3*(c*(d*sec(e + f*x))**p)**n*(4*n*p + 1)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(-n**2*p**2 + 1)*sqrt(sin(e + f*x)**2)) + a**3*(c*(d*sec(e + f*x))**p)**n*(2*n*p + 5)*tan(e + f*x)/(f*(n*p + 1)*(n*p + 2)) + a**3*(c*(d*sec(e + f*x))**p)**n*(4*n*p + 7)*sin(e + f*x)*hyper((sympy.S.Half, -n*p/2), (-n*p/2 + 1,), cos(e + f*x)**2)/(f*n*p*(n*p + 2)*sqrt(sin(e + f*x)**2)) + (c*(d*sec(e + f*x))**p)**n*(a**3*sec(e + f*x) + a**3)*tan(e + f*x)/(f*(n*p + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_226():
    f = (c*(d*sec(e + f*x))**p)**n*(a*sec(e + f*x) + a)**2
    F = -a**2*(c*(d*sec(e + f*x))**p)**n*(2*n*p + 1)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(-n**2*p**2 + 1)*sqrt(sin(e + f*x)**2)) + a**2*(c*(d*sec(e + f*x))**p)**n*tan(e + f*x)/(f*(n*p + 1)) + 2*a**2*(c*(d*sec(e + f*x))**p)**n*sin(e + f*x)*hyper((sympy.S.Half, -n*p/2), (-n*p/2 + 1,), cos(e + f*x)**2)/(f*n*p*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_227():
    f = (c*(d*sec(e + f*x))**p)**n*(a*sec(e + f*x) + a)
    F = -a*(c*(d*sec(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(-n*p + 1)*sqrt(sin(e + f*x)**2)) + a*(c*(d*sec(e + f*x))**p)**n*sin(e + f*x)*hyper((sympy.S.Half, -n*p/2), (-n*p/2 + 1,), cos(e + f*x)**2)/(f*n*p*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_228():
    f = (c*(d*sec(e + f*x))**p)**n/(a*sec(e + f*x) + a)
    F = (c*(d*sec(e + f*x))**p)**n*sin(e + f*x)/(f*(a*sec(e + f*x) + a)) + (c*(d*sec(e + f*x))**p)**n*(-n*p + 1)*sin(e + f*x)*cos(e + f*x)**2*hyper((sympy.S.Half, -n*p/2 + 1), (-n*p/2 + 2,), cos(e + f*x)**2)/(a*f*(-n*p + 2)*sqrt(sin(e + f*x)**2)) - (c*(d*sec(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(a*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_229():
    f = (c*(d*sec(e + f*x))**p)**n/(a*sec(e + f*x) + a)**2
    F = -(c*(d*sec(e + f*x))**p)**n*tan(e + f*x)/(3*f*(a*sec(e + f*x) + a)**2) - (c*(d*sec(e + f*x))**p)**n*(-2*n*p + 3)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(3*a**2*f*sqrt(sin(e + f*x)**2)) + (c*(d*sec(e + f*x))**p)**n*(-2*n*p + 4)*sin(e + f*x)*hyper((sympy.S.Half, -n*p/2), (-n*p/2 + 1,), cos(e + f*x)**2)/(3*a**2*f*sqrt(sin(e + f*x)**2)) - (c*(d*sec(e + f*x))**p)**n*(-2*n*p + 4)*tan(e + f*x)/(3*a**2*f*(sec(e + f*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_230():
    f = (c*(d*sec(e + f*x))**p)**n*(a + b*sec(e + f*x))**m
    F = (((Symbol('c') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('p'))))**(Symbol('n')) * sympy.Function('Unintegrable')((((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Symbol('n') * Symbol('p'))) * ((Symbol('a') + (Symbol('b') * sympy.sec((Symbol('e') + (Symbol('f') * x))))))**(Symbol('m'))), x)) * (((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Symbol('n') * Symbol('p'))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_231():
    f = (c*(d*sec(e + f*x))**p)**n*(a + b*sec(e + f*x))**3
    F = a*b**2*(c*(d*sec(e + f*x))**p)**n*(2*n*p + 5)*tan(e + f*x)/(f*(n*p + 1)*(n*p + 2)) - a*(c*(d*sec(e + f*x))**p)**n*(a**2*(n*p + 1) + 3*b**2*n*p)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(-n**2*p**2 + 1)*sqrt(sin(e + f*x)**2)) + b**2*(c*(d*sec(e + f*x))**p)**n*(a + b*sec(e + f*x))*tan(e + f*x)/(f*(n*p + 2)) + b*(c*(d*sec(e + f*x))**p)**n*(3*a**2*(n*p + 2) + b**2*(n*p + 1))*sin(e + f*x)*hyper((sympy.S.Half, -n*p/2), (-n*p/2 + 1,), cos(e + f*x)**2)/(f*n*p*(n*p + 2)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_232():
    f = (c*(d*sec(e + f*x))**p)**n*(a + b*sec(e + f*x))**2
    F = 2*a*b*(c*(d*sec(e + f*x))**p)**n*sin(e + f*x)*hyper((sympy.S.Half, -n*p/2), (-n*p/2 + 1,), cos(e + f*x)**2)/(f*n*p*sqrt(sin(e + f*x)**2)) + b**2*(c*(d*sec(e + f*x))**p)**n*tan(e + f*x)/(f*(n*p + 1)) - (c*(d*sec(e + f*x))**p)**n*(a**2*(n*p + 1) + b**2*n*p)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(-n**2*p**2 + 1)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_233():
    f = (c*(d*sec(e + f*x))**p)**n*(a + b*sec(e + f*x))
    F = -a*(c*(d*sec(e + f*x))**p)**n*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, -n*p/2 + sympy.S.Half), (-n*p/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(-n*p + 1)*sqrt(sin(e + f*x)**2)) + b*(c*(d*sec(e + f*x))**p)**n*sin(e + f*x)*hyper((sympy.S.Half, -n*p/2), (-n*p/2 + 1,), cos(e + f*x)**2)/(f*n*p*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_234():
    f = (c*(d*sec(e + f*x))**p)**n/(a + b*sec(e + f*x))
    F = a*(c*(d*sec(e + f*x))**p)**n*(cos(e + f*x)**2)**(n*p/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 1, n*p/2 + sympy.S(-1)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)) - b*(c*(d*sec(e + f*x))**p)**n*(cos(e + f*x)**2)**(n*p/2)*sin(e + f*x)*appellf1(sympy.S.Half, 1, n*p/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_5_Secant_4_5_2_1_a_plus_b_sec_pow_m_c_plus_d_sec_pow_n_235():
    f = (c*(d*sec(e + f*x))**p)**n/(a + b*sec(e + f*x))**2
    F = a**2*(c*(d*sec(e + f*x))**p)**n*(cos(e + f*x)**2)**(n*p/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, n*p/2 + sympy.S(-3)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2) - 2*a*b*(c*(d*sec(e + f*x))**p)**n*(cos(e + f*x)**2)**(n*p/2)*sin(e + f*x)*appellf1(sympy.S.Half, 2, n*p/2 - 1, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2) + b**2*(c*(d*sec(e + f*x))**p)**n*(cos(e + f*x)**2)**(n*p/2 + sympy.S(-1)/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 2, n*p/2 + sympy.S(-1)/2, sympy.S(3)/2, a**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)**2)
    assert integrate(f, x) == F

