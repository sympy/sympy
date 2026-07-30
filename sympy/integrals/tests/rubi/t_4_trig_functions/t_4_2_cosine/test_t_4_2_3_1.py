"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.2 Cosine/4.2.3.1 (a+b cos)^m (c+d cos)^n (A+B cos).m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

A, B, a, b, c, d, e, f, m, n = symbols('A B a b c d e f m n')

def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_1():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*cos(c + d*x)**3
    F = B*a*sin(c + d*x)*cos(c + d*x)**4/(5*d) + 3*a*x*(A + B)/8 + a*(A + B)*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*(A + B)*sin(c + d*x)*cos(c + d*x)/(8*d) - a*(5*A + 4*B)*sin(c + d*x)**3/(15*d) + a*(5*A + 4*B)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_2():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*cos(c + d*x)**2
    F = B*a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + a*x*(4*A + 3*B)/8 - a*(A + B)*sin(c + d*x)**3/(3*d) + a*(A + B)*sin(c + d*x)/d + a*(4*A + 3*B)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_3():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*cos(c + d*x)
    F = B*a*sin(c + d*x)*cos(c + d*x)**2/(3*d) + a*x*(A + B)/2 + a*(A + B)*sin(c + d*x)*cos(c + d*x)/(2*d) + a*(3*A + 2*B)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_4():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)
    F = B*a*sin(c + d*x)*cos(c + d*x)/(2*d) + a*x*(2*A + B)/2 + a*(A + B)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_5():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)
    F = A*a*atanh(sin(c + d*x))/d + B*a*sin(c + d*x)/d + a*x*(A + B)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_6():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)**2
    F = A*a*tan(c + d*x)/d + B*a*x + a*(A + B)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_7():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)**3
    F = A*a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*(A + B)*tan(c + d*x)/d + a*(A + 2*B)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_8():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)**4
    F = A*a*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a*(A + B)*tan(c + d*x)*sec(c + d*x)/(2*d) + a*(A + B)*atanh(sin(c + d*x))/(2*d) + a*(2*A + 3*B)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_9():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)**5
    F = A*a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a*(A + B)*tan(c + d*x)**3/(3*d) + a*(A + B)*tan(c + d*x)/d + a*(3*A + 4*B)*tan(c + d*x)*sec(c + d*x)/(8*d) + a*(3*A + 4*B)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_10():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*cos(c + d*x)**3
    F = B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**4/(6*d) + a**2*x*(12*A + 11*B)/16 + a**2*(6*A + 7*B)*sin(c + d*x)*cos(c + d*x)**4/(30*d) - a**2*(9*A + 8*B)*sin(c + d*x)**3/(15*d) + a**2*(9*A + 8*B)*sin(c + d*x)/(5*d) + a**2*(12*A + 11*B)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + a**2*(12*A + 11*B)*sin(c + d*x)*cos(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_11():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*cos(c + d*x)**2
    F = B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**3/(5*d) + a**2*x*(7*A + 6*B)/8 + a**2*(5*A + 6*B)*sin(c + d*x)*cos(c + d*x)**3/(20*d) + a**2*(7*A + 6*B)*sin(c + d*x)*cos(c + d*x)/(8*d) - a**2*(10*A + 9*B)*sin(c + d*x)**3/(15*d) + a**2*(10*A + 9*B)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_12():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*cos(c + d*x)
    F = B*(a*cos(c + d*x) + a)**3*sin(c + d*x)/(4*a*d) + a**2*x*(8*A + 7*B)/8 + a**2*(8*A + 7*B)*sin(c + d*x)*cos(c + d*x)/(24*d) + a**2*(8*A + 7*B)*sin(c + d*x)/(6*d) + (4*A - B)*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(12*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_13():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2
    F = B*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(3*d) + a**2*x*(3*A + 2*B)/2 + a**2*(3*A + 2*B)*sin(c + d*x)*cos(c + d*x)/(6*d) + 2*a**2*(3*A + 2*B)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_14():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)
    F = A*a**2*atanh(sin(c + d*x))/d + B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(2*d) + a**2*x*(4*A + 3*B)/2 + a**2*(2*A + 3*B)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_15():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)**2
    F = A*(a**2*cos(c + d*x) + a**2)*tan(c + d*x)/d + a**2*x*(A + 2*B) - a**2*(A - B)*sin(c + d*x)/d + a**2*(2*A + B)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_16():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)**3
    F = A*(a**2*cos(c + d*x) + a**2)*tan(c + d*x)*sec(c + d*x)/(2*d) + B*a**2*x + a**2*(3*A + 2*B)*tan(c + d*x)/(2*d) + a**2*(3*A + 4*B)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_17():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)**4
    F = A*(a**2*cos(c + d*x) + a**2)*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a**2*(2*A + 3*B)*atanh(sin(c + d*x))/(2*d) + a**2*(4*A + 3*B)*tan(c + d*x)*sec(c + d*x)/(6*d) + a**2*(5*A + 6*B)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_18():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)**5
    F = A*(a**2*cos(c + d*x) + a**2)*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a**2*(4*A + 5*B)*tan(c + d*x)/(3*d) + a**2*(5*A + 4*B)*tan(c + d*x)*sec(c + d*x)**2/(12*d) + a**2*(7*A + 8*B)*tan(c + d*x)*sec(c + d*x)/(8*d) + a**2*(7*A + 8*B)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_19():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*cos(c + d*x)**2
    F = B*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)**3/(6*d) + a**3*x*(26*A + 23*B)/16 - a**3*(19*A + 17*B)*sin(c + d*x)**3/(15*d) + a**3*(19*A + 17*B)*sin(c + d*x)/(5*d) + a**3*(22*A + 21*B)*sin(c + d*x)*cos(c + d*x)**3/(40*d) + a**3*(26*A + 23*B)*sin(c + d*x)*cos(c + d*x)/(16*d) + (3*A + 4*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)*cos(c + d*x)**3/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_20():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*cos(c + d*x)
    F = B*(a*cos(c + d*x) + a)**4*sin(c + d*x)/(5*a*d) + a**3*x*(15*A + 13*B)/8 - a**3*(15*A + 13*B)*sin(c + d*x)**3/(60*d) + 3*a**3*(15*A + 13*B)*sin(c + d*x)*cos(c + d*x)/(40*d) + a**3*(15*A + 13*B)*sin(c + d*x)/(5*d) + (5*A - B)*(a*cos(c + d*x) + a)**3*sin(c + d*x)/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_21():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3
    F = B*(a*cos(c + d*x) + a)**3*sin(c + d*x)/(4*d) + 5*a**3*x*(4*A + 3*B)/8 - a**3*(4*A + 3*B)*sin(c + d*x)**3/(12*d) + 3*a**3*(4*A + 3*B)*sin(c + d*x)*cos(c + d*x)/(8*d) + a**3*(4*A + 3*B)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_22():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)
    F = A*a**3*atanh(sin(c + d*x))/d + B*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(3*d) + a**3*x*(7*A + 5*B)/2 + 5*a**3*(A + B)*sin(c + d*x)/(2*d) + (3*A + 5*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_23():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**2
    F = A*a*(a*cos(c + d*x) + a)**2*tan(c + d*x)/d + 5*B*a**3*sin(c + d*x)/(2*d) + a**3*x*(6*A + 7*B)/2 + a**3*(3*A + B)*atanh(sin(c + d*x))/d - (2*A - B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_24():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**3
    F = -5*A*a**3*sin(c + d*x)/(2*d) + A*a*(a*cos(c + d*x) + a)**2*tan(c + d*x)*sec(c + d*x)/(2*d) + a**3*x*(A + 3*B) + a**3*(7*A + 6*B)*atanh(sin(c + d*x))/(2*d) + (2*A + B)*(a**3*cos(c + d*x) + a**3)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_25():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**4
    F = A*a*(a*cos(c + d*x) + a)**2*tan(c + d*x)*sec(c + d*x)**2/(3*d) + B*a**3*x + 5*a**3*(A + B)*tan(c + d*x)/(2*d) + a**3*(5*A + 7*B)*atanh(sin(c + d*x))/(2*d) + (5*A + 3*B)*(a**3*cos(c + d*x) + a**3)*tan(c + d*x)*sec(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_26():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**5
    F = A*a*(a*cos(c + d*x) + a)**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 5*a**3*(3*A + 4*B)*atanh(sin(c + d*x))/(8*d) + a**3*(9*A + 11*B)*tan(c + d*x)/(3*d) + a**3*(27*A + 28*B)*tan(c + d*x)*sec(c + d*x)/(24*d) + (3*A + 2*B)*(a**3*cos(c + d*x) + a**3)*tan(c + d*x)*sec(c + d*x)**2/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_27():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**6
    F = A*a*(a*cos(c + d*x) + a)**2*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a**3*(13*A + 15*B)*tan(c + d*x)*sec(c + d*x)/(8*d) + a**3*(13*A + 15*B)*atanh(sin(c + d*x))/(8*d) + a**3*(38*A + 45*B)*tan(c + d*x)/(15*d) + a**3*(43*A + 45*B)*tan(c + d*x)*sec(c + d*x)**2/(60*d) + (7*A + 5*B)*(a**3*cos(c + d*x) + a**3)*tan(c + d*x)*sec(c + d*x)**3/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_28():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*cos(c + d*x)**2
    F = B*a*(a*cos(c + d*x) + a)**3*sin(c + d*x)*cos(c + d*x)**3/(7*d) + a**4*x*(49*A + 44*B)/16 + a**4*(49*A + 44*B)*sin(c + d*x)*cos(c + d*x)/(16*d) - a**4*(252*A + 227*B)*sin(c + d*x)**3/(105*d) + a**4*(252*A + 227*B)*sin(c + d*x)/(35*d) + a**4*(301*A + 276*B)*sin(c + d*x)*cos(c + d*x)**3/(280*d) + (7*A + 7*B)*(a**4*cos(c + d*x) + a**4)*sin(c + d*x)*cos(c + d*x)**3/(15*d) + (7*A + 10*B)*(a**2*cos(c + d*x) + a**2)**2*sin(c + d*x)*cos(c + d*x)**3/(42*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_29():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*cos(c + d*x)
    F = B*(a*cos(c + d*x) + a)**5*sin(c + d*x)/(6*a*d) + 7*a**4*x*(8*A + 7*B)/16 - 2*a**4*(8*A + 7*B)*sin(c + d*x)**3/(15*d) + a**4*(8*A + 7*B)*sin(c + d*x)*cos(c + d*x)**3/(40*d) + 27*a**4*(8*A + 7*B)*sin(c + d*x)*cos(c + d*x)/(80*d) + 4*a**4*(8*A + 7*B)*sin(c + d*x)/(5*d) + (6*A - B)*(a*cos(c + d*x) + a)**4*sin(c + d*x)/(30*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_30():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4
    F = B*(a*cos(c + d*x) + a)**4*sin(c + d*x)/(5*d) + 7*a**4*x*(5*A + 4*B)/8 - 4*a**4*(5*A + 4*B)*sin(c + d*x)**3/(15*d) + a**4*(5*A + 4*B)*sin(c + d*x)*cos(c + d*x)**3/(20*d) + 27*a**4*(5*A + 4*B)*sin(c + d*x)*cos(c + d*x)/(40*d) + 8*a**4*(5*A + 4*B)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_31():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*sec(c + d*x)
    F = A*a**4*atanh(sin(c + d*x))/d + B*a*(a*cos(c + d*x) + a)**3*sin(c + d*x)/(4*d) + a**4*x*(48*A + 35*B)/8 + 5*a**4*(8*A + 7*B)*sin(c + d*x)/(8*d) + (4*A + 7*B)*(a**2*cos(c + d*x) + a**2)**2*sin(c + d*x)/(12*d) + (32*A + 35*B)*(a**4*cos(c + d*x) + a**4)*sin(c + d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_32():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*sec(c + d*x)**2
    F = A*a*(a*cos(c + d*x) + a)**3*tan(c + d*x)/d + a**4*x*(13*A + 12*B)/2 + 5*a**4*(A + 2*B)*sin(c + d*x)/(2*d) + a**4*(4*A + B)*atanh(sin(c + d*x))/d - (3*A - 8*B)*(a**4*cos(c + d*x) + a**4)*sin(c + d*x)/(6*d) - (3*A - B)*(a**2*cos(c + d*x) + a**2)**2*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_33():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*sec(c + d*x)**3
    F = A*a*(a*cos(c + d*x) + a)**3*tan(c + d*x)*sec(c + d*x)/(2*d) + a**4*x*(8*A + 13*B)/2 - 5*a**4*(A - B)*sin(c + d*x)/(2*d) + a**4*(13*A + 8*B)*atanh(sin(c + d*x))/(2*d) + (5*A + 2*B)*(a**2*cos(c + d*x) + a**2)**2*tan(c + d*x)/(2*d) - (6*A + B)*(a**4*cos(c + d*x) + a**4)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_34():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*sec(c + d*x)**4
    F = A*a*(a*cos(c + d*x) + a)**3*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a**4*x*(A + 4*B) - 5*a**4*(2*A + B)*sin(c + d*x)/(2*d) + a**4*(12*A + 13*B)*atanh(sin(c + d*x))/(2*d) + (2*A + B)*(a**2*cos(c + d*x) + a**2)**2*tan(c + d*x)*sec(c + d*x)/(2*d) + (11*A + 9*B)*(a**4*cos(c + d*x) + a**4)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_35():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*sec(c + d*x)**5
    F = A*a*(a*cos(c + d*x) + a)**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + B*a**4*x + 5*a**4*(7*A + 8*B)*tan(c + d*x)/(8*d) + a**4*(35*A + 48*B)*atanh(sin(c + d*x))/(8*d) + (7*A + 4*B)*(a**2*cos(c + d*x) + a**2)**2*tan(c + d*x)*sec(c + d*x)**2/(12*d) + (35*A + 32*B)*(a**4*cos(c + d*x) + a**4)*tan(c + d*x)*sec(c + d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_36():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*sec(c + d*x)**6
    F = A*a*(a*cos(c + d*x) + a)**3*tan(c + d*x)*sec(c + d*x)**4/(5*d) + 7*a**4*(4*A + 5*B)*atanh(sin(c + d*x))/(8*d) + a**4*(83*A + 100*B)*tan(c + d*x)/(15*d) + a**4*(244*A + 275*B)*tan(c + d*x)*sec(c + d*x)/(120*d) + (8*A + 5*B)*(a**2*cos(c + d*x) + a**2)**2*tan(c + d*x)*sec(c + d*x)**3/(20*d) + (26*A + 25*B)*(a**4*cos(c + d*x) + a**4)*tan(c + d*x)*sec(c + d*x)**2/(30*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_37():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**4*sec(c + d*x)**7
    F = A*a*(a*cos(c + d*x) + a)**3*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 7*a**4*(7*A + 8*B)*tan(c + d*x)*sec(c + d*x)/(16*d) + 7*a**4*(7*A + 8*B)*atanh(sin(c + d*x))/(16*d) + a**4*(72*A + 83*B)*tan(c + d*x)/(15*d) + a**4*(159*A + 176*B)*tan(c + d*x)*sec(c + d*x)**2/(120*d) + (3*A + 2*B)*(a**2*cos(c + d*x) + a**2)**2*tan(c + d*x)*sec(c + d*x)**4/(10*d) + (73*A + 72*B)*(a**4*cos(c + d*x) + a**4)*tan(c + d*x)*sec(c + d*x)**3/(120*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_38():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a*cos(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**4/(d*(a*cos(c + d*x) + a)) - x*(12*A - 15*B)/(8*a) - (4*A - 5*B)*sin(c + d*x)*cos(c + d*x)**3/(4*a*d) - (4*A - 4*B)*sin(c + d*x)**3/(3*a*d) + (4*A - 4*B)*sin(c + d*x)/(a*d) - (12*A - 15*B)*sin(c + d*x)*cos(c + d*x)/(8*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_39():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a*cos(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**3/(d*(a*cos(c + d*x) + a)) + x*(3*A - 3*B)/(2*a) + (3*A - 4*B)*sin(c + d*x)**3/(3*a*d) - (3*A - 4*B)*sin(c + d*x)/(a*d) + (3*A - 3*B)*sin(c + d*x)*cos(c + d*x)/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_40():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a*cos(c + d*x) + a)
    F = B*sin(c + d*x)/(a*d) + x*(A - B)/a - (A - B)*sin(c + d*x)/(a*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_41():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)
    F = B*x/a + (A - B)*sin(c + d*x)/(d*(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_42():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a*cos(c + d*x) + a)
    F = A*atanh(sin(c + d*x))/(a*d) - (A - B)*sin(c + d*x)/(d*(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_43():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a*cos(c + d*x) + a)
    F = -(A - B)*tan(c + d*x)/(d*(a*cos(c + d*x) + a)) - (A - B)*atanh(sin(c + d*x))/(a*d) + (2*A - B)*tan(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_44():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a*cos(c + d*x) + a)
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)/(d*(a*cos(c + d*x) + a)) - (2*A - 2*B)*tan(c + d*x)/(a*d) + (3*A - 2*B)*tan(c + d*x)*sec(c + d*x)/(2*a*d) + (3*A - 2*B)*atanh(sin(c + d*x))/(2*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_45():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**4/(a*cos(c + d*x) + a)
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)**2/(d*(a*cos(c + d*x) + a)) - (3*A - 3*B)*tan(c + d*x)*sec(c + d*x)/(2*a*d) - (3*A - 3*B)*atanh(sin(c + d*x))/(2*a*d) + (4*A - 3*B)*tan(c + d*x)**3/(3*a*d) + (4*A - 3*B)*tan(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_46():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a*cos(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**4/(3*d*(a*cos(c + d*x) + a)**2) + x*(7*A - 10*B)/(2*a**2) + (7*A - 10*B)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) + (7*A - 10*B)*sin(c + d*x)*cos(c + d*x)**3/(3*a**2*d*(cos(c + d*x) + 1)) + (8*A - 12*B)*sin(c + d*x)**3/(3*a**2*d) - (8*A - 12*B)*sin(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_47():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a*cos(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**3/(3*d*(a*cos(c + d*x) + a)**2) - x*(4*A - 7*B)/(2*a**2) - (4*A - 7*B)*sin(c + d*x)*cos(c + d*x)/(2*a**2*d) + (5*A - 8*B)*sin(c + d*x)*cos(c + d*x)**2/(3*a**2*d*(cos(c + d*x) + 1)) + (10*A - 16*B)*sin(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_48():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a*cos(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**2/(3*d*(a*cos(c + d*x) + a)**2) + x*(A - 2*B)/a**2 - (A - 4*B)*sin(c + d*x)/(3*a**2*d) - (A - 2*B)*sin(c + d*x)/(a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_49():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a*cos(c + d*x) + a)**2
    F = B*x/a**2 - (A - B)*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) + (2*A - 5*B)*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_50():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) + (A + 2*B)*sin(c + d*x)/(3*d*(a**2*cos(c + d*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_51():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a*cos(c + d*x) + a)**2
    F = A*atanh(sin(c + d*x))/(a**2*d) - (A - B)*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) - (4*A - B)*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_52():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a*cos(c + d*x) + a)**2
    F = -(A - B)*tan(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) - (2*A - B)*atanh(sin(c + d*x))/(a**2*d) - (2*A - B)*tan(c + d*x)/(a**2*d*(cos(c + d*x) + 1)) + (10*A - 4*B)*tan(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_53():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a*cos(c + d*x) + a)**2
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2) + (7*A - 4*B)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) + (7*A - 4*B)*atanh(sin(c + d*x))/(2*a**2*d) - (8*A - 5*B)*tan(c + d*x)*sec(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1)) - (16*A - 10*B)*tan(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_54():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**4/(a*cos(c + d*x) + a)**2
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)**2/(3*d*(a*cos(c + d*x) + a)**2) - (10*A - 7*B)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) - (10*A - 7*B)*atanh(sin(c + d*x))/(2*a**2*d) - (10*A - 7*B)*tan(c + d*x)*sec(c + d*x)**2/(3*a**2*d*(cos(c + d*x) + 1)) + (12*A - 8*B)*tan(c + d*x)**3/(3*a**2*d) + (12*A - 8*B)*tan(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_55():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**5/(a*cos(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**5/(5*d*(a*cos(c + d*x) + a)**3) + (13*A - 23*B)*sin(c + d*x)*cos(c + d*x)**3/(3*d*(a**3*cos(c + d*x) + a**3)) + (8*A - 13*B)*sin(c + d*x)*cos(c + d*x)**4/(15*a*d*(a*cos(c + d*x) + a)**2) + x*(13*A - 23*B)/(2*a**3) + (13*A - 23*B)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) + (76*A - 136*B)*sin(c + d*x)**3/(15*a**3*d) - (76*A - 136*B)*sin(c + d*x)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_56():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a*cos(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**4/(5*d*(a*cos(c + d*x) + a)**3) + (36*A - 76*B)*sin(c + d*x)*cos(c + d*x)**2/(15*d*(a**3*cos(c + d*x) + a**3)) + (6*A - 11*B)*sin(c + d*x)*cos(c + d*x)**3/(15*a*d*(a*cos(c + d*x) + a)**2) - x*(6*A - 13*B)/(2*a**3) - (6*A - 13*B)*sin(c + d*x)*cos(c + d*x)/(2*a**3*d) + (72*A - 152*B)*sin(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_57():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a*cos(c + d*x) + a)**3
    F = -(A - 3*B)*sin(c + d*x)/(d*(a**3*cos(c + d*x) + a**3)) + (A - B)*sin(c + d*x)*cos(c + d*x)**3/(5*d*(a*cos(c + d*x) + a)**3) + (4*A - 9*B)*sin(c + d*x)*cos(c + d*x)**2/(15*a*d*(a*cos(c + d*x) + a)**2) + x*(A - 3*B)/a**3 - (7*A - 27*B)*sin(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_58():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a*cos(c + d*x) + a)**3
    F = B*x/a**3 + (A - B)*sin(c + d*x)*cos(c + d*x)**2/(5*d*(a*cos(c + d*x) + a)**3) + (4*A - 29*B)*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) - (2*A - 7*B)*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_59():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a*cos(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) + (3*A + 7*B)*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) + (3*A - 8*B)*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_60():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) + (2*A + 3*B)*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) + (2*A + 3*B)*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_61():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a*cos(c + d*x) + a)**3
    F = A*atanh(sin(c + d*x))/(a**3*d) - (A - B)*sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) - (22*A - 2*B)*sin(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) - (7*A - 2*B)*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_62():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a*cos(c + d*x) + a)**3
    F = -(A - B)*tan(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) - (3*A - B)*tan(c + d*x)/(d*(a**3*cos(c + d*x) + a**3)) - (9*A - 4*B)*tan(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2) - (3*A - B)*atanh(sin(c + d*x))/(a**3*d) + (72*A - 22*B)*tan(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_63():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a*cos(c + d*x) + a)**3
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3) - (76*A - 36*B)*tan(c + d*x)*sec(c + d*x)/(15*d*(a**3*cos(c + d*x) + a**3)) - (11*A - 6*B)*tan(c + d*x)*sec(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2) + (13*A - 6*B)*tan(c + d*x)*sec(c + d*x)/(2*a**3*d) + (13*A - 6*B)*atanh(sin(c + d*x))/(2*a**3*d) - (152*A - 72*B)*tan(c + d*x)/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_64():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**5/(a*cos(c + d*x) + a)**4
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**5/(7*d*(a*cos(c + d*x) + a)**4) + (A - 2*B)*sin(c + d*x)*cos(c + d*x)**4/(5*a*d*(a*cos(c + d*x) + a)**3) - x*(8*A - 21*B)/(2*a**4) - (8*A - 21*B)*sin(c + d*x)*cos(c + d*x)/(2*a**4*d) + (52*A - 129*B)*sin(c + d*x)*cos(c + d*x)**3/(105*a**4*d*(cos(c + d*x) + 1)**2) + (332*A - 864*B)*sin(c + d*x)*cos(c + d*x)**2/(105*a**4*d*(cos(c + d*x) + 1)) + (664*A - 1728*B)*sin(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_65():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a*cos(c + d*x) + a)**4
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**4/(7*d*(a*cos(c + d*x) + a)**4) + (5*A - 12*B)*sin(c + d*x)*cos(c + d*x)**3/(35*a*d*(a*cos(c + d*x) + a)**3) + x*(A - 4*B)/a**4 - (A - 4*B)*sin(c + d*x)/(a**4*d*(cos(c + d*x) + 1)) + (25*A - 88*B)*sin(c + d*x)*cos(c + d*x)**2/(105*a**4*d*(cos(c + d*x) + 1)**2) - (55*A - 244*B)*sin(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_66():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a*cos(c + d*x) + a)**4
    F = B*x/a**4 + (A - B)*sin(c + d*x)*cos(c + d*x)**3/(7*d*(a*cos(c + d*x) + a)**4) + (3*A - 10*B)*sin(c + d*x)*cos(c + d*x)**2/(35*a*d*(a*cos(c + d*x) + a)**3) - (6*A - 55*B)*sin(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1)**2) + (12*A - 215*B)*sin(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_67():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a*cos(c + d*x) + a)**4
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**2/(7*d*(a*cos(c + d*x) + a)**4) - (A - 8*B)*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3) - (2*A + 54*B)*sin(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1)**2) + (13*A + 36*B)*sin(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_68():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a*cos(c + d*x) + a)**4
    F = -(A - B)*sin(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) + (8*A + 13*B)*sin(c + d*x)/(105*d*(a**4*cos(c + d*x) + a**4)) + (8*A + 13*B)*sin(c + d*x)/(105*d*(a**2*cos(c + d*x) + a**2)**2) + (4*A - 11*B)*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_69():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)**4
    F = (A - B)*sin(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) + (6*A + 8*B)*sin(c + d*x)/(105*d*(a**4*cos(c + d*x) + a**4)) + (6*A + 8*B)*sin(c + d*x)/(105*d*(a**2*cos(c + d*x) + a**2)**2) + (3*A + 4*B)*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_70():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a*cos(c + d*x) + a)**4
    F = A*atanh(sin(c + d*x))/(a**4*d) - (A - B)*sin(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) - (10*A - 3*B)*sin(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3) - (55*A - 6*B)*sin(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1)**2) - (160*A - 6*B)*sin(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_71():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a*cos(c + d*x) + a)**4
    F = -(A - B)*tan(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) - (12*A - 5*B)*tan(c + d*x)/(35*a*d*(a*cos(c + d*x) + a)**3) - (4*A - B)*atanh(sin(c + d*x))/(a**4*d) - (4*A - B)*tan(c + d*x)/(a**4*d*(cos(c + d*x) + 1)) - (88*A - 25*B)*tan(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1)**2) + (664*A - 160*B)*tan(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_72():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a*cos(c + d*x) + a)**4
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)/(7*d*(a*cos(c + d*x) + a)**4) - (2*A - B)*tan(c + d*x)*sec(c + d*x)/(5*a*d*(a*cos(c + d*x) + a)**3) + (21*A - 8*B)*tan(c + d*x)*sec(c + d*x)/(2*a**4*d) + (21*A - 8*B)*atanh(sin(c + d*x))/(2*a**4*d) - (129*A - 52*B)*tan(c + d*x)*sec(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1)**2) - (864*A - 332*B)*tan(c + d*x)*sec(c + d*x)/(105*a**4*d*(cos(c + d*x) + 1)) - (1728*A - 664*B)*tan(c + d*x)/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_73():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**3
    F = 2*B*a*sin(c + d*x)*cos(c + d*x)**4/(9*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(9*A + 8*B)*sin(c + d*x)*cos(c + d*x)**3/(63*d*sqrt(a*cos(c + d*x) + a)) + 4*a*(9*A + 8*B)*sin(c + d*x)/(45*d*sqrt(a*cos(c + d*x) + a)) - (72*A + 64*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(315*d) + (36*A + 32*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(105*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_74():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**2
    F = 2*B*a*sin(c + d*x)*cos(c + d*x)**3/(7*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(7*A + 6*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) - (28*A + 24*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(105*d) + (14*A + 12*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_75():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)
    F = 2*B*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*a*d) + 2*a*(5*A + 7*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) + (10*A - 4*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_76():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)
    F = 2*B*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d) + 2*a*(3*A + B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_77():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)
    F = 2*A*sqrt(a)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 2*B*a*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_78():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**2
    F = A*a*tan(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)) + sqrt(a)*(A + 2*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_79():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**3
    F = A*a*tan(c + d*x)*sec(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)) + sqrt(a)*(3*A + 4*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a*(3*A + 4*B)*tan(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_80():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**4
    F = A*a*tan(c + d*x)*sec(c + d*x)**2/(3*d*sqrt(a*cos(c + d*x) + a)) + sqrt(a)*(5*A + 6*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a*(5*A + 6*B)*tan(c + d*x)*sec(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)) + a*(5*A + 6*B)*tan(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_81():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3
    F = 2*B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**4/(11*d) + 2*a**2*(11*A + 12*B)*sin(c + d*x)*cos(c + d*x)**4/(99*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(187*A + 168*B)*sin(c + d*x)*cos(c + d*x)**3/(693*d*sqrt(a*cos(c + d*x) + a)) + 4*a**2*(187*A + 168*B)*sin(c + d*x)/(495*d*sqrt(a*cos(c + d*x) + a)) - 8*a*(187*A + 168*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3465*d) + (748*A + 672*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(1155*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_82():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2
    F = 2*B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**3/(9*d) + 2*a**2*(9*A + 10*B)*sin(c + d*x)*cos(c + d*x)**3/(63*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(39*A + 34*B)*sin(c + d*x)/(45*d*sqrt(a*cos(c + d*x) + a)) - 4*a*(39*A + 34*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(315*d) + (78*A + 68*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(105*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_83():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)
    F = 2*B*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(7*a*d) + 8*a**2*(21*A + 19*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(21*A + 19*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(105*d) + (14*A - 4*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_84():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*B*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 8*a**2*(5*A + 3*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(5*A + 3*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_85():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)
    F = 2*A*a**(sympy.S(3)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 2*B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d) + 2*a**2*(3*A + 4*B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_86():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2
    F = A*a*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)/d + a**(sympy.S(3)/2)*(3*A + 2*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d - a**2*(A - 2*B)*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_87():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3
    F = A*a*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)/(2*d) + a**(sympy.S(3)/2)*(7*A + 12*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a**2*(5*A + 4*B)*tan(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_88():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4
    F = A*a*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a**(sympy.S(3)/2)*(11*A + 14*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a**2*(7*A + 6*B)*tan(c + d*x)*sec(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)) + a**2*(11*A + 14*B)*tan(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_89():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**5
    F = A*a*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a**(sympy.S(3)/2)*(75*A + 88*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(64*d) + a**2*(9*A + 8*B)*tan(c + d*x)*sec(c + d*x)**2/(24*d*sqrt(a*cos(c + d*x) + a)) + a**2*(75*A + 88*B)*tan(c + d*x)*sec(c + d*x)/(96*d*sqrt(a*cos(c + d*x) + a)) + a**2*(75*A + 88*B)*tan(c + d*x)/(64*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_90():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2
    F = 2*B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**3/(11*d) + 2*a**3*(209*A + 194*B)*sin(c + d*x)*cos(c + d*x)**3/(693*d*sqrt(a*cos(c + d*x) + a)) + 2*a**3*(803*A + 710*B)*sin(c + d*x)/(495*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(11*A + 14*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**3/(99*d) - 4*a**2*(803*A + 710*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3465*d) + 2*a*(803*A + 710*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(1155*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_91():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)
    F = 2*B*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sin(c + d*x)/(9*a*d) + 64*a**3*(15*A + 13*B)*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)) + 16*a**2*(15*A + 13*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(315*d) + 2*a*(15*A + 13*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(105*d) + (18*A - 4*B)*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_92():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*B*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sin(c + d*x)/(7*d) + 64*a**3*(7*A + 5*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)) + 16*a**2*(7*A + 5*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(105*d) + 2*a*(7*A + 5*B)*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_93():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)
    F = 2*A*a**(sympy.S(5)/2)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 2*B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + 2*a**3*(35*A + 32*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(5*A + 8*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_94():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**2
    F = A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)/d + a**(sympy.S(5)/2)*(5*A + 2*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + a**3*(3*A + 14*B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)) - a**2*(3*A - 2*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_95():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**3
    F = A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)*sec(c + d*x)/(2*d) + a**(sympy.S(5)/2)*(19*A + 20*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) - a**3*(9*A - 4*B)*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)) + a**2*(7*A + 4*B)*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_96():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4
    F = A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a**(sympy.S(5)/2)*(25*A + 38*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a**3*(49*A + 54*B)*tan(c + d*x)/(24*d*sqrt(a*cos(c + d*x) + a)) + a**2*(3*A + 2*B)*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_97():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**5
    F = A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a**(sympy.S(5)/2)*(163*A + 200*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(64*d) + a**3*(95*A + 104*B)*tan(c + d*x)*sec(c + d*x)/(96*d*sqrt(a*cos(c + d*x) + a)) + a**3*(163*A + 200*B)*tan(c + d*x)/(64*d*sqrt(a*cos(c + d*x) + a)) + a**2*(11*A + 8*B)*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**2/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_98():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**6
    F = A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a**(sympy.S(5)/2)*(283*A + 326*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(128*d) + a**3*(157*A + 170*B)*tan(c + d*x)*sec(c + d*x)**2/(240*d*sqrt(a*cos(c + d*x) + a)) + a**3*(283*A + 326*B)*tan(c + d*x)*sec(c + d*x)/(192*d*sqrt(a*cos(c + d*x) + a)) + a**3*(283*A + 326*B)*tan(c + d*x)/(128*d*sqrt(a*cos(c + d*x) + a)) + a**2*(13*A + 10*B)*sqrt(a*cos(c + d*x) + a)*tan(c + d*x)*sec(c + d*x)**3/(40*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_99():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/sqrt(a*cos(c + d*x) + a)
    F = 2*B*sin(c + d*x)*cos(c + d*x)**3/(7*d*sqrt(a*cos(c + d*x) + a)) + (14*A - 2*B)*sin(c + d*x)*cos(c + d*x)**2/(35*d*sqrt(a*cos(c + d*x) + a)) + (196*A - 148*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)) - (14*A - 62*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(105*a*d) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_100():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/sqrt(a*cos(c + d*x) + a)
    F = 2*B*sin(c + d*x)*cos(c + d*x)**2/(5*d*sqrt(a*cos(c + d*x) + a)) - (20*A - 28*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)) + (10*A - 2*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*a*d) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_101():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/sqrt(a*cos(c + d*x) + a)
    F = 2*B*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*a*d) + (6*A - 4*B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_102():
    f = (A + B*cos(c + d*x))/sqrt(a*cos(c + d*x) + a)
    F = 2*B*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_103():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/sqrt(a*cos(c + d*x) + a)
    F = 2*A*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_104():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/sqrt(a*cos(c + d*x) + a)
    F = A*tan(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)) - (A - 2*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_105():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/sqrt(a*cos(c + d*x) + a)
    F = A*tan(c + d*x)*sec(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)) - (A - 4*B)*tan(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(sqrt(a)*d) + (7*A - 4*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_106():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**4/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - (7*A - 11*B)*sin(c + d*x)*cos(c + d*x)**3/(14*a*d*sqrt(a*cos(c + d*x) + a)) + (63*A - 67*B)*sin(c + d*x)*cos(c + d*x)**2/(70*a*d*sqrt(a*cos(c + d*x) + a)) + (651*A - 799*B)*sin(c + d*x)/(105*a*d*sqrt(a*cos(c + d*x) + a)) - (273*A - 397*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(210*a**2*d) - sqrt(2)*(15*A - 19*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_107():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**3/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - (5*A - 9*B)*sin(c + d*x)*cos(c + d*x)**2/(10*a*d*sqrt(a*cos(c + d*x) + a)) - (65*A - 93*B)*sin(c + d*x)/(15*a*d*sqrt(a*cos(c + d*x) + a)) + (35*A - 39*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(30*a**2*d) + sqrt(2)*(11*A - 15*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_108():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**2/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (9*A - 13*B)*sin(c + d*x)/(3*a*d*sqrt(a*cos(c + d*x) + a)) - (3*A - 7*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(6*a**2*d) - sqrt(2)*(7*A - 11*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_109():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*B*sin(c + d*x)/(a*d*sqrt(a*cos(c + d*x) + a)) - (A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A - 7*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_110():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(A + 3*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_111():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*A*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - (A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(5*A - B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_112():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*tan(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (3*A - B)*tan(c + d*x)/(2*a*d*sqrt(a*cos(c + d*x) + a)) - (3*A - 2*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + sqrt(2)*(9*A - 5*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_113():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (2*A - B)*tan(c + d*x)*sec(c + d*x)/(2*a*d*sqrt(a*cos(c + d*x) + a)) - (7*A - 6*B)*tan(c + d*x)/(4*a*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(13*A - 9*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d) + (19*A - 12*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_114():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**4/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (13*A - 21*B)*sin(c + d*x)*cos(c + d*x)**3/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - (85*A - 157*B)*sin(c + d*x)*cos(c + d*x)**2/(80*a**2*d*sqrt(a*cos(c + d*x) + a)) - (985*A - 1729*B)*sin(c + d*x)/(120*a**2*d*sqrt(a*cos(c + d*x) + a)) + (475*A - 787*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(240*a**3*d) + sqrt(2)*(163*A - 283*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_115():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**3/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (9*A - 17*B)*sin(c + d*x)*cos(c + d*x)**2/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (93*A - 197*B)*sin(c + d*x)/(24*a**2*d*sqrt(a*cos(c + d*x) + a)) - (39*A - 95*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(48*a**3*d) - sqrt(2)*(75*A - 163*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_116():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**2/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (5*A - 13*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - (A - 9*B)*sin(c + d*x)/(4*a**2*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(19*A - 75*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_117():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (5*A - 13*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(5*A + 19*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_118():
    f = (A + B*cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (3*A + 5*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A + 5*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_119():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*A*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - (A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (11*A - 3*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*(43*A - 3*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_120():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*tan(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (15*A - 7*B)*tan(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (35*A - 11*B)*tan(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) - (5*A - 2*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + sqrt(2)*(115*A - 43*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_121():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*tan(c + d*x)*sec(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (19*A - 11*B)*tan(c + d*x)*sec(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (31*A - 15*B)*tan(c + d*x)*sec(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) - (63*A - 35*B)*tan(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) + (39*A - 20*B)*atanh(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*a**(sympy.S(5)/2)*d) - sqrt(2)*(219*A - 115*B)*atanh(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_122():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*B*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + 2*a*(A + B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 10*a*(A + B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 10*a*(A + B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*a*(9*A + 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 2*a*(9*A + 7*B)*elliptic_e(c/2 + d*x/2, 2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_123():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + 2*a*(A + B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*a*(A + B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a*(7*A + 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 2*a*(7*A + 5*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_124():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))
    F = 2*B*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a*(A + B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*(A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*a*(5*A + 3*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_125():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)/sqrt(cos(c + d*x))
    F = 2*B*a*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*a*(A + B)*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*(3*A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_126():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*a*(A - B)*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*(A + B)*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_127():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(A + B)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*a*(A + B)*elliptic_e(c/2 + d*x/2, 2)/d + 2*a*(A + 3*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_128():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a*(A + B)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*a*(3*A + 5*B)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 2*a*(3*A + 5*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_129():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(9*d) + 4*a**2*(6*A + 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 4*a**2*(6*A + 5*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 4*a**2*(9*A + 8*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 4*a**2*(9*A + 8*B)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 2*a**2*(9*A + 11*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_130():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))
    F = 2*B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(7*d) + 4*a**2*(4*A + 3*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*(7*A + 6*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 4*a**2*(7*A + 6*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*a**2*(7*A + 9*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_131():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2/sqrt(cos(c + d*x))
    F = 2*B*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 4*a**2*(2*A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 4*a**2*(5*A + 4*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a**2*(5*A + 7*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_132():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 4*B*a**2*elliptic_e(c/2 + d*x/2, 2)/d - 2*a**2*(3*A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 4*a**2*(3*A + 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_133():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(5)/2)
    F = -4*A*a**2*elliptic_e(c/2 + d*x/2, 2)/d + 2*A*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(2*A + 3*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*a**2*(5*A + 3*B)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_134():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 4*a**2*(A + 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 4*a**2*(4*A + 5*B)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 4*a**2*(4*A + 5*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 2*a**2*(7*A + 5*B)*sin(c + d*x)/(15*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_135():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*(a**2*cos(c + d*x) + a**2)*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 4*a**2*(3*A + 4*B)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 4*a**2*(3*A + 4*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**2*(6*A + 7*B)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(6*A + 7*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 2*a**2*(9*A + 7*B)*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_136():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(11*d) + 4*a**3*(17*A + 15*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + 4*a**3*(17*A + 15*B)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 20*a**3*(22*A + 21*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(693*d) + 4*a**3*(121*A + 105*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(231*d) + 4*a**3*(121*A + 105*B)*elliptic_f(c/2 + d*x/2, 2)/(231*d) + (22*A + 30*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(99*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_137():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sqrt(cos(c + d*x))
    F = 2*B*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(9*d) + 4*a**3*(13*A + 11*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + 4*a**3*(13*A + 11*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 4*a**3*(21*A + 17*B)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 4*a**3*(24*A + 23*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(105*d) + (18*A + 26*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_138():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3/sqrt(cos(c + d*x))
    F = 2*B*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)*sqrt(cos(c + d*x))/(7*d) + 4*a**3*(9*A + 7*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**3*(21*A + 13*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 4*a**3*(42*A + 41*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(105*d) + (14*A + 22*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)*sqrt(cos(c + d*x))/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_139():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 4*a**3*(5*A - 6*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d) + 4*a**3*(5*A + 3*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 4*a**3*(5*A + 9*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) - (10*A - 2*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_140():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) - 4*a**3*(A - B)*elliptic_e(c/2 + d*x/2, 2)/d + 20*a**3*(A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) - 4*a**3*(4*A + B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + (14*A + 6*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_141():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 4*a**3*(3*A + 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*d) - 4*a**3*(9*A + 5*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**3*(21*A + 20*B)*sin(c + d*x)/(15*d*sqrt(cos(c + d*x))) + (18*A + 10*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)/(15*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_142():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 4*a**3*(7*A + 9*B)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - 4*a**3*(7*A + 9*B)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + 4*a**3*(13*A + 21*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 4*a**3*(41*A + 42*B)*sin(c + d*x)/(105*d*cos(c + d*x)**(sympy.S(3)/2)) + (22*A + 14*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_143():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3/cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**2*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2)) + 4*a**3*(11*A + 13*B)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(3)/2)) + 4*a**3*(11*A + 13*B)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 4*a**3*(17*A + 21*B)*sin(c + d*x)/(15*d*sqrt(cos(c + d*x))) - 4*a**3*(17*A + 21*B)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + 4*a**3*(23*A + 24*B)*sin(c + d*x)/(105*d*cos(c + d*x)**(sympy.S(5)/2)) + (26*A + 18*B)*(a**3*cos(c + d*x) + a**3)*sin(c + d*x)/(63*d*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_144():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(d*(a*cos(c + d*x) + a)) + (-15*A + 21*B)*elliptic_e(c/2 + d*x/2, 2)/(5*a*d) - (5*A - 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*a*d) + (5*A - 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) + (5*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_145():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(d*(a*cos(c + d*x) + a)) - (3*A - 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d) - (3*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a*d) + (3*A - 3*B)*elliptic_e(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_146():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)
    F = (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*cos(c + d*x) + a)) - (A - 3*B)*elliptic_e(c/2 + d*x/2, 2)/(a*d) + (A - B)*elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_147():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a*cos(c + d*x) + a)) + (A - B)*elliptic_e(c/2 + d*x/2, 2)/(a*d) + (A + B)*elliptic_f(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_148():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(d*(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - (A - B)*elliptic_f(c/2 + d*x/2, 2)/(a*d) + (3*A - B)*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) - (3*A - B)*elliptic_e(c/2 + d*x/2, 2)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_149():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(d*(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - (3*A - 3*B)*sin(c + d*x)/(a*d*sqrt(cos(c + d*x))) + (3*A - 3*B)*elliptic_e(c/2 + d*x/2, 2)/(a*d) + (5*A - 3*B)*sin(c + d*x)/(3*a*d*cos(c + d*x)**(sympy.S(3)/2)) + (5*A - 3*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_150():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(3*d*(a*cos(c + d*x) + a)**2) + (-35*A + 56*B)*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d) + (2*A - 3*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(a**2*d*(cos(c + d*x) + 1)) + (10*A - 15*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d) + (10*A - 15*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - (35*A - 56*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_151():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*d*(a*cos(c + d*x) + a)**2) + (4*A - 7*B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + (4*A - 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(cos(c + d*x) + 1)) - (5*A - 10*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d) - (5*A - 10*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_152():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**2
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d*(a*cos(c + d*x) + a)**2) - (A - 4*B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + (2*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) + (2*A - 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a**2*d*(cos(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_153():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**2
    F = -B*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) + B*sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d*(cos(c + d*x) + 1)) + (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*cos(c + d*x) + a)**2) + (A + 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_154():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x)))
    F = A*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - A*sin(c + d*x)*sqrt(cos(c + d*x))/(a**2*d*(cos(c + d*x) + 1)) - (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*(a*cos(c + d*x) + a)**2) + (2*A + B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_155():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))) + (4*A - B)*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) - (4*A - B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - (5*A - 2*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) - (5*A - 2*B)*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_156():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(3*d*(a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) - (7*A - 4*B)*sin(c + d*x)/(a**2*d*sqrt(cos(c + d*x))) + (7*A - 4*B)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d) - (7*A - 4*B)*sin(c + d*x)/(3*a**2*d*(cos(c + d*x) + 1)*cos(c + d*x)**(sympy.S(3)/2)) + (10*A - 5*B)*sin(c + d*x)/(3*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + (10*A - 5*B)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_157():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(9)/2)/(a*cos(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(9)/2)/(5*d*(a*cos(c + d*x) + a)**3) + (33*A - 63*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(10*d*(a**3*cos(c + d*x) + a**3)) + (7*A - 12*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(15*a*d*(a*cos(c + d*x) + a)**2) + (-119*A + 231*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + (11*A - 21*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*a**3*d) + (11*A - 21*B)*elliptic_f(c/2 + d*x/2, 2)/(2*a**3*d) - (119*A - 231*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(30*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_158():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(5*d*(a*cos(c + d*x) + a)**3) + (49*A - 119*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(30*d*(a**3*cos(c + d*x) + a**3)) + (A - 2*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*a*d*(a*cos(c + d*x) + a)**2) - (13*A - 33*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(6*a**3*d) - (13*A - 33*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) + (49*A - 119*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_159():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(5*d*(a*cos(c + d*x) + a)**3) + (3*A - 13*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a**3*cos(c + d*x) + a**3)) + (3*A - 8*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(15*a*d*(a*cos(c + d*x) + a)**2) + (3*A - 13*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) - (9*A - 49*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_160():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**3
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d*(a*cos(c + d*x) + a)**3) + (A + 9*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) + (A - 6*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) + (A + 3*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) - (A + 9*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_161():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) + (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*cos(c + d*x) + a)**3) + (A + 4*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) + (A - B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d) + (A + B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_162():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**3*sqrt(cos(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d*(a*cos(c + d*x) + a)**3) - (9*A + B)*sin(c + d*x)*sqrt(cos(c + d*x))/(10*d*(a**3*cos(c + d*x) + a**3)) - (6*A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*a*d*(a*cos(c + d*x) + a)**2) + (3*A + B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) + (9*A + B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_163():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3*sqrt(cos(c + d*x))) - (13*A - 3*B)*sin(c + d*x)/(6*d*(a**3*cos(c + d*x) + a**3)*sqrt(cos(c + d*x))) - (8*A - 3*B)*sin(c + d*x)/(15*a*d*(a*cos(c + d*x) + a)**2*sqrt(cos(c + d*x))) - (13*A - 3*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) + (49*A - 9*B)*sin(c + d*x)/(10*a**3*d*sqrt(cos(c + d*x))) - (49*A - 9*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_164():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(5*d*(a*cos(c + d*x) + a)**3*cos(c + d*x)**(sympy.S(3)/2)) - (119*A - 49*B)*sin(c + d*x)/(30*d*(a**3*cos(c + d*x) + a**3)*cos(c + d*x)**(sympy.S(3)/2)) - (2*A - B)*sin(c + d*x)/(3*a*d*(a*cos(c + d*x) + a)**2*cos(c + d*x)**(sympy.S(3)/2)) + (33*A - 13*B)*sin(c + d*x)/(6*a**3*d*cos(c + d*x)**(sympy.S(3)/2)) + (33*A - 13*B)*elliptic_f(c/2 + d*x/2, 2)/(6*a**3*d) - (119*A - 49*B)*sin(c + d*x)/(10*a**3*d*sqrt(cos(c + d*x))) + (119*A - 49*B)*elliptic_e(c/2 + d*x/2, 2)/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_165():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)
    F = B*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(4*d*sqrt(a*cos(c + d*x) + a)) + 5*sqrt(a)*(8*A + 7*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(64*d) + a*(8*A + 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(24*d*sqrt(a*cos(c + d*x) + a)) + 5*a*(8*A + 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(96*d*sqrt(a*cos(c + d*x) + a)) + 5*a*(8*A + 7*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(64*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_166():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)
    F = B*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) + sqrt(a)*(6*A + 5*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a*(6*A + 5*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(12*d*sqrt(a*cos(c + d*x) + a)) + a*(6*A + 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_167():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))
    F = B*a*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(a*cos(c + d*x) + a)) + sqrt(a)*(4*A + 3*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a*(4*A + 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_168():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)/sqrt(cos(c + d*x))
    F = B*a*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(a*cos(c + d*x) + a)) + sqrt(a)*(2*A + B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_169():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*B*sqrt(a)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_170():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(2*A + 3*B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_171():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 4*a*(4*A + 5*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a*(4*A + 5*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_172():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sin(c + d*x)/(7*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + 16*a*(6*A + 7*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 8*a*(6*A + 7*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(6*A + 7*B)*sin(c + d*x)/(35*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_173():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d) + a**(sympy.S(3)/2)*(88*A + 75*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(64*d) + a**2*(8*A + 9*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(24*d*sqrt(a*cos(c + d*x) + a)) + a**2*(88*A + 75*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(96*d*sqrt(a*cos(c + d*x) + a)) + a**2*(88*A + 75*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(64*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_174():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(3*d) + a**(sympy.S(3)/2)*(14*A + 11*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a**2*(6*A + 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(12*d*sqrt(a*cos(c + d*x) + a)) + a**2*(14*A + 11*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(8*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_175():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d) + a**(sympy.S(3)/2)*(12*A + 7*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) + a**2*(4*A + 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_176():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + a**(sympy.S(3)/2)*(2*A + 3*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d - a**2*(2*A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_177():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*B*a**(sympy.S(3)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 2*a**2*(4*A + 3*B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_178():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(6*A + 5*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(18*A + 25*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_179():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 2*a**2*(8*A + 7*B)*sin(c + d*x)/(35*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 4*a**2*(52*A + 63*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*(52*A + 63*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_180():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2)) + 2*a**2*(10*A + 9*B)*sin(c + d*x)/(63*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + 16*a**2*(34*A + 39*B)*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 8*a**2*(34*A + 39*B)*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(34*A + 39*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_181():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(5*d) + a**(sympy.S(5)/2)*(326*A + 283*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(128*d) + a**3*(170*A + 157*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(240*d*sqrt(a*cos(c + d*x) + a)) + a**3*(326*A + 283*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(192*d*sqrt(a*cos(c + d*x) + a)) + a**3*(326*A + 283*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(128*d*sqrt(a*cos(c + d*x) + a)) + a**2*(10*A + 13*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(40*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_182():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(4*d) + a**(sympy.S(5)/2)*(200*A + 163*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(64*d) + a**3*(104*A + 95*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(96*d*sqrt(a*cos(c + d*x) + a)) + a**3*(200*A + 163*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(64*d*sqrt(a*cos(c + d*x) + a)) + a**2*(8*A + 11*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_183():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + a**(sympy.S(5)/2)*(38*A + 25*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(8*d) + a**3*(54*A + 49*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(24*d*sqrt(a*cos(c + d*x) + a)) + a**2*(2*A + 3*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_184():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + a**(sympy.S(5)/2)*(20*A + 19*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*d) - a**3*(4*A - 9*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a)) - a**2*(4*A - B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_185():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + a**(sympy.S(5)/2)*(2*A + 5*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d - a**3*(14*A + 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(2*A + B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_186():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*B*a**(sympy.S(5)/2)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/d + 2*a**3*(32*A + 35*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*(8*A + 5*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(15*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_187():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 2*a**3*(10*A + 11*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**3*(230*A + 301*B)*sin(c + d*x)/(105*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**2*(10*A + 7*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_188():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2)) + 2*a**3*(124*A + 135*B)*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 4*a**3*(292*A + 345*B)*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 2*a**3*(292*A + 345*B)*sin(c + d*x)/(315*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(4*A + 3*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_189():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(13)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(11*d*cos(c + d*x)**(sympy.S(11)/2)) + 2*a**3*(194*A + 209*B)*sin(c + d*x)/(693*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2)) + 16*a**3*(710*A + 803*B)*sin(c + d*x)/(3465*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + 8*a**3*(710*A + 803*B)*sin(c + d*x)/(3465*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**3*(710*A + 803*B)*sin(c + d*x)/(1155*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(14*A + 11*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(99*d*cos(c + d*x)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_190():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/sqrt(a*cos(c + d*x) + a)
    F = B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*sqrt(a*cos(c + d*x) + a)) + (4*A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d) - (4*A - 7*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_191():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/sqrt(a*cos(c + d*x) + a)
    F = B*sin(c + d*x)*sqrt(cos(c + d*x))/(d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d) + (2*A - B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_192():
    f = (A + B*cos(c + d*x))/(sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x)))
    F = 2*B*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(sqrt(a)*d) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_193():
    f = (A + B*cos(c + d*x))/(sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*A*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_194():
    f = (A + B*cos(c + d*x))/(sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*A*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - (2*A - 6*B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_195():
    f = (A + B*cos(c + d*x))/(sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(7)/2))
    F = 2*A*sin(c + d*x)/(5*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(5)/2)) - (2*A - 10*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) + (26*A - 10*B)*sin(c + d*x)/(15*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(A - B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_196():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - (A - 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*a*d*sqrt(a*cos(c + d*x) + a)) + (2*A - 3*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) - sqrt(2)*(5*A - 9*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_197():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*B*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(3)/2)*d) + (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(A - 5*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_198():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A + B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_199():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + (5*A - B)*sin(c + d*x)/(2*a*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(7*A - 3*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_200():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + (7*A - 3*B)*sin(c + d*x)/(6*a*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - (19*A - 15*B)*sin(c + d*x)/(6*a*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(11*A - 7*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_201():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (7*A - 15*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - (11*A - 35*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) + (2*A - 5*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) - sqrt(2)*(43*A - 115*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_202():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*B*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(5)/2)*d) + (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (3*A - 11*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(3*A - 43*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_203():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (A + 7*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(5*A + 3*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_204():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (9*A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(19*A + 5*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_205():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))) - (13*A - 5*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + (49*A - 9*B)*sin(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(75*A - 19*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_206():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)) - (17*A - 9*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + (95*A - 39*B)*sin(c + d*x)/(48*a**2*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - (299*A - 147*B)*sin(c + d*x)/(48*a**2*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(163*A - 75*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_207():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) + (3*A - 7*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (79*A - 259*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) - (49*A - 189*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(64*a**3*d*sqrt(a*cos(c + d*x) + a)) + (2*A - 7*B)*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(7)/2)*d) - sqrt(2)*(177*A - 637*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_208():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = 2*B*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))/(a**(sympy.S(7)/2)*d) + (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) + (5*A - 17*B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (5*A - 49*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(5*A - 177*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_209():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = (A - B)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) + (A - 13*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) + (17*A + 67*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(7*A + 5*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_210():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = (A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) + (A + 3*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (5*A - 17*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(13*A + 7*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_211():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(cos(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - (5*A - B)*sin(c + d*x)*sqrt(cos(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (103*A + 5*B)*sin(c + d*x)*sqrt(cos(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*(63*A + 13*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_212():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(cos(c + d*x))) - (19*A - 7*B)*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(cos(c + d*x))) - (199*A - 43*B)*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(cos(c + d*x))) + (691*A - 103*B)*sin(c + d*x)/(192*a**3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) - sqrt(2)*(363*A - 63*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_213():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = -(A - B)*sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**(sympy.S(3)/2)) - (23*A - 11*B)*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)) - (109*A - 41*B)*sin(c + d*x)/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)) + (579*A - 199*B)*sin(c + d*x)/(192*a**3*d*sqrt(a*cos(c + d*x) + a)*cos(c + d*x)**(sympy.S(3)/2)) - (1887*A - 691*B)*sin(c + d*x)/(192*a**3*d*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))) + sqrt(2)*(1015*A - 363*B)*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_214():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*cos(c + d*x)**2
    F = B*b*sin(c + d*x)*cos(c + d*x)**3/(4*d) + x*(A*a/2 + 3*B*b/8) + (4*A*a + 3*B*b)*sin(c + d*x)*cos(c + d*x)/(8*d) - (A*b + B*a)*sin(c + d*x)**3/(3*d) + (A*b + B*a)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_215():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*cos(c + d*x)
    F = B*b*sin(c + d*x)*cos(c + d*x)**2/(3*d) + x*(A*b/2 + B*a/2) + (3*A*a + 2*B*b)*sin(c + d*x)/(3*d) + (A*b + B*a)*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_216():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))
    F = B*b*sin(c + d*x)*cos(c + d*x)/(2*d) + x*(A*a + B*b/2) + (A*b + B*a)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_217():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)
    F = A*a*atanh(sin(c + d*x))/d + B*b*sin(c + d*x)/d + x*(A*b + B*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_218():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)**2
    F = A*a*tan(c + d*x)/d + B*b*x + (A*b + B*a)*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_219():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)**3
    F = A*a*tan(c + d*x)*sec(c + d*x)/(2*d) + (A*a + 2*B*b)*atanh(sin(c + d*x))/(2*d) + (A*b + B*a)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_220():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)**4
    F = A*a*tan(c + d*x)*sec(c + d*x)**2/(3*d) + (2*A*a + 3*B*b)*tan(c + d*x)/(3*d) + (A*b + B*a)*tan(c + d*x)*sec(c + d*x)/(2*d) + (A*b + B*a)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_221():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)**5
    F = A*a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + (3*A*a + 4*B*b)*tan(c + d*x)*sec(c + d*x)/(8*d) + (3*A*a + 4*B*b)*atanh(sin(c + d*x))/(8*d) + (A*b + B*a)*tan(c + d*x)**3/(3*d) + (A*b + B*a)*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_222():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*cos(c + d*x)**2
    F = B*b*(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**3/(5*d) + b*(5*A*b + 6*B*a)*sin(c + d*x)*cos(c + d*x)**3/(20*d) + x*(A*a**2/2 + 3*A*b**2/8 + 3*B*a*b/4) - (4*B*b**2 + 5*a*(2*A*b + B*a))*sin(c + d*x)**3/(15*d) + (4*B*b**2 + 5*a*(2*A*b + B*a))*sin(c + d*x)/(5*d) + (4*A*a**2 + 3*A*b**2 + 6*B*a*b)*sin(c + d*x)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_223():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*cos(c + d*x)
    F = B*(a + b*cos(c + d*x))**3*sin(c + d*x)/(4*b*d) + x*(A*a*b + B*a**2/2 + 3*B*b**2/8) + (8*A*a*b - 2*B*a**2 + 9*B*b**2)*sin(c + d*x)*cos(c + d*x)/(24*d) + (a + b*cos(c + d*x))**2*(4*A*b - B*a)*sin(c + d*x)/(12*b*d) + (4*A*a**2*b + 4*A*b**3 - B*a**3 + 8*B*a*b**2)*sin(c + d*x)/(6*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_224():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2
    F = B*(a + b*cos(c + d*x))**2*sin(c + d*x)/(3*d) + b*(3*A*b + 2*B*a)*sin(c + d*x)*cos(c + d*x)/(6*d) + x*(A*a**2 + A*b**2/2 + B*a*b) + (6*A*a*b + 2*B*a**2 + 2*B*b**2)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_225():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)
    F = A*a**2*atanh(sin(c + d*x))/d + B*b*(a + b*cos(c + d*x))*sin(c + d*x)/(2*d) + b*(2*A*b + 3*B*a)*sin(c + d*x)/(2*d) + x*(4*A*a*b + 2*B*a**2 + B*b**2)/2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_226():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)**2
    F = A*a**2*tan(c + d*x)/d + B*b**2*sin(c + d*x)/d + a*(2*A*b + B*a)*atanh(sin(c + d*x))/d + b*x*(A*b + 2*B*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_227():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)**3
    F = A*a**2*tan(c + d*x)*sec(c + d*x)/(2*d) + B*b**2*x + a*(2*A*b + B*a)*tan(c + d*x)/d + (A*a**2 + 2*A*b**2 + 4*B*a*b)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_228():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)**4
    F = A*a**2*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a*(2*A*b + B*a)*tan(c + d*x)*sec(c + d*x)/(2*d) + (2*A*a**2 + 3*A*b**2 + 6*B*a*b)*tan(c + d*x)/(3*d) + (2*A*a*b + B*a**2 + 2*B*b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_229():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)**5
    F = A*a**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a*(2*A*b + B*a)*tan(c + d*x)*sec(c + d*x)**2/(3*d) + (3*A*a**2 + 4*A*b**2 + 8*B*a*b)*tan(c + d*x)*sec(c + d*x)/(8*d) + (3*A*a**2 + 4*A*b**2 + 8*B*a*b)*atanh(sin(c + d*x))/(8*d) + (4*A*a*b + 2*B*a**2 + 3*B*b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_230():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*cos(c + d*x)**2
    F = B*b*(a + b*cos(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**3/(6*d) + b**2*(3*A*b + 4*B*a)*sin(c + d*x)*cos(c + d*x)**4/(15*d) + b*(18*A*a*b + 14*B*a**2 + 5*B*b**2)*sin(c + d*x)*cos(c + d*x)**3/(24*d) + x*(A*a**3/2 + 9*A*a*b**2/8 + 9*B*a**2*b/8 + 5*B*b**3/16) + (8*A*a**3 + 18*A*a*b**2 + 18*B*a**2*b + 5*B*b**3)*sin(c + d*x)*cos(c + d*x)/(16*d) - (15*A*a**2*b + 4*A*b**3 + 5*B*a**3 + 12*B*a*b**2)*sin(c + d*x)**3/(15*d) + (15*A*a**2*b + 4*A*b**3 + 5*B*a**3 + 12*B*a*b**2)*sin(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_231():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*cos(c + d*x)
    F = B*(a + b*cos(c + d*x))**4*sin(c + d*x)/(5*b*d) + x*(3*A*a**2*b/2 + 3*A*b**3/8 + B*a**3/2 + 9*B*a*b**2/8) + (30*A*a**2*b + 45*A*b**3 - 6*B*a**3 + 71*B*a*b**2)*sin(c + d*x)*cos(c + d*x)/(120*d) + (a + b*cos(c + d*x))**3*(5*A*b - B*a)*sin(c + d*x)/(20*b*d) + (a + b*cos(c + d*x))**2*(15*A*a*b - 3*B*a**2 + 16*B*b**2)*sin(c + d*x)/(60*b*d) + (15*A*a**3*b + 60*A*a*b**3 - 3*B*a**4 + 52*B*a**2*b**2 + 16*B*b**4)*sin(c + d*x)/(30*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_232():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3
    F = B*(a + b*cos(c + d*x))**3*sin(c + d*x)/(4*d) + b*(20*A*a*b + 6*B*a**2 + 9*B*b**2)*sin(c + d*x)*cos(c + d*x)/(24*d) + x*(A*a**3 + 3*A*a*b**2/2 + 3*B*a**2*b/2 + 3*B*b**3/8) + (a + b*cos(c + d*x))**2*(4*A*b + 3*B*a)*sin(c + d*x)/(12*d) + (16*A*a**2*b + 4*A*b**3 + 3*B*a**3 + 12*B*a*b**2)*sin(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_233():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)
    F = A*a**3*atanh(sin(c + d*x))/d + B*b*(a + b*cos(c + d*x))**2*sin(c + d*x)/(3*d) + b**2*(3*A*b + 5*B*a)*sin(c + d*x)*cos(c + d*x)/(6*d) + b*(9*A*a*b + 8*B*a**2 + 2*B*b**2)*sin(c + d*x)/(3*d) + x*(3*A*a**2*b + A*b**3/2 + B*a**3 + 3*B*a*b**2/2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_234():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**2
    F = A*a*(a + b*cos(c + d*x))**2*tan(c + d*x)/d + a**2*(3*A*b + B*a)*atanh(sin(c + d*x))/d - b**2*(2*A*a - B*b)*sin(c + d*x)*cos(c + d*x)/(2*d) + b*x*(6*A*a*b + 6*B*a**2 + B*b**2)/2 - b*(2*A*a**2 - A*b**2 - 3*B*a*b)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_235():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**3
    F = A*a*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)/(2*d) + a**2*(2*A*b + B*a)*tan(c + d*x)/d + a*(A*a**2 + 6*A*b**2 + 6*B*a*b)*atanh(sin(c + d*x))/(2*d) + b**2*x*(A*b + 3*B*a) - b**2*(A*a - 2*B*b)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_236():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**4
    F = A*a*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)**2/(3*d) + B*b**3*x + a**2*(5*A*b + 3*B*a)*tan(c + d*x)*sec(c + d*x)/(6*d) + a*(2*A*a**2 + 8*A*b**2 + 9*B*a*b)*tan(c + d*x)/(3*d) + (3*A*a**2*b + 2*A*b**3 + B*a**3 + 6*B*a*b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_237():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**5
    F = A*a*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)**3/(4*d) + a**2*(3*A*b + 2*B*a)*tan(c + d*x)*sec(c + d*x)**2/(6*d) + a*(3*A*a**2 + 10*A*b**2 + 12*B*a*b)*tan(c + d*x)*sec(c + d*x)/(8*d) + (3*A*a**3 + 12*A*a*b**2 + 12*B*a**2*b + 8*B*b**3)*atanh(sin(c + d*x))/(8*d) + (6*A*a**2*b + 3*A*b**3 + 2*B*a**3 + 9*B*a*b**2)*tan(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_238():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**6
    F = A*a*(a + b*cos(c + d*x))**2*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a**2*(7*A*b + 5*B*a)*tan(c + d*x)*sec(c + d*x)**3/(20*d) + a*(4*A*a**2 + 12*A*b**2 + 15*B*a*b)*tan(c + d*x)*sec(c + d*x)**2/(15*d) + (8*A*a**3 + 30*A*a*b**2 + 30*B*a**2*b + 15*B*b**3)*tan(c + d*x)/(15*d) + (9*A*a**2*b + 4*A*b**3 + 3*B*a**3 + 12*B*a*b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + (9*A*a**2*b + 4*A*b**3 + 3*B*a**3 + 12*B*a*b**2)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_239():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*cos(c + d*x)**2
    F = B*b*(a + b*cos(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**3/(7*d) + b**2*(49*A*a*b + 31*B*a**2 + 18*B*b**2)*sin(c + d*x)*cos(c + d*x)**4/(105*d) + b*(a + b*cos(c + d*x))**2*(7*A*b + 10*B*a)*sin(c + d*x)*cos(c + d*x)**3/(42*d) + b*(224*A*a**2*b + 35*A*b**3 + 104*B*a**3 + 140*B*a*b**2)*sin(c + d*x)*cos(c + d*x)**3/(168*d) + x*(A*a**4/2 + 9*A*a**2*b**2/4 + 5*A*b**4/16 + 3*B*a**3*b/2 + 5*B*a*b**3/4) + (8*A*a**4 + 36*A*a**2*b**2 + 5*A*b**4 + 24*B*a**3*b + 20*B*a*b**3)*sin(c + d*x)*cos(c + d*x)/(16*d) - (140*A*a**3*b + 112*A*a*b**3 + 35*B*a**4 + 168*B*a**2*b**2 + 24*B*b**4)*sin(c + d*x)**3/(105*d) + (140*A*a**3*b + 112*A*a*b**3 + 35*B*a**4 + 168*B*a**2*b**2 + 24*B*b**4)*sin(c + d*x)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_240():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*cos(c + d*x)
    F = B*(a + b*cos(c + d*x))**5*sin(c + d*x)/(6*b*d) + x*(2*A*a**3*b + 3*A*a*b**3/2 + B*a**4/2 + 9*B*a**2*b**2/4 + 5*B*b**4/16) + (48*A*a**3*b + 232*A*a*b**3 - 8*B*a**4 + 178*B*a**2*b**2 + 75*B*b**4)*sin(c + d*x)*cos(c + d*x)/(240*d) + (a + b*cos(c + d*x))**4*(6*A*b - B*a)*sin(c + d*x)/(30*b*d) + (a + b*cos(c + d*x))**3*(24*A*a*b - 4*B*a**2 + 25*B*b**2)*sin(c + d*x)/(120*b*d) + (a + b*cos(c + d*x))**2*(24*A*a**2*b + 32*A*b**3 - 4*B*a**3 + 53*B*a*b**2)*sin(c + d*x)/(120*b*d) + (24*A*a**4*b + 224*A*a**2*b**3 + 32*A*b**5 - 4*B*a**5 + 121*B*a**3*b**2 + 128*B*a*b**4)*sin(c + d*x)/(60*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_241():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4
    F = B*(a + b*cos(c + d*x))**4*sin(c + d*x)/(5*d) + b*(130*A*a**2*b + 45*A*b**3 + 24*B*a**3 + 116*B*a*b**2)*sin(c + d*x)*cos(c + d*x)/(120*d) + x*(A*a**4 + 3*A*a**2*b**2 + 3*A*b**4/8 + 2*B*a**3*b + 3*B*a*b**3/2) + (a + b*cos(c + d*x))**3*(5*A*b + 4*B*a)*sin(c + d*x)/(20*d) + (a + b*cos(c + d*x))**2*(35*A*a*b + 12*B*a**2 + 16*B*b**2)*sin(c + d*x)/(60*d) + (95*A*a**3*b + 80*A*a*b**3 + 12*B*a**4 + 112*B*a**2*b**2 + 16*B*b**4)*sin(c + d*x)/(30*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_242():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*sec(c + d*x)
    F = A*a**4*atanh(sin(c + d*x))/d + B*b*(a + b*cos(c + d*x))**3*sin(c + d*x)/(4*d) + b**2*(32*A*a*b + 26*B*a**2 + 9*B*b**2)*sin(c + d*x)*cos(c + d*x)/(24*d) + b*(a + b*cos(c + d*x))**2*(4*A*b + 7*B*a)*sin(c + d*x)/(12*d) + b*(34*A*a**2*b + 4*A*b**3 + 19*B*a**3 + 16*B*a*b**2)*sin(c + d*x)/(6*d) + x*(4*A*a**3*b + 2*A*a*b**3 + B*a**4 + 3*B*a**2*b**2 + 3*B*b**4/8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_243():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*sec(c + d*x)**2
    F = A*a*(a + b*cos(c + d*x))**3*tan(c + d*x)/d + a**3*(4*A*b + B*a)*atanh(sin(c + d*x))/d - b**2*(6*A*a**2 - 3*A*b**2 - 8*B*a*b)*sin(c + d*x)*cos(c + d*x)/(6*d) + b*x*(12*A*a**2*b + A*b**3 + 8*B*a**3 + 4*B*a*b**2)/2 - b*(a + b*cos(c + d*x))**2*(3*A*a - B*b)*sin(c + d*x)/(3*d) - b*(6*A*a**3 - 12*A*a*b**2 - 17*B*a**2*b - 2*B*b**3)*sin(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_244():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*sec(c + d*x)**3
    F = A*a*(a + b*cos(c + d*x))**3*tan(c + d*x)*sec(c + d*x)/(2*d) + a**2*(A*a**2 + 12*A*b**2 + 8*B*a*b)*atanh(sin(c + d*x))/(2*d) + a*(a + b*cos(c + d*x))**2*(5*A*b + 2*B*a)*tan(c + d*x)/(2*d) + b**2*x*(8*A*a*b + 12*B*a**2 + B*b**2)/2 - b**2*(6*A*a*b + 2*B*a**2 - B*b**2)*sin(c + d*x)*cos(c + d*x)/(2*d) - b*(13*A*a**2*b - 2*A*b**3 + 4*B*a**3 - 8*B*a*b**2)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_245():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*sec(c + d*x)**4
    F = A*a*(a + b*cos(c + d*x))**3*tan(c + d*x)*sec(c + d*x)**2/(3*d) + a**2*(2*A*a**2 + 9*A*b**2 + 9*B*a*b)*tan(c + d*x)/(3*d) + a*(a + b*cos(c + d*x))**2*(2*A*b + B*a)*tan(c + d*x)*sec(c + d*x)/(2*d) + a*(4*A*a**2*b + 8*A*b**3 + B*a**3 + 12*B*a*b**2)*atanh(sin(c + d*x))/(2*d) + b**3*x*(A*b + 4*B*a) - b**2*(8*A*a*b + 3*B*a**2 - 6*B*b**2)*sin(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_246():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*sec(c + d*x)**5
    F = A*a*(a + b*cos(c + d*x))**3*tan(c + d*x)*sec(c + d*x)**3/(4*d) + B*b**4*x + a**2*(9*A*a**2 + 26*A*b**2 + 32*B*a*b)*tan(c + d*x)*sec(c + d*x)/(24*d) + a*(a + b*cos(c + d*x))**2*(7*A*b + 4*B*a)*tan(c + d*x)*sec(c + d*x)**2/(12*d) + a*(16*A*a**2*b + 19*A*b**3 + 4*B*a**3 + 34*B*a*b**2)*tan(c + d*x)/(6*d) + (3*A*a**4 + 24*A*a**2*b**2 + 8*A*b**4 + 16*B*a**3*b + 32*B*a*b**3)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_247():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*sec(c + d*x)**6
    F = A*a*(a + b*cos(c + d*x))**3*tan(c + d*x)*sec(c + d*x)**4/(5*d) + a**2*(8*A*a**2 + 18*A*b**2 + 25*B*a*b)*tan(c + d*x)*sec(c + d*x)**2/(30*d) + a*(a + b*cos(c + d*x))**2*(8*A*b + 5*B*a)*tan(c + d*x)*sec(c + d*x)**3/(20*d) + a*(60*A*a**2*b + 56*A*b**3 + 15*B*a**3 + 110*B*a*b**2)*tan(c + d*x)*sec(c + d*x)/(40*d) + (8*A*a**4 + 60*A*a**2*b**2 + 15*A*b**4 + 40*B*a**3*b + 60*B*a*b**3)*tan(c + d*x)/(15*d) + (12*A*a**3*b + 16*A*a*b**3 + 3*B*a**4 + 24*B*a**2*b**2 + 8*B*b**4)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_248():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**4*sec(c + d*x)**7
    F = A*a*(a + b*cos(c + d*x))**3*tan(c + d*x)*sec(c + d*x)**5/(6*d) + a**2*(25*A*a**2 + 48*A*b**2 + 72*B*a*b)*tan(c + d*x)*sec(c + d*x)**3/(120*d) + a*(a + b*cos(c + d*x))**2*(3*A*b + 2*B*a)*tan(c + d*x)*sec(c + d*x)**4/(10*d) + a*(16*A*a**2*b + 13*A*b**3 + 4*B*a**3 + 27*B*a*b**2)*tan(c + d*x)*sec(c + d*x)**2/(15*d) + (5*A*a**4 + 36*A*a**2*b**2 + 8*A*b**4 + 24*B*a**3*b + 32*B*a*b**3)*tan(c + d*x)*sec(c + d*x)/(16*d) + (5*A*a**4 + 36*A*a**2*b**2 + 8*A*b**4 + 24*B*a**3*b + 32*B*a*b**3)*atanh(sin(c + d*x))/(16*d) + (32*A*a**3*b + 40*A*a*b**3 + 8*B*a**4 + 60*B*a**2*b**2 + 15*B*b**4)*tan(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_249():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))
    F = B*sin(c + d*x)*cos(c + d*x)**2/(3*b*d) - 2*a**3*(A*b - B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*sqrt(a - b)*sqrt(a + b)) + (A*b - B*a)*sin(c + d*x)*cos(c + d*x)/(2*b**2*d) - (3*A*a*b - 3*B*a**2 - 2*B*b**2)*sin(c + d*x)/(3*b**3*d) + x*(2*a**2 + b**2)*(A*b - B*a)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_250():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))
    F = B*sin(c + d*x)*cos(c + d*x)/(2*b*d) + 2*a**2*(A*b - B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*sqrt(a - b)*sqrt(a + b)) + (A*b - B*a)*sin(c + d*x)/(b**2*d) - x*(2*A*a*b - 2*B*a**2 - B*b**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_251():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))
    F = B*sin(c + d*x)/(b*d) - 2*a*(A*b - B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*sqrt(a - b)*sqrt(a + b)) + x*(A*b - B*a)/b**2
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_252():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))
    F = B*x/b + (2*A*b - 2*B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b*d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_253():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))
    F = A*atanh(sin(c + d*x))/(a*d) - (2*A*b - 2*B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_254():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))
    F = A*tan(c + d*x)/(a*d) + 2*b*(A*b - B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*sqrt(a - b)*sqrt(a + b)) - (A*b - B*a)*atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_255():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))
    F = A*tan(c + d*x)*sec(c + d*x)/(2*a*d) - (A*b - B*a)*tan(c + d*x)/(a**2*d) - 2*b**2*(A*b - B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*sqrt(a - b)*sqrt(a + b)) + (A*a**2 + 2*A*b**2 - 2*B*a*b)*atanh(sin(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_256():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**4/(a + b*cos(c + d*x))
    F = A*tan(c + d*x)*sec(c + d*x)**2/(3*a*d) - (A*b - B*a)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) + (2*A*a**2 + 3*A*b**2 - 3*B*a*b)*tan(c + d*x)/(3*a**3*d) + 2*b**3*(A*b - B*a)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*sqrt(a - b)*sqrt(a + b)) - (a**2 + 2*b**2)*(A*b - B*a)*atanh(sin(c + d*x))/(2*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_257():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))**2
    F = 2*a**2*(2*A*a**2*b - 3*A*b**3 - 3*B*a**3 + 4*B*a*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**2/(b*d*(a + b*cos(c + d*x))*(a**2 - b**2)) - (2*A*a*b - 3*B*a**2 + B*b**2)*sin(c + d*x)*cos(c + d*x)/(2*b**2*d*(a**2 - b**2)) + (2*A*a**2*b - A*b**3 - 3*B*a**3 + 2*B*a*b**2)*sin(c + d*x)/(b**3*d*(a**2 - b**2)) - x*(4*A*a*b - 6*B*a**2 - B*b**2)/(2*b**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_258():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))**2
    F = B*sin(c + d*x)/(b**2*d) - a**2*(A*b - B*a)*sin(c + d*x)/(b**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)) - 2*a*(A*a**2*b - 2*A*b**3 - 2*B*a**3 + 3*B*a*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + x*(A*b - 2*B*a)/b**3
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_259():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))**2
    F = B*x/b**2 + a*(A*b - B*a)*sin(c + d*x)/(b*d*(a + b*cos(c + d*x))*(a**2 - b**2)) - (2*A*b**3 + 2*B*a**3 - 4*B*a*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_260():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))**2
    F = -(A*b - B*a)*sin(c + d*x)/(d*(a + b*cos(c + d*x))*(a**2 - b**2)) + (2*A*a - 2*B*b)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_261():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**2
    F = A*atanh(sin(c + d*x))/(a**2*d) + b*(A*b - B*a)*sin(c + d*x)/(a*d*(a + b*cos(c + d*x))*(a**2 - b**2)) - (4*A*a**2*b - 2*A*b**3 - 2*B*a**3)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_262():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))**2
    F = b*(A*b - B*a)*tan(c + d*x)/(a*d*(a + b*cos(c + d*x))*(a**2 - b**2)) + (A*a**2 - 2*A*b**2 + B*a*b)*tan(c + d*x)/(a**2*d*(a**2 - b**2)) + 2*b*(3*A*a**2*b - 2*A*b**3 - 2*B*a**3 + B*a*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) - (2*A*b - B*a)*atanh(sin(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_263():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))**2
    F = b*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)/(a*d*(a + b*cos(c + d*x))*(a**2 - b**2)) + (A*a**2 - 3*A*b**2 + 2*B*a*b)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d*(a**2 - b**2)) - (2*A*a**2*b - 3*A*b**3 - B*a**3 + 2*B*a*b**2)*tan(c + d*x)/(a**3*d*(a**2 - b**2)) - 2*b**2*(4*A*a**2*b - 3*A*b**3 - 3*B*a**3 + 2*B*a*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(3)/2)*(a + b)**(sympy.S(3)/2)) + (A*a**2 + 6*A*b**2 - 4*B*a*b)*atanh(sin(c + d*x))/(2*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_264():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a + b*cos(c + d*x))**3
    F = a**2*(6*A*a**4*b - 15*A*a**2*b**3 + 12*A*b**5 - 12*B*a**5 + 29*B*a**3*b**2 - 20*B*a*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**3/(2*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + a*(2*A*a**2*b - 5*A*b**3 - 4*B*a**3 + 7*B*a*b**2)*sin(c + d*x)*cos(c + d*x)**2/(2*b**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - (3*A*a**3*b - 6*A*a*b**3 - 6*B*a**4 + 10*B*a**2*b**2 - B*b**4)*sin(c + d*x)*cos(c + d*x)/(2*b**3*d*(a**2 - b**2)**2) + (6*A*a**4*b - 11*A*a**2*b**3 + 2*A*b**5 - 12*B*a**5 + 21*B*a**3*b**2 - 6*B*a*b**4)*sin(c + d*x)/(2*b**4*d*(a**2 - b**2)**2) - x*(6*A*a*b - 12*B*a**2 - B*b**2)/(2*b**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_265():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))**3
    F = -a**2*(A*a**2*b - 4*A*b**3 - 3*B*a**3 + 6*B*a*b**2)*sin(c + d*x)/(2*b**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) + a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**2/(2*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) - a*(2*A*a**4*b - 5*A*a**2*b**3 + 6*A*b**5 - 6*B*a**5 + 15*B*a**3*b**2 - 12*B*a*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - (A*a*b - 3*B*a**2 + 2*B*b**2)*sin(c + d*x)/(2*b**3*d*(a**2 - b**2)) + x*(A*b - 3*B*a)/b**4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_266():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))**3
    F = B*x/b**3 - a**2*(A*b - B*a)*sin(c + d*x)/(2*b**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + a*(A*a**2*b - 4*A*b**3 - 3*B*a**3 + 6*B*a*b**2)*sin(c + d*x)/(2*b**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) + (A*a**2*b**3 + 2*A*b**5 - 2*B*a**5 + 5*B*a**3*b**2 - 6*B*a*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_267():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))**3
    F = a*(A*b - B*a)*sin(c + d*x)/(2*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) - (3*A*a*b - B*a**2 - 2*B*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + (A*a**2*b + 2*A*b**3 + B*a**3 - 4*B*a*b**2)*sin(c + d*x)/(2*b*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_268():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))**3
    F = -(3*A*a*b - B*a**2 - 2*B*b**2)*sin(c + d*x)/(2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - (A*b - B*a)*sin(c + d*x)/(d*(a + b*cos(c + d*x))**2*(2*a**2 - 2*b**2)) + (2*A*a**2 + A*b**2 - 3*B*a*b)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_269():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**3
    F = A*atanh(sin(c + d*x))/(a**3*d) + b*(A*b - B*a)*sin(c + d*x)/(2*a*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + b*(5*A*a**2*b - 2*A*b**3 - 3*B*a**3)*sin(c + d*x)/(2*a**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) - (6*A*a**4*b - 5*A*a**2*b**3 + 2*A*b**5 - 2*B*a**5 - B*a**3*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_270():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))**3
    F = b*(A*b - B*a)*tan(c + d*x)/(2*a*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + b*(6*A*a**2*b - 3*A*b**3 - 4*B*a**3 + B*a*b**2)*tan(c + d*x)/(2*a**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) + (2*A*a**4 - 11*A*a**2*b**2 + 6*A*b**4 + 5*B*a**3*b - 2*B*a*b**3)*tan(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) + b*(12*A*a**4*b - 15*A*a**2*b**3 + 6*A*b**5 - 6*B*a**5 + 5*B*a**3*b**2 - 2*B*a*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) - (3*A*b - B*a)*atanh(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_271():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))**3
    F = b*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)/(2*a*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)) + b*(7*A*a**2*b - 4*A*b**3 - 5*B*a**3 + 2*B*a*b**2)*tan(c + d*x)*sec(c + d*x)/(2*a**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**2) + (A*a**4 - 10*A*a**2*b**2 + 6*A*b**4 + 6*B*a**3*b - 3*B*a*b**3)*tan(c + d*x)*sec(c + d*x)/(2*a**3*d*(a**2 - b**2)**2) - (6*A*a**4*b - 21*A*a**2*b**3 + 12*A*b**5 - 2*B*a**5 + 11*B*a**3*b**2 - 6*B*a*b**4)*tan(c + d*x)/(2*a**4*d*(a**2 - b**2)**2) - b**2*(20*A*a**4*b - 29*A*a**2*b**3 + 12*A*b**5 - 12*B*a**5 + 15*B*a**3*b**2 - 6*B*a*b**4)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(5)/2)*(a + b)**(sympy.S(5)/2)) + (A*a**2 + 12*A*b**2 - 6*B*a*b)*atanh(sin(c + d*x))/(2*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_272():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a + b*cos(c + d*x))**4
    F = -a**2*(A*a**4*b - 2*A*a**2*b**3 + 6*A*b**5 - 4*B*a**5 + 11*B*a**3*b**2 - 12*B*a*b**4)*sin(c + d*x)/(2*b**4*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) + a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**3/(3*b*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + a*(A*a**2*b - 6*A*b**3 - 4*B*a**3 + 9*B*a*b**2)*sin(c + d*x)*cos(c + d*x)**2/(6*b**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) - a*(2*A*a**6*b - 7*A*a**4*b**3 + 8*A*a**2*b**5 - 8*A*b**7 - 8*B*a**7 + 28*B*a**5*b**2 - 35*B*a**3*b**4 + 20*B*a*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - (3*A*a**3*b - 8*A*a*b**3 - 12*B*a**4 + 23*B*a**2*b**2 - 6*B*b**4)*sin(c + d*x)/(6*b**4*d*(a**2 - b**2)**2) + x*(A*b - 4*B*a)/b**5
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_273():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))**4
    F = B*x/b**4 + a**2*(5*A*b**3 + 3*B*a**3 - 8*B*a*b**2)*sin(c + d*x)/(6*b**3*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**2/(3*b*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) - a*(A*a**2*b**3 - 16*A*b**5 + 9*B*a**5 - 28*B*a**3*b**2 + 34*B*a*b**4)*sin(c + d*x)/(6*b**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - (3*A*a**2*b**5 + 2*A*b**7 + 2*B*a**7 - 7*B*a**5*b**2 + 8*B*a**3*b**4 - 8*B*a*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_274():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))**4
    F = -a**2*(A*b - B*a)*sin(c + d*x)/(3*b**2*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + a*(A*a**2*b - 6*A*b**3 - 4*B*a**3 + 9*B*a*b**2)*sin(c + d*x)/(6*b**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + (A*a**3 + 4*A*a*b**2 - 3*B*a**2*b - 2*B*b**3)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (A*a**4*b - 10*A*a**2*b**3 - 6*A*b**5 + 2*B*a**5 - 5*B*a**3*b**2 + 18*B*a*b**4)*sin(c + d*x)/(6*b**2*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_275():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))**4
    F = a*(A*b - B*a)*sin(c + d*x)/(3*b*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) - (4*A*a**2*b + A*b**3 - B*a**3 - 4*B*a*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (2*A*a**3*b + 13*A*a*b**3 + B*a**4 - 10*B*a**2*b**2 - 6*B*b**4)*sin(c + d*x)/(6*b*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) + (2*A*a**2*b + 3*A*b**3 + B*a**3 - 6*B*a*b**2)*sin(c + d*x)/(6*b*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_276():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))**4
    F = -(11*A*a**2*b + 4*A*b**3 - 2*B*a**3 - 13*B*a*b**2)*sin(c + d*x)/(6*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - (5*A*a*b - 2*B*a**2 - 3*B*b**2)*sin(c + d*x)/(6*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) - (A*b - B*a)*sin(c + d*x)/(d*(a + b*cos(c + d*x))**3*(3*a**2 - 3*b**2)) + (2*A*a**3 + 3*A*a*b**2 - 4*B*a**2*b - B*b**3)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_277():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**4
    F = A*atanh(sin(c + d*x))/(a**4*d) + b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + b*(8*A*a**2*b - 3*A*b**3 - 5*B*a**3)*sin(c + d*x)/(6*a**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + b*(26*A*a**4*b - 17*A*a**2*b**3 + 6*A*b**5 - 11*B*a**5 - 4*B*a**3*b**2)*sin(c + d*x)/(6*a**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) - (8*A*a**6*b - 8*A*a**4*b**3 + 7*A*a**2*b**5 - 2*A*b**7 - 2*B*a**7 - 3*B*a**5*b**2)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**4*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_278():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))**4
    F = b*(A*b - B*a)*tan(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + b*(9*A*a**2*b - 4*A*b**3 - 6*B*a**3 + B*a*b**2)*tan(c + d*x)/(6*a**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + b*(12*A*a**4*b - 11*A*a**2*b**3 + 4*A*b**5 - 6*B*a**5 + 2*B*a**3*b**2 - B*a*b**4)*tan(c + d*x)/(2*a**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) + (6*A*a**6 - 65*A*a**4*b**2 + 68*A*a**2*b**4 - 24*A*b**6 + 26*B*a**5*b - 17*B*a**3*b**3 + 6*B*a*b**5)*tan(c + d*x)/(6*a**4*d*(a**2 - b**2)**3) + b*(20*A*a**6*b - 35*A*a**4*b**3 + 28*A*a**2*b**5 - 8*A*b**7 - 8*B*a**7 + 8*B*a**5*b**2 - 7*B*a**3*b**4 + 2*B*a*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**5*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) - (4*A*b - B*a)*atanh(sin(c + d*x))/(a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_279():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))**4
    F = b*(A*b - B*a)*tan(c + d*x)*sec(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**3*(a**2 - b**2)) + b*(10*A*a**2*b - 5*A*b**3 - 7*B*a**3 + 2*B*a*b**2)*tan(c + d*x)*sec(c + d*x)/(6*a**2*d*(a + b*cos(c + d*x))**2*(a**2 - b**2)**2) + b*(48*A*a**4*b - 53*A*a**2*b**3 + 20*A*b**5 - 27*B*a**5 + 20*B*a**3*b**2 - 8*B*a*b**4)*tan(c + d*x)*sec(c + d*x)/(6*a**3*d*(a + b*cos(c + d*x))*(a**2 - b**2)**3) + (A*a**6 - 23*A*a**4*b**2 + 27*A*a**2*b**4 - 10*A*b**6 + 12*B*a**5*b - 11*B*a**3*b**3 + 4*B*a*b**5)*tan(c + d*x)*sec(c + d*x)/(2*a**4*d*(a**2 - b**2)**3) - (24*A*a**6*b - 146*A*a**4*b**3 + 167*A*a**2*b**5 - 60*A*b**7 - 6*B*a**7 + 65*B*a**5*b**2 - 68*B*a**3*b**4 + 24*B*a*b**6)*tan(c + d*x)/(6*a**5*d*(a**2 - b**2)**3) - b**2*(40*A*a**6*b - 84*A*a**4*b**3 + 69*A*a**2*b**5 - 20*A*b**7 - 20*B*a**7 + 35*B*a**5*b**2 - 28*B*a**3*b**4 + 8*B*a*b**6)*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**6*d*(a - b)**(sympy.S(7)/2)*(a + b)**(sympy.S(7)/2)) + (A*a**2 + 20*A*b**2 - 8*B*a*b)*atanh(sin(c + d*x))/(2*a**6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_280():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))
    F = -B*sin(c + d*x)**3/(3*d) + B*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_281():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))
    F = B*x/2 + B*sin(c + d*x)*cos(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_282():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))
    F = B*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_283():
    f = (B*a + B*b*cos(c + d*x))/(a + b*cos(c + d*x))
    F = B*x
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_284():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))
    F = B*atanh(sin(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_285():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))
    F = B*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_286():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))
    F = B*tan(c + d*x)*sec(c + d*x)/(2*d) + B*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_287():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**4/(a + b*cos(c + d*x))
    F = B*tan(c + d*x)**3/(3*d) + B*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_288():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))**2
    F = -2*B*a**3*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**3*d*sqrt(a - b)*sqrt(a + b)) - B*a*sin(c + d*x)/(b**2*d) + B*sin(c + d*x)*cos(c + d*x)/(2*b*d) + B*x*(2*a**2 + b**2)/(2*b**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_289():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))**2
    F = 2*B*a**2*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b**2*d*sqrt(a - b)*sqrt(a + b)) - B*a*x/b**2 + B*sin(c + d*x)/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_290():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))**2
    F = -2*B*a*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(b*d*sqrt(a - b)*sqrt(a + b)) + B*x/b
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_291():
    f = (B*a + B*b*cos(c + d*x))/(a + b*cos(c + d*x))**2
    F = 2*B*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(d*sqrt(a - b)*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_292():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**2
    F = -2*B*b*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a*d*sqrt(a - b)*sqrt(a + b)) + B*atanh(sin(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_293():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))**2
    F = B*tan(c + d*x)/(a*d) + 2*B*b**2*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**2*d*sqrt(a - b)*sqrt(a + b)) - B*b*atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_294():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))**2
    F = B*tan(c + d*x)*sec(c + d*x)/(2*a*d) - B*b*tan(c + d*x)/(a**2*d) - 2*B*b**3*atan(sqrt(a - b)*tan(c/2 + d*x/2)/sqrt(a + b))/(a**3*d*sqrt(a - b)*sqrt(a + b)) + B*(a**2 + 2*b**2)*atanh(sin(c + d*x))/(2*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_295():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*cos(c + d*x)**3
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)**2/(9*b*d) + (a + b*cos(c + d*x))**(sympy.S(3)/2)*(6*A*b - 4*B*a)*sin(c + d*x)*cos(c + d*x)/(21*b**2*d) - (a + b*cos(c + d*x))**(sympy.S(3)/2)*(72*A*a*b - 48*B*a**2 - 98*B*b**2)*sin(c + d*x)/(315*b**3*d) + sqrt(a + b*cos(c + d*x))*(48*A*a**2*b + 150*A*b**3 - 32*B*a**3 - 72*B*a*b**2)*sin(c + d*x)/(315*b**3*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(24*A*a**2*b + 75*A*b**3 - 16*B*a**3 - 36*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(315*b**4*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(48*A*a**3*b + 114*A*a*b**3 - 32*B*a**4 - 48*B*a**2*b**2 + 294*B*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(315*b**4*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_296():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*cos(c + d*x)**2
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)/(7*b*d) + (a + b*cos(c + d*x))**(sympy.S(3)/2)*(14*A*b - 8*B*a)*sin(c + d*x)/(35*b**2*d) - sqrt(a + b*cos(c + d*x))*(28*A*a*b - 16*B*a**2 - 50*B*b**2)*sin(c + d*x)/(105*b**2*d) + sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(14*A*a*b - 8*B*a**2 - 25*B*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(105*b**3*d*sqrt(a + b*cos(c + d*x))) - sqrt(a + b*cos(c + d*x))*(28*A*a**2*b - 126*A*b**3 - 16*B*a**3 - 38*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(105*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_297():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*cos(c + d*x)
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*b*d) + sqrt(a + b*cos(c + d*x))*(10*A*b - 4*B*a)*sin(c + d*x)/(15*b*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(5*A*b - 2*B*a)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b**2*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(10*A*a*b - 4*B*a**2 + 18*B*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_298():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))
    F = 2*B*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*d) - B*sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(6*A*b + 2*B*a)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_299():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)
    F = ((Integer(2) * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('A') * Symbol('b') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_300():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)**2
    F = (Integer(-1) * ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Symbol('a') * Symbol('A')) + (Integer(2) * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_301():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)**3
    F = (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(4) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_302():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)**4
    F = (Integer(-1) * ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(6) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(18) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(2) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(6) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Symbol('A') * Symbol('b')) + (Integer(6) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_303():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)*cos(c + d*x)/(9*b*d) + (a + b*cos(c + d*x))**(sympy.S(5)/2)*(18*A*b - 8*B*a)*sin(c + d*x)/(63*b**2*d) - (a + b*cos(c + d*x))**(sympy.S(3)/2)*(36*A*a*b - 16*B*a**2 - 98*B*b**2)*sin(c + d*x)/(315*b**2*d) - sqrt(a + b*cos(c + d*x))*(36*A*a**2*b - 150*A*b**3 - 16*B*a**3 - 78*B*a*b**2)*sin(c + d*x)/(315*b**2*d) + sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(18*A*a**2*b - 75*A*b**3 - 8*B*a**3 - 39*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(315*b**3*d*sqrt(a + b*cos(c + d*x))) - sqrt(a + b*cos(c + d*x))*(36*A*a**3*b - 492*A*a*b**3 - 16*B*a**4 - 66*B*a**2*b**2 - 294*B*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(315*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_304():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*b*d) + (a + b*cos(c + d*x))**(sympy.S(3)/2)*(14*A*b - 4*B*a)*sin(c + d*x)/(35*b*d) + sqrt(a + b*cos(c + d*x))*(42*A*a*b - 12*B*a**2 + 50*B*b**2)*sin(c + d*x)/(105*b*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(21*A*a*b - 6*B*a**2 + 25*B*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(105*b**2*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(42*A*a**2*b + 126*A*b**3 - 12*B*a**3 + 164*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(105*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_305():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(5*d) + sqrt(a + b*cos(c + d*x))*(10*A*b + 6*B*a)*sin(c + d*x)/(15*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(5*A*b + 3*B*a)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(40*A*a*b + 6*B*a**2 + 18*B*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_306():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)
    F = ((Integer(2) * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_307():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2
    F = (Integer(-1) * ((((Symbol('a') * Symbol('A')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_308():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**3
    F = (Integer(-1) * ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(7) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(8) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(12) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_309():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**4
    F = (Integer(-1) * ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(30) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(17) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(42) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(12) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(30) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((((Integer(7) * Symbol('A') * Symbol('b')) + (Integer(6) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_310():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)*cos(c + d*x)/(11*b*d) + (a + b*cos(c + d*x))**(sympy.S(7)/2)*(22*A*b - 8*B*a)*sin(c + d*x)/(99*b**2*d) - (a + b*cos(c + d*x))**(sympy.S(5)/2)*(44*A*a*b - 16*B*a**2 - 162*B*b**2)*sin(c + d*x)/(693*b**2*d) - (a + b*cos(c + d*x))**(sympy.S(3)/2)*(220*A*a**2*b - 1078*A*b**3 - 80*B*a**3 - 670*B*a*b**2)*sin(c + d*x)/(3465*b**2*d) - sqrt(a + b*cos(c + d*x))*(220*A*a**3*b - 2508*A*a*b**3 - 80*B*a**4 - 570*B*a**2*b**2 - 1350*B*b**4)*sin(c + d*x)/(3465*b**2*d) + sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(110*A*a**3*b - 1254*A*a*b**3 - 40*B*a**4 - 285*B*a**2*b**2 - 675*B*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3465*b**3*d*sqrt(a + b*cos(c + d*x))) - sqrt(a + b*cos(c + d*x))*(220*A*a**4*b - 6138*A*a**2*b**3 - 3234*A*b**5 - 80*B*a**5 - 510*B*a**3*b**2 - 7410*B*a*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3465*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_311():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(9*b*d) + (a + b*cos(c + d*x))**(sympy.S(5)/2)*(18*A*b - 4*B*a)*sin(c + d*x)/(63*b*d) + (a + b*cos(c + d*x))**(sympy.S(3)/2)*(90*A*a*b - 20*B*a**2 + 98*B*b**2)*sin(c + d*x)/(315*b*d) + sqrt(a + b*cos(c + d*x))*(90*A*a**2*b + 150*A*b**3 - 20*B*a**3 + 228*B*a*b**2)*sin(c + d*x)/(315*b*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(45*A*a**2*b + 75*A*b**3 - 10*B*a**3 + 114*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(315*b**2*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(90*A*a**3*b + 870*A*a*b**3 - 20*B*a**4 + 558*B*a**2*b**2 + 294*B*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(315*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_312():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*B*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(7*d) + (a + b*cos(c + d*x))**(sympy.S(3)/2)*(14*A*b + 10*B*a)*sin(c + d*x)/(35*d) + sqrt(a + b*cos(c + d*x))*(112*A*a*b + 30*B*a**2 + 50*B*b**2)*sin(c + d*x)/(105*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*a**2 - 2*b**2)*(56*A*a*b + 15*B*a**2 + 25*B*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(105*b*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(322*A*a**2*b + 126*A*b**3 + 30*B*a**3 + 290*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(105*b*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_313():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)
    F = ((Integer(2) * ((Integer(35) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(23) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(9) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(15) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * ((Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(8) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(15) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_314():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(6) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('A')) + (Integer(12) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B')) + (Integer(2) * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (((Symbol('a'))**(Integer(2)) * ((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(3) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Integer(2) * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_315():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**3
    F = (Integer(-1) * ((((Integer(9) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(16) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(20) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('a') * ((Integer(7) * Symbol('A') * Symbol('b')) + (Integer(4) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_316():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**4
    F = (Integer(-1) * ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(33) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(54) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(16) * (Symbol('a'))**(Integer(3)) * Symbol('A')) + (Integer(59) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(66) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B')) + (Integer(48) * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(20) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(30) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(8) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(16) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(33) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(54) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_317():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**5
    F = (Integer(-1) * ((((Integer(284) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(128) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(264) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(192) * Symbol('a') * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(356) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(133) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(128) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(472) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(192) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(48) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(120) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(4)))) + (Integer(160) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(40) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(64) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((((Integer(284) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(128) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(264) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((((Integer(36) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(59) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(104) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(96) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * ((Integer(11) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(2)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**(Integer(3)) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_318():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/sqrt(a + b*cos(c + d*x))
    F = 2*B*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**2/(7*b*d) + sqrt(a + b*cos(c + d*x))*(14*A*b - 12*B*a)*sin(c + d*x)*cos(c + d*x)/(35*b**2*d) - sqrt(a + b*cos(c + d*x))*(56*A*a*b - 48*B*a**2 - 50*B*b**2)*sin(c + d*x)/(105*b**3*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(112*A*a**3*b + 98*A*a*b**3 - 96*B*a**4 - 64*B*a**2*b**2 - 50*B*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(105*b**4*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(112*A*a**2*b + 126*A*b**3 - 96*B*a**3 - 88*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(105*b**4*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_319():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/sqrt(a + b*cos(c + d*x))
    F = 2*B*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)/(5*b*d) + sqrt(a + b*cos(c + d*x))*(10*A*b - 8*B*a)*sin(c + d*x)/(15*b**2*d) + sqrt((a + b*cos(c + d*x))/(a + b))*(20*A*a**2*b + 10*A*b**3 - 16*B*a**3 - 14*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b**3*d*sqrt(a + b*cos(c + d*x))) - sqrt(a + b*cos(c + d*x))*(20*A*a*b - 16*B*a**2 - 18*B*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_320():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/sqrt(a + b*cos(c + d*x))
    F = 2*B*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*b*d) - sqrt((a + b*cos(c + d*x))/(a + b))*(6*A*a*b - 4*B*a**2 - 2*B*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(6*A*b - 4*B*a)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_321():
    f = (A + B*cos(c + d*x))/sqrt(a + b*cos(c + d*x))
    F = 2*B*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt((a + b*cos(c + d*x))/(a + b))) + sqrt((a + b*cos(c + d*x))/(a + b))*(2*A*b - 2*B*a)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_322():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/sqrt(a + b*cos(c + d*x))
    F = ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_323():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/sqrt(a + b*cos(c + d*x))
    F = (Integer(-1) * ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_324():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/sqrt(a + b*cos(c + d*x))
    F = ((((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_325():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**2/(b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - sqrt(a + b*cos(c + d*x))*(10*A*a*b - 12*B*a**2 + 2*B*b**2)*sin(c + d*x)*cos(c + d*x)/(5*b**2*d*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(40*A*a**2*b - 10*A*b**3 - 48*B*a**3 + 18*B*a*b**2)*sin(c + d*x)/(15*b**3*d*(a**2 - b**2)) + sqrt((a + b*cos(c + d*x))/(a + b))*(80*A*a**2*b + 10*A*b**3 - 96*B*a**3 - 24*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b**4*d*sqrt(a + b*cos(c + d*x))) - sqrt(a + b*cos(c + d*x))*(80*A*a**3*b - 50*A*a*b**3 - 96*B*a**4 + 48*B*a**2*b**2 + 18*B*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b**4*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_326():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*b**2*d) - 2*a**2*(A*b - B*a)*sin(c + d*x)/(b**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - sqrt((a + b*cos(c + d*x))/(a + b))*(12*A*a*b - 16*B*a**2 - 2*B*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt(a + b*cos(c + d*x))) + sqrt(a + b*cos(c + d*x))*(12*A*a**2*b - 6*A*b**3 - 16*B*a**3 + 10*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_327():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*a*(A*b - B*a)*sin(c + d*x)/(b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt((a + b*cos(c + d*x))/(a + b))*(2*A*b - 4*B*a)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(b**2*d*sqrt(a + b*cos(c + d*x))) - sqrt(a + b*cos(c + d*x))*(2*A*a*b - 4*B*a**2 + 2*B*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(b**2*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_328():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt(a + b*cos(c + d*x))) - (2*A*b - 2*B*a)*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(2*A*b - 2*B*a)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(b*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_329():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_330():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = (Integer(-1) * (((((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('b') * (((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('A') * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_331():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = ((((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('A') * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_332():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**4/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**3/(3*b*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*a*(5*A*a**2*b - 9*A*b**3 - 8*B*a**3 + 12*B*a*b**2)*sin(c + d*x)*cos(c + d*x)**2/(3*b**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - sqrt(a + b*cos(c + d*x))*(60*A*a**3*b - 100*A*a*b**3 - 96*B*a**4 + 142*B*a**2*b**2 - 6*B*b**4)*sin(c + d*x)*cos(c + d*x)/(15*b**3*d*(a**2 - b**2)**2) + sqrt(a + b*cos(c + d*x))*(80*A*a**4*b - 130*A*a**2*b**3 + 10*A*b**5 - 128*B*a**5 + 196*B*a**3*b**2 - 28*B*a*b**4)*sin(c + d*x)/(15*b**4*d*(a**2 - b**2)**2) + sqrt((a + b*cos(c + d*x))/(a + b))*(160*A*a**4*b - 160*A*a**2*b**3 - 10*A*b**5 - 256*B*a**5 + 232*B*a**3*b**2 + 34*B*a*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(15*b**5*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - sqrt(a + b*cos(c + d*x))*(160*A*a**5*b - 280*A*a**3*b**3 + 80*A*a*b**5 - 256*B*a**6 + 424*B*a**4*b**2 - 110*B*a**2*b**4 - 18*B*b**6)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(15*b**5*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_333():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**3/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*(3*A*a**2*b - 7*A*b**3 - 6*B*a**3 + 10*B*a*b**2)*sin(c + d*x)/(3*b**3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + 2*a*(A*b - B*a)*sin(c + d*x)*cos(c + d*x)**2/(3*b*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - sqrt(a + b*cos(c + d*x))*(2*A*a*b - 4*B*a**2 + 2*B*b**2)*sin(c + d*x)/(3*b**3*d*(a**2 - b**2)) - sqrt((a + b*cos(c + d*x))/(a + b))*(16*A*a**3*b - 18*A*a*b**3 - 32*B*a**4 + 32*B*a**2*b**2 + 2*B*b**4)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**4*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(16*A*a**4*b - 30*A*a**2*b**3 + 6*A*b**5 - 32*B*a**5 + 56*B*a**3*b**2 - 16*B*a*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**4*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_334():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*a**2*(A*b - B*a)*sin(c + d*x)/(3*b**2*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*a*(2*A*a**2*b - 6*A*b**3 - 5*B*a**3 + 9*B*a*b**2)*sin(c + d*x)/(3*b**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt((a + b*cos(c + d*x))/(a + b))*(4*A*a**2*b - 6*A*b**3 - 16*B*a**3 + 18*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - sqrt(a + b*cos(c + d*x))*(4*A*a**3*b - 12*A*a*b**3 - 16*B*a**4 + 30*B*a**2*b**2 - 6*B*b**4)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**3*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_335():
    f = (A + B*cos(c + d*x))*cos(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*(A*b - B*a)*sin(c + d*x)/(3*b*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + (2*A*a**2*b + 6*A*b**3 + 4*B*a**3 - 12*B*a*b**2)*sin(c + d*x)/(3*b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt((a + b*cos(c + d*x))/(a + b))*(2*A*a*b + 4*B*a**2 - 6*B*b**2)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - sqrt(a + b*cos(c + d*x))*(2*A*a**2*b + 6*A*b**3 + 4*B*a**3 - 12*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b**2*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_336():
    f = (A + B*cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -(8*A*a*b - 2*B*a**2 - 6*B*b**2)*sin(c + d*x)/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - (2*A*b - 2*B*a)*sin(c + d*x)/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) - sqrt((a + b*cos(c + d*x))/(a + b))*(2*A*b - 2*B*a)*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(8*A*a*b - 2*B*a**2 - 6*B*b**2)*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(3*b*d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_337():
    f = (A + B*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(3)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_338():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**2/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(14) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * (((Symbol('a'))**(Integer(3)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(26) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(14) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('A') * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_339():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**3/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = ((((Integer(33) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(170) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(104) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(60) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(27) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(20) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Integer(27) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(33) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(170) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(105) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(12) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(104) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(60) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(4) * Symbol('a') * Symbol('B')))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('A') * sympy.sec((Symbol('c') + (Symbol('d') * x))) * sympy.tan((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_340():
    f = (B*a + B*b*cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt((a + b*cos(c + d*x))/(a + b))*elliptic_f(c/2 + d*x/2, 2*b/(a + b))/(d*sqrt(a + b*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_341():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = (Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_342():
    f = (B*a + B*b*cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = -2*B*b*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + 2*B*sqrt(a + b*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2*b/(a + b))/(d*sqrt((a + b*cos(c + d*x))/(a + b))*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_343():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('B') * sympy.sqrt(((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.Function('EllipticPi')(Integer(2), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), ((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_344():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*B*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(9*d) + (18*A*a + 14*B*b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + (18*A*a + 14*B*b)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + (2*A*b + 2*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + (10*A*b + 10*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (10*A*b + 10*B*a)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_345():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(7*d) + (14*A*a + 10*B*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (14*A*a + 10*B*b)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + (2*A*b + 2*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + (6*A*b + 6*B*a)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_346():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sqrt(cos(c + d*x))
    F = 2*B*b*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + (10*A*a + 6*B*b)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (2*A*b + 2*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + (2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_347():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))/sqrt(cos(c + d*x))
    F = 2*B*b*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + (6*A*a + 2*B*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (2*A*b + 2*B*a)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_348():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - (2*A*a - 2*B*b)*elliptic_e(c/2 + d*x/2, 2)/d + (2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_349():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + (2*A*a + 6*B*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d) + (2*A*b + 2*B*a)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - (2*A*b + 2*B*a)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_350():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + (6*A*a + 10*B*b)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - (6*A*a + 10*B*b)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (2*A*b + 2*B*a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + (2*A*b + 2*B*a)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_351():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2)
    F = 2*B*b*(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(11*d) + 2*b*(11*A*b + 13*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(99*d) + (18*B*b**2 + 22*a*(2*A*b + B*a))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(77*d) + (90*B*b**2 + 110*a*(2*A*b + B*a))*sin(c + d*x)*sqrt(cos(c + d*x))/(231*d) + (90*B*b**2 + 110*a*(2*A*b + B*a))*elliptic_f(c/2 + d*x/2, 2)/(231*d) + (18*A*a**2 + 14*A*b**2 + 28*B*a*b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + (18*A*a**2 + 14*A*b**2 + 28*B*a*b)*elliptic_e(c/2 + d*x/2, 2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_352():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(9*d) + 2*b*(9*A*b + 11*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d) + (14*B*b**2 + 18*a*(2*A*b + B*a))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + (14*B*b**2 + 18*a*(2*A*b + B*a))*elliptic_e(c/2 + d*x/2, 2)/(15*d) + (14*A*a**2 + 10*A*b**2 + 20*B*a*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (14*A*a**2 + 10*A*b**2 + 20*B*a*b)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_353():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sqrt(cos(c + d*x))
    F = 2*B*b*(a + b*cos(c + d*x))*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*b*(7*A*b + 9*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + (10*B*b**2 + 14*a*(2*A*b + B*a))*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (10*B*b**2 + 14*a*(2*A*b + B*a))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + (10*A*a**2 + 6*A*b**2 + 12*B*a*b)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_354():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2/sqrt(cos(c + d*x))
    F = 2*B*b*(a + b*cos(c + d*x))*sin(c + d*x)*sqrt(cos(c + d*x))/(5*d) + 2*b*(5*A*b + 7*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(15*d) + (6*B*b**2 + 10*a*(2*A*b + B*a))*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (6*A*a**2 + 2*A*b**2 + 4*B*a*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_355():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + 2*B*b**2*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) - (2*A*a**2 - 2*A*b**2 - 4*B*a*b)*elliptic_e(c/2 + d*x/2, 2)/d + (12*A*a*b + 6*B*a**2 + 2*B*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_356():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(2*A*b + B*a)*sin(c + d*x)/(d*sqrt(cos(c + d*x))) + (2*A*a**2 + 6*A*b**2 + 12*B*a*b)*elliptic_f(c/2 + d*x/2, 2)/(3*d) - (4*A*a*b + 2*B*a**2 - 2*B*b**2)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_357():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a*(2*A*b + B*a)*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + (6*A*a**2 + 10*A*b**2 + 20*B*a*b)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - (6*A*a**2 + 10*A*b**2 + 20*B*a*b)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (4*A*a*b + 2*B*a**2 + 6*B*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_358():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*cos(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a + b*cos(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(11*d) + 2*b**2*(11*A*b + 15*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(7)/2)/(99*d) + 2*b*(33*A*a*b + 26*B*a**2 + 9*B*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(77*d) + (154*A*a**3 + 330*A*a*b**2 + 330*B*a**2*b + 90*B*b**3)*sin(c + d*x)*sqrt(cos(c + d*x))/(231*d) + (154*A*a**3 + 330*A*a*b**2 + 330*B*a**2*b + 90*B*b**3)*elliptic_f(c/2 + d*x/2, 2)/(231*d) + (54*A*a**2*b + 14*A*b**3 + 18*B*a**3 + 42*B*a*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + (54*A*a**2*b + 14*A*b**3 + 18*B*a**3 + 42*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_359():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sqrt(cos(c + d*x))
    F = 2*B*b*(a + b*cos(c + d*x))**2*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(9*d) + 2*b**2*(9*A*b + 13*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(5)/2)/(63*d) + 2*b*(27*A*a*b + 22*B*a**2 + 7*B*b**2)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(45*d) + (30*A*a**3 + 54*A*a*b**2 + 54*B*a**2*b + 14*B*b**3)*elliptic_e(c/2 + d*x/2, 2)/(15*d) + (42*A*a**2*b + 10*A*b**3 + 14*B*a**3 + 30*B*a*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (42*A*a**2*b + 10*A*b**3 + 14*B*a**3 + 30*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_360():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3/sqrt(cos(c + d*x))
    F = 2*B*b*(a + b*cos(c + d*x))**2*sin(c + d*x)*sqrt(cos(c + d*x))/(7*d) + 2*b**2*(7*A*b + 11*B*a)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(35*d) + 2*b*(21*A*a*b + 18*B*a**2 + 5*B*b**2)*sin(c + d*x)*sqrt(cos(c + d*x))/(21*d) + (42*A*a**3 + 42*A*a*b**2 + 42*B*a**2*b + 10*B*b**3)*elliptic_f(c/2 + d*x/2, 2)/(21*d) + (30*A*a**2*b + 6*A*b**3 + 10*B*a**3 + 18*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_361():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3/cos(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**2*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*b**2*(5*A*a - B*b)*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) - 2*b*(6*A*a**2 - A*b**2 - 3*B*a*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) - (10*A*a**3 - 30*A*a*b**2 - 30*B*a**2*b - 6*B*b**3)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (18*A*a**2*b + 2*A*b**3 + 6*B*a**3 + 6*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_362():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**2*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(7*A*b + 3*B*a)*sin(c + d*x)/(3*d*sqrt(cos(c + d*x))) - 2*b**2*(A*a - B*b)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + (2*A*a**3 + 18*A*a*b**2 + 18*B*a**2*b + 2*B*b**3)*elliptic_f(c/2 + d*x/2, 2)/(3*d) - (6*A*a**2*b - 2*A*b**3 + 2*B*a**3 - 6*B*a*b**2)*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_363():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**2*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + 2*a**2*(9*A*b + 5*B*a)*sin(c + d*x)/(15*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*a*(3*A*a**2 + 14*A*b**2 + 15*B*a*b)*sin(c + d*x)/(5*d*sqrt(cos(c + d*x))) - (6*A*a**3 + 30*A*a*b**2 + 30*B*a**2*b - 10*B*b**3)*elliptic_e(c/2 + d*x/2, 2)/(5*d) + (6*A*a**2*b + 6*A*b**3 + 2*B*a**3 + 18*B*a*b**2)*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_364():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))
    F = (Integer(-1) * ((Integer(2) * ((Integer(5) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(5) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(4)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_365():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_366():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))
    F = ((Integer(2) * Symbol('B') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_367():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = ((Integer(2) * Symbol('B') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_368():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(2) * Symbol('A') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('A') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_369():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_370():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**2
    F = ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(2) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_371():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((((Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_372():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**2
    F = ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a') * Symbol('A') * Symbol('b')) + ((Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_373():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**2*sqrt(cos(c + d*x)))
    F = (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_374():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_375():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2))
    F = ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_376():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(33) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(38) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(35) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_377():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((((Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('B')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5)))) + (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_378():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**3
    F = ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(5)))) + ((Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_379():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**3*sqrt(cos(c + d*x)))
    F = (Integer(-1) * ((((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_380():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**3*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_381():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**3*cos(c + d*x)**(sympy.S(5)/2))
    F = ((((Integer(24) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(33) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Integer(63) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(86) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(38) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(61) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(33) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(24) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(65) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(8) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(13) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_382():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))
    F = 2*B*sin(c + d*x)*cos(c + d*x)**(sympy.S(3)/2)/(5*d) + 6*B*elliptic_e(c/2 + d*x/2, 2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_383():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))
    F = 2*B*sin(c + d*x)*sqrt(cos(c + d*x))/(3*d) + 2*B*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_384():
    f = (B*a + B*b*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))
    F = 2*B*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_385():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = 2*B*elliptic_f(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_386():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*B*sin(c + d*x)/(d*sqrt(cos(c + d*x))) - 2*B*elliptic_e(c/2 + d*x/2, 2)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_387():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*B*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*B*elliptic_f(c/2 + d*x/2, 2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_388():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('B') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Symbol('b'))**(Integer(2))) * Symbol('B') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_389():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**2
    F = ((Integer(2) * Symbol('B') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('B') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('B') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_390():
    f = (B*a + B*b*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**2
    F = ((Integer(2) * Symbol('B') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * Symbol('B') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_391():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))**2*sqrt(cos(c + d*x)))
    F = (Integer(2) * Symbol('B') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_392():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((Integer(2) * Symbol('B') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('B') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_393():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))**2*cos(c + d*x)**(sympy.S(5)/2))
    F = ((Integer(2) * Symbol('b') * Symbol('B') * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B') * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_394():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(6) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * ((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B'))) + (Integer(8) * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Integer(6) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_395():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + ((Symbol('a') + (Integer(2) * Symbol('b'))) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_396():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))/sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A')) + Symbol('B')) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_397():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('A') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * (Symbol('A') + (Integer(-1) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_398():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(A - 3*B)*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(A*b + 3*B*a)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_399():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(2*A*b + 10*B*a)*sin(c + d*x)/(15*a*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a + 2*A*b - 5*B*a)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a**2 - 2*A*b**2 + 5*B*a*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_400():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*cos(c + d*x))*(2*A*b + 14*B*a)*sin(c + d*x)/(35*a*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(50*A*a**2 - 8*A*b**2 + 14*B*a*b)*sin(c + d*x)/(105*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(8*A*b**2 + a**2*(25*A - 63*B) + 2*a*b*(3*A - 7*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**3*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(19*A*a**2*b + 8*A*b**3 + 63*B*a**3 - 14*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_401():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(156) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(4) * Symbol('A')) + Symbol('B')))) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(3)) * ((Integer(16) * Symbol('A')) + (Integer(9) * Symbol('B'))))) + (Integer(-1) * (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(28) * Symbol('A')) + (Integer(39) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(96) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(-1) * (Integer(24) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(48) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(156) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(8) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(12) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(32) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_402():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(12) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(14) * Symbol('a') * Symbol('b') * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(7) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_403():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * Symbol('a') * Symbol('A')) + (Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B')) + (Integer(2) * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(12) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_404():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * (Symbol('A') + (Integer(-1) * Symbol('B')))) + (Integer(-1) * (Symbol('b') * ((Integer(4) * Symbol('A')) + Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_405():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + ((Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(-1) * (Integer(3) * Symbol('B'))))) + (Integer(-1) * (Symbol('a') * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * Symbol('b') * Symbol('B'))))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_406():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(12*A*b + 10*B*a)*sin(c + d*x)/(15*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a - 3*A*b - 5*B*a + 15*B*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a**2 + 3*A*b**2 + 20*B*a*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_407():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*cos(c + d*x))*(16*A*b + 14*B*a)*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(50*A*a**2 + 6*A*b**2 + 84*B*a*b)*sin(c + d*x)/(105*a*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(6*A*b**2 - a**2*(25*A - 63*B) + 3*a*b*(19*A - 7*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(82*A*a**2*b - 6*A*b**3 + 63*B*a**3 + 21*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_408():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2)) + sqrt(a + b*cos(c + d*x))*(20*A*b + 18*B*a)*sin(c + d*x)/(63*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*cos(c + d*x))*(98*A*a**2 + 6*A*b**2 + 144*B*a*b)*sin(c + d*x)/(315*a*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(176*A*a**2*b - 8*A*b**3 + 150*B*a**3 + 18*B*a*b**2)*sin(c + d*x)/(315*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(8*A*b**3 - a**3*(147*A - 75*B) + 3*a**2*b*(13*A - 57*B) + 6*a*b**2*(A - 3*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**3*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*A*a**4 + 33*A*a**2*b**2 + 8*A*b**4 + 246*B*a**3*b - 18*B*a*b**3)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_409():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(150) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(2840) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(45) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(1692) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(1024) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(1920) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(45) * (Symbol('a'))**(Integer(4)) * Symbol('B')) + (Integer(-1) * (Integer(30) * (Symbol('a'))**(Integer(3)) * Symbol('b') * ((Integer(5) * Symbol('A')) + Symbol('B')))) + (Integer(-1) * (Integer(16) * (Symbol('b'))**(Integer(4)) * ((Integer(45) * Symbol('A')) + (Integer(64) * Symbol('B'))))) + (Integer(-1) * (Integer(8) * Symbol('a') * (Symbol('b'))**(Integer(3)) * ((Integer(355) * Symbol('A')) + (Integer(193) * Symbol('B'))))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(295) * Symbol('A')) + (Integer(423) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(1920) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(10) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(240) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(96) * Symbol('A') * (Symbol('b'))**(Integer(5)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(-1) * (Integer(40) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(240) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(128) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Integer(150) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(2840) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(45) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(1692) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(1024) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1920) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(50) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(120) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(172) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(320) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((((Integer(50) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(64) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(240) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((((Integer(10) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(40) * Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_410():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(8) * (Symbol('b'))**(Integer(3)) * ((Integer(16) * Symbol('A')) + (Integer(9) * Symbol('B')))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(132) * Symbol('A')) + (Integer(59) * Symbol('B')))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(52) * Symbol('A')) + (Integer(71) * Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(40) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(160) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(120) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(48) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(24) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(12) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(32) * Symbol('d')))**(Integer(-1))) + ((((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(11) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_411():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/sqrt(cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Integer(3) * Symbol('A')) + (Integer(4) * Symbol('B')))) + ((Symbol('a'))**(Integer(2)) * ((Integer(48) * Symbol('A')) + (Integer(33) * Symbol('B')))) + (Symbol('a') * ((Integer(54) * Symbol('A') * Symbol('b')) + (Integer(26) * Symbol('b') * Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(30) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_412():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(3)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(-1) * Symbol('B')))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('A')) + Symbol('B')))) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('b') * ((Integer(8) * Symbol('A')) + (Integer(3) * Symbol('B')))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(20) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(4) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_413():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(5)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(14) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * Symbol('b') * ((Integer(7) * Symbol('A')) + (Integer(-1) * (Integer(9) * Symbol('B'))))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(-1) * (Integer(3) * Symbol('B')))))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(6) * Symbol('A')) + Symbol('B'))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(14) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_414():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(7)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(23) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(35) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(15) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(23) * Symbol('A')) + (Integer(-1) * (Integer(45) * Symbol('B')))))) + ((Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(17) * Symbol('A')) + (Integer(-1) * (Integer(35) * Symbol('B'))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * ((Integer(9) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('B'))))))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(15) * Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_415():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(7*d*cos(c + d*x)**(sympy.S(7)/2)) + 2*a*sqrt(a + b*cos(c + d*x))*(10*A*b + 7*B*a)*sin(c + d*x)/(35*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(50*A*a**2 + 90*A*b**2 + 154*B*a*b)*sin(c + d*x)/(105*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(a**2*(25*A - 63*B) - 8*a*b*(15*A - 7*B) + 15*b**2*(A - 7*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(145*A*a**2*b + 15*A*b**3 + 63*B*a**3 + 161*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_416():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(9*d*cos(c + d*x)**(sympy.S(9)/2)) + 2*a*sqrt(a + b*cos(c + d*x))*(4*A*b + 3*B*a)*sin(c + d*x)/(21*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*cos(c + d*x))*(98*A*a**2 + 150*A*b**2 + 270*B*a*b)*sin(c + d*x)/(315*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(326*A*a**2*b + 10*A*b**3 + 150*B*a**3 + 270*B*a*b**2)*sin(c + d*x)/(315*a*d*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(10*A*b**3 + 3*a**3*(49*A - 25*B) - 6*a**2*b*(19*A - 60*B) + 15*a*b**2*(11*A - 3*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*A*a**4 + 279*A*a**2*b**2 - 10*A*b**4 + 435*B*a**3*b + 45*B*a*b**3)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_417():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/cos(c + d*x)**(sympy.S(13)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(11*d*cos(c + d*x)**(sympy.S(11)/2)) + 2*a*sqrt(a + b*cos(c + d*x))*(14*A*b + 11*B*a)*sin(c + d*x)/(99*d*cos(c + d*x)**(sympy.S(9)/2)) + sqrt(a + b*cos(c + d*x))*(162*A*a**2 + 226*A*b**2 + 418*B*a*b)*sin(c + d*x)/(693*d*cos(c + d*x)**(sympy.S(7)/2)) + sqrt(a + b*cos(c + d*x))*(2290*A*a**2*b + 30*A*b**3 + 1078*B*a**3 + 1650*B*a*b**2)*sin(c + d*x)/(3465*a*d*cos(c + d*x)**(sympy.S(5)/2)) + sqrt(a + b*cos(c + d*x))*(1350*A*a**4 + 2050*A*a**2*b**2 - 40*A*b**4 + 3586*B*a**3*b + 110*B*a*b**3)*sin(c + d*x)/(3465*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(40*A*b**4 + 3*a**4*(225*A - 539*B) - 6*a**3*b*(505*A - 209*B) + 15*a**2*b**2*(19*A - 121*B) + 10*a*b**3*(3*A - 11*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3465*a**3*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3705*A*a**4*b + 255*A*a**2*b**3 + 40*A*b**5 + 1617*B*a**5 + 3069*B*a**3*b**2 - 110*B*a*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3465*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_418():
    f = (a + b*cos(c + d*x))**(sympy.S(5)/2)*(B*cos(c + d*x) + 3*B*b/(2*a))/cos(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('a'))**(Integer(2)) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a') + (Integer(-1) * (Integer(3) * Symbol('b')))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Symbol('a') * Symbol('b'))) + (Integer(3) * (Symbol('b'))**(Integer(2)))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(5) * Symbol('a')) + ((Integer(3) * (Symbol('b'))**(Integer(2))) * (Symbol('a'))**(Integer(-1)))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (Symbol('d'))**(Integer(-1)))) + ((Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_419():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B'))) + (Integer(2) * Symbol('b') * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_420():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/sqrt(a + b*cos(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_421():
    f = (A + B*cos(c + d*x))/(sqrt(a + b*cos(c + d*x))*sqrt(cos(c + d*x)))
    F = ((Integer(2) * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_422():
    f = (A + B*cos(c + d*x))/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*A*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d) - 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(A - B)*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_423():
    f = (A + B*cos(c + d*x))/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(3*a*d*cos(c + d*x)**(sympy.S(3)/2)) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*A*b + a*(A - 3*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(2*A*b - 3*B*a)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_424():
    f = (A + B*cos(c + d*x))/(sqrt(a + b*cos(c + d*x))*cos(c + d*x)**(sympy.S(7)/2))
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)/(5*a*d*cos(c + d*x)**(sympy.S(5)/2)) - sqrt(a + b*cos(c + d*x))*(8*A*b - 10*B*a)*sin(c + d*x)/(15*a**2*d*cos(c + d*x)**(sympy.S(3)/2)) - 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(8*A*b**2 + a**2*(9*A - 5*B) - 2*a*b*(A + 5*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**3*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a**2 + 8*A*b**2 - 10*B*a*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_425():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (((Integer(3) * Symbol('a')) + Symbol('b')) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_426():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_427():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = -(2*A*b - 2*B*a)*sin(c + d*x)/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*A + 2*B)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*A*b - 2*B*a)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_428():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*sqrt(cos(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(4*A*b + 2*a*(A - B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*A*a**2 - 4*A*b**2 + 2*B*a*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**3*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_429():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*cos(c + d*x))*(2*A*a**2 - 8*A*b**2 + 6*B*a*b)*sin(c + d*x)/(3*a**2*d*(a**2 - b**2)*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a + 4*b)*(4*A*b + a*(A - 3*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(a + b)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(10*A*a**2*b - 16*A*b**3 - 6*B*a**3 + 12*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*sqrt(a + b))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_430():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = ((((Integer(6) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(2) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(12) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B'))) + (Integer(21) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(4)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.cos((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(6) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_431():
    f = (A + B*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = ((Integer(2) * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B'))) + (Integer(6) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_432():
    f = (A + B*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = (6*A*a**2 + 2*A*b**2 - 8*B*a*b)*sin(c + d*x)/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) - (2*A*b - 2*B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*A*a - 2*A*b + 2*B*a - 6*B*b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*A*a**2 + 2*A*b**2 - 8*B*a*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_433():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(cos(c + d*x)))
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sqrt(cos(c + d*x))/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) - (12*A*a**2*b - 4*A*b**3 - 6*B*a**3 - 2*B*a*b**2)*sin(c + d*x)/(3*a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(4*A*b**2 - 6*a**2*(A + B) + 2*a*b*(3*A + B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*sqrt(a + b)*(a**2 - b**2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(12*A*a**2*b - 4*A*b**3 - 6*B*a**3 - 2*B*a*b**2)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_434():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(cos(c + d*x))) + 2*b*(8*A*a**2*b - 4*A*b**3 - 5*B*a**3 + B*a*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*sqrt(cos(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(16*A*b**3 - 6*a**3*(A - B) - 6*a**2*b*(3*A + B) + 4*a*b**2*(3*A - B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(a + b)*(a**2 - b**2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*A*a**4 - 30*A*a**2*b**2 + 16*A*b**4 + 12*B*a**3*b - 4*B*a*b**3)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_435():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**(sympy.S(5)/2))
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*cos(c + d*x)**(sympy.S(3)/2)) + 2*b*(10*A*a**2*b - 6*A*b**3 - 7*B*a**3 + 3*B*a*b**2)*sin(c + d*x)/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2*cos(c + d*x)**(sympy.S(3)/2)) + sqrt(a + b*cos(c + d*x))*(2*A*a**4 - 26*A*a**2*b**2 + 16*A*b**4 + 16*B*a**3*b - 8*B*a*b**3)*sin(c + d*x)/(3*a**3*d*(a**2 - b**2)**2*cos(c + d*x)**(sympy.S(3)/2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(32*A*b**4 - 2*a**4*(A - 3*B) - 18*a**3*b*(A - B) - 4*a**2*b**2*(8*A + 3*B) + 8*a*b**3*(3*A - 2*B))*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*sqrt(a + b)*(a**2 - b**2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(16*A*a**4*b - 56*A*a**2*b**3 + 32*A*b**5 - 6*B*a**5 + 30*B*a**3*b**2 - 16*B*a*b**4)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**5*d*(a - b)*(a + b)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_436():
    f = (B*a + B*b*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * Symbol('d')))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('a') * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_437():
    f = (B*a + B*b*cos(c + d*x))*sqrt(cos(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.cot((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_438():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(cos(c + d*x)))
    F = 2*B*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_439():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -2*B*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*cot(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d) + B*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*cot(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_440():
    f = (cos(c + d*x) + 1)/(sqrt(3*cos(c + d*x) + 2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(1 - sec(c + d*x))*sqrt(-sec(c + d*x) - 1)*cot(c + d*x)*elliptic_e(asin(sqrt(5)*sqrt(3*cos(c + d*x) + 2)/(5*sqrt(cos(c + d*x)))), 5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_441():
    f = (cos(c + d*x) + 1)/(sqrt(3*cos(c + d*x) - 2)*cos(c + d*x)**(sympy.S(3)/2))
    F = -sqrt(5)*sqrt(sec(c + d*x) - 1)*sqrt(sec(c + d*x) + 1)*cot(c + d*x)*elliptic_e(asin(sqrt(3*cos(c + d*x) - 2)/sqrt(cos(c + d*x))), sympy.S(1)/5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_442():
    f = (cos(c + d*x) + 1)/(sqrt(2 - 3*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = sqrt(5)*sqrt(-cos(c + d*x))*sqrt(sec(c + d*x) - 1)*sqrt(sec(c + d*x) + 1)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(2 - 3*cos(c + d*x))/sqrt(-cos(c + d*x))), sympy.S(1)/5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_443():
    f = (cos(c + d*x) + 1)/(sqrt(-3*cos(c + d*x) - 2)*cos(c + d*x)**(sympy.S(3)/2))
    F = sqrt(-cos(c + d*x))*sqrt(1 - sec(c + d*x))*sqrt(-sec(c + d*x) - 1)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(5)*sqrt(-3*cos(c + d*x) - 2)/(5*sqrt(-cos(c + d*x)))), 5)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_444():
    f = (cos(c + d*x) + 1)/(sqrt(2*cos(c + d*x) + 3)*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*sqrt(1 - sec(c + d*x))*sqrt(sec(c + d*x) + 1)*cot(c + d*x)*elliptic_e(asin(sqrt(5)*sqrt(2*cos(c + d*x) + 3)/(5*sqrt(cos(c + d*x)))), -5)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_445():
    f = (cos(c + d*x) + 1)/(sqrt(3 - 2*cos(c + d*x))*cos(c + d*x)**(sympy.S(3)/2))
    F = 2*sqrt(5)*sqrt(1 - sec(c + d*x))*sqrt(sec(c + d*x) + 1)*cot(c + d*x)*elliptic_e(asin(sqrt(3 - 2*cos(c + d*x))/sqrt(cos(c + d*x))), sympy.S(-1)/5)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_446():
    f = (cos(c + d*x) + 1)/(sqrt(2*cos(c + d*x) - 3)*cos(c + d*x)**(sympy.S(3)/2))
    F = -2*sqrt(5)*sqrt(-cos(c + d*x))*sqrt(1 - sec(c + d*x))*sqrt(sec(c + d*x) + 1)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(2*cos(c + d*x) - 3)/sqrt(-cos(c + d*x))), sympy.S(-1)/5)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_447():
    f = (cos(c + d*x) + 1)/(sqrt(-2*cos(c + d*x) - 3)*cos(c + d*x)**(sympy.S(3)/2))
    F = -2*sqrt(-cos(c + d*x))*sqrt(1 - sec(c + d*x))*sqrt(sec(c + d*x) + 1)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(5)*sqrt(-2*cos(c + d*x) - 3)/(5*sqrt(-cos(c + d*x)))), -5)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_448():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**n
    F = sympy.Function('Unintegrable')((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Symbol('n')) * (Symbol('A') + (Symbol('B') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_449():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**4
    F = B*b*(c*cos(e + f*x))**(m + 1)*(a + b*cos(e + f*x))**3*sin(e + f*x)/(c*f*(m + 5)) + b**2*(c*cos(e + f*x))**(m + 1)*(2*A*a*b*(m + 5)**2 + B*a**2*(m**2 + 11*m + 36) + B*b**2*(m + 4)**2)*sin(e + f*x)*cos(e + f*x)/(c*f*(m + 3)*(m + 4)*(m + 5)) + b*(c*cos(e + f*x))**(m + 1)*(a + b*cos(e + f*x))**2*(A*b*(m + 5) + B*a*(m + 8))*sin(e + f*x)/(c*f*(m + 4)*(m + 5)) + b*(c*cos(e + f*x))**(m + 1)*(A*a**2*b*(5*m**2 + 47*m + 110) + A*b**3*(m**2 + 8*m + 15) + 2*B*a**3*(m**2 + 10*m + 28) + 4*B*a*b**2*(m**2 + 8*m + 15))*sin(e + f*x)/(c*f*(m + 2)*(m + 4)*(m + 5)) - (c*cos(e + f*x))**(m + 1)*(A*a**4*(m**2 + 6*m + 8) + 6*A*a**2*b**2*(m**2 + 5*m + 4) + A*b**4*(m**2 + 4*m + 3) + 4*B*a**3*b*(m**2 + 5*m + 4) + 4*B*a*b**3*(m**2 + 4*m + 3))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(c*f*(m + 1)*(m + 2)*(m + 4)*sqrt(sin(e + f*x)**2)) - (c*cos(e + f*x))**(m + 2)*(4*A*a**3*b*(m**2 + 8*m + 15) + 4*A*a*b**3*(m**2 + 7*m + 10) + B*a**4*(m**2 + 8*m + 15) + 6*B*a**2*b**2*(m**2 + 7*m + 10) + B*b**4*(m**2 + 6*m + 8))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(e + f*x)**2)/(c**2*f*(m + 2)*(m + 3)*(m + 5)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_450():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**3
    F = B*b*(c*cos(e + f*x))**(m + 1)*(a + b*cos(e + f*x))**2*sin(e + f*x)/(c*f*(m + 4)) + b**2*(c*cos(e + f*x))**(m + 1)*(A*b*(m + 4) + B*a*(m + 6))*sin(e + f*x)*cos(e + f*x)/(c*f*(m + 3)*(m + 4)) + b*(c*cos(e + f*x))**(m + 1)*(3*A*a*b*(m + 4) + 2*B*a**2*(m + 5) + B*b**2*(m + 3))*sin(e + f*x)/(c*f*(m + 2)*(m + 4)) - (c*cos(e + f*x))**(m + 1)*(a**2*(m + 2)*(A*a*(m + 4) + B*b*(m + 1)) + b*(m + 1)*(3*A*a*b*(m + 4) + 2*B*a**2*(m + 5) + B*b**2*(m + 3)))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(c*f*(m + 1)*(m + 2)*(m + 4)*sqrt(sin(e + f*x)**2)) - (c*cos(e + f*x))**(m + 2)*(3*A*a**2*b*(m + 3) + A*b**3*(m + 2) + B*a**3*(m + 3) + 3*B*a*b**2*(m + 2))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(e + f*x)**2)/(c**2*f*(m + 2)*(m + 3)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_451():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**2
    F = B*b*(c*cos(e + f*x))**(m + 1)*(a + b*cos(e + f*x))*sin(e + f*x)/(c*f*(m + 3)) + b*(c*cos(e + f*x))**(m + 1)*(A*b*(m + 3) + B*a*(m + 4))*sin(e + f*x)/(c*f*(m + 2)*(m + 3)) - (c*cos(e + f*x))**(m + 1)*(A*a**2*(m + 2) + A*b**2*(m + 1) + 2*B*a*b*(m + 1))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(c*f*(m + 1)*(m + 2)*sqrt(sin(e + f*x)**2)) - (c*cos(e + f*x))**(m + 2)*(B*b**2*(m + 2) + a*(m + 3)*(2*A*b + B*a))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(e + f*x)**2)/(c**2*f*(m + 2)*(m + 3)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_452():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))
    F = B*b*(c*cos(e + f*x))**(m + 1)*sin(e + f*x)/(c*f*(m + 2)) - (c*cos(e + f*x))**(m + 1)*(A*a*(m + 2) + B*b*(m + 1))*sin(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(c*f*(m + 1)*(m + 2)*sqrt(sin(e + f*x)**2)) - (c*cos(e + f*x))**(m + 2)*(A*b + B*a)*sin(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (m/2 + 2,), cos(e + f*x)**2)/(c**2*f*(m + 2)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_453():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))/(a + b*cos(e + f*x))
    F = -B*(c*cos(e + f*x))**(m + 1)*sin(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(b*c*f*(m + 1)*sqrt(sin(e + f*x)**2)) + a*c*(c*cos(e + f*x))**(m - 1)*(A*b - B*a)*(cos(e + f*x)**2)**(sympy.S.Half - m/2)*sin(e + f*x)*appellf1(sympy.S.Half, 1, sympy.S.Half - m/2, sympy.S(3)/2, -b**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(b*f*(a**2 - b**2)) - (c*cos(e + f*x))**m*(A*b - B*a)*sin(e + f*x)*appellf1(sympy.S.Half, 1, -m/2, sympy.S(3)/2, -b**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(f*(a**2 - b**2)*(cos(e + f*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_454():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('b') * Symbol('B') * ((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**((Integer(1) + Symbol('m'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('c') * Symbol('f') * (Integer(5) + (Integer(2) * Symbol('m')))))**(Integer(-1))) + ((Integer(2) * sympy.Function('Unintegrable')(((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * (((Integer(2))**(Integer(-1)) * Symbol('a') * Symbol('c') * ((Integer(2) * Symbol('b') * Symbol('B') * (Integer(1) + Symbol('m'))) + (Integer(2) * Symbol('a') * Symbol('A') * ((Integer(5) * (Integer(2))**(Integer(-1))) + Symbol('m'))))) + ((Integer(2))**(Integer(-1)) * Symbol('c') * (((Symbol('b'))**(Integer(2)) * Symbol('B') * (Integer(3) + (Integer(2) * Symbol('m')))) + (Symbol('a') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * (Integer(5) + (Integer(2) * Symbol('m'))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * ((Integer(2) * Symbol('a') * Symbol('B') * (Integer(3) + Symbol('m'))) + (Symbol('A') * Symbol('b') * (Integer(5) + (Integer(2) * Symbol('m'))))) * (sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))), x)) * ((Symbol('c') * (Integer(5) + (Integer(2) * Symbol('m')))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_455():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))*sqrt(a + b*cos(e + f*x))
    F = sympy.Function('Unintegrable')((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (Symbol('A') + (Symbol('B') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_456():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))/sqrt(a + b*cos(e + f*x))
    F = sympy.Function('Unintegrable')(((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * (Symbol('A') + (Symbol('B') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_457():
    f = (c*cos(e + f*x))**m*(A + B*cos(e + f*x))/(a + b*cos(e + f*x))**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * ((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**((Integer(1) + Symbol('m'))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('c') * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * sympy.Function('Unintegrable')(((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * (((Integer(2))**(Integer(-1)) * Symbol('c') * ((Symbol('a') * ((Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B'))))) + (Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * ((Integer(2))**(Integer(-1)) + Symbol('m'))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * Symbol('c') * (Integer(3) + (Integer(2) * Symbol('m'))) * (sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))))) * (sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))), x)) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_458():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a*(A + B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*a*(3*A + 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 2*a*(3*A + 5*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_459():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*(A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*(A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_460():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*a*(A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_461():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))
    F = 2*B*a*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*a*(3*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_462():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)/sqrt(sec(c + d*x))
    F = 2*B*a*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a*(A + B)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*a*(5*A + 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_463():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*a*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*a*(A + B)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*a*(A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a*(7*A + 5*B)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 2*a*(7*A + 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_464():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 4*a**2*(A + 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(4*A + 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 4*a**2*(4*A + 5*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a**2*(7*A + 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_465():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(5)/2)
    F = -4*A*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*A*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(2*A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 2*a**2*(5*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_466():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 4*B*a**2*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 2*B*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a**2*(3*A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(3*A + 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_467():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2*sqrt(sec(c + d*x))
    F = 2*B*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**2*(2*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 4*a**2*(5*A + 4*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 2*a**2*(5*A + 7*B)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_468():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**2/sqrt(sec(c + d*x))
    F = 2*B*(a**2*sec(c + d*x) + a**2)*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 4*a**2*(4*A + 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**2*(7*A + 6*B)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 4*a**2*(7*A + 6*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 2*a**2*(7*A + 9*B)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_469():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(7*d) + 4*a**3*(7*A + 9*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - 4*a**3*(7*A + 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(13*A + 21*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**3*(41*A + 42*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*d) + (22*A + 14*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_470():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(3*A + 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - 4*a**3*(9*A + 5*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(21*A + 20*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) + (18*A + 10*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_471():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) - 4*a**3*(A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + 20*a**3*(A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 4*a**3*(4*A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + (2*A - 2*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_472():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 4*a**3*(5*A - 6*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d) + 4*a**3*(5*A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + 4*a**3*(5*A + 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (10*A + 18*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_473():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3*sqrt(sec(c + d*x))
    F = 2*B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 4*a**3*(9*A + 7*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + 4*a**3*(21*A + 13*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**3*(42*A + 41*B)*sin(c + d*x)/(105*d*sqrt(sec(c + d*x))) + (14*A + 22*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_474():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**3/sqrt(sec(c + d*x))
    F = 2*B*a*(a*sec(c + d*x) + a)**2*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 4*a**3*(13*A + 11*B)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + 4*a**3*(13*A + 11*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + 4*a**3*(21*A + 17*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + 4*a**3*(24*A + 23*B)*sin(c + d*x)/(105*d*sec(c + d*x)**(sympy.S(3)/2)) + (18*A + 26*B)*(a**3*sec(c + d*x) + a**3)*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_475():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(d*(a*sec(c + d*x) + a)) - (3*A - 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d) + (3*A - 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (5*A - 3*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d) + (5*A - 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_476():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(d*(a*sec(c + d*x) + a)) - (A - B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (3*A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d) - (3*A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_477():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) + (A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_478():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    F = (A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*(a*sec(c + d*x) + a)) - (A - 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d) + (A - B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_479():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = (A - B)*sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sqrt(sec(c + d*x))) - (3*A - 5*B)*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) - (3*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d) + (3*A - 3*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_480():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)/(d*(a*sec(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) - (5*A - 7*B)*sin(c + d*x)/(5*a*d*sec(c + d*x)**(sympy.S(3)/2)) + (5*A - 5*B)*sin(c + d*x)/(3*a*d*sqrt(sec(c + d*x))) + (5*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a*d) - (15*A - 21*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_481():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**2
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3*d*(a*sec(c + d*x) + a)**2) + (4*A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d) - (4*A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - (5*A - 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) - (5*A - 2*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_482():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**2
    F = A*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - A*sin(c + d*x)*sqrt(sec(c + d*x))/(a**2*d*(sec(c + d*x) + 1)) - (A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*(a*sec(c + d*x) + a)**2) + (2*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_483():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**2*sqrt(sec(c + d*x)))
    F = -B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) - (A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*(a*sec(c + d*x) + a)**2) + (A + 2*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) + (A + 2*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_484():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(3)/2))
    F = (A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*(a*sec(c + d*x) + a)**2) - (A - 4*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + (2*A - 5*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d) + (2*A - 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*(sec(c + d*x) + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_485():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**2*sec(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)/(3*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) + (4*A - 7*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(a**2*d) + (4*A - 7*B)*sin(c + d*x)/(3*a**2*d*(sec(c + d*x) + 1)*sqrt(sec(c + d*x))) - (5*A - 10*B)*sin(c + d*x)/(3*a**2*d*sqrt(sec(c + d*x))) - (5*A - 10*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_486():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(5*d*(a*sec(c + d*x) + a)**3) - (13*A - 3*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*d*(a**3*sec(c + d*x) + a**3)) - (8*A - 3*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(15*a*d*(a*sec(c + d*x) + a)**2) - (13*A - 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) + (49*A - 9*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(10*a**3*d) - (49*A - 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_487():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**3
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*(a*sec(c + d*x) + a)**3) - (9*A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(10*d*(a**3*sec(c + d*x) + a**3)) - (6*A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*a*d*(a*sec(c + d*x) + a)**2) + (3*A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) + (9*A + B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_488():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**3*sqrt(sec(c + d*x)))
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d*(a*sec(c + d*x) + a)**3) + (A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) - (4*A + B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + (A - B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d) + (A + B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_489():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(3)/2))
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*(a*sec(c + d*x) + a)**3) + (A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) + (2*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + (A + 3*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) - (A + 9*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_490():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d*(a*sec(c + d*x) + a)**3) + (3*A - 13*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a**3*sec(c + d*x) + a**3)) + (3*A - 8*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*a*d*(a*sec(c + d*x) + a)**2) + (3*A - 13*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) - (9*A - 49*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_491():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**3*sec(c + d*x)**(sympy.S(7)/2))
    F = (A - B)*sin(c + d*x)/(5*d*(a*sec(c + d*x) + a)**3*sqrt(sec(c + d*x))) + (49*A - 119*B)*sin(c + d*x)/(30*d*(a**3*sec(c + d*x) + a**3)*sqrt(sec(c + d*x))) + (A - 2*B)*sin(c + d*x)/(3*a*d*(a*sec(c + d*x) + a)**2*sqrt(sec(c + d*x))) - (13*A - 33*B)*sin(c + d*x)/(6*a**3*d*sqrt(sec(c + d*x))) - (13*A - 33*B)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(6*a**3*d) + (49*A - 119*B)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(10*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_492():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(8*A + 9*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d*sqrt(a*cos(c + d*x) + a)) + 4*a*(8*A + 9*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) + 16*a*(8*A + 9*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*d*sqrt(a*cos(c + d*x) + a)) + 32*a*(8*A + 9*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_493():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(6*A + 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d*sqrt(a*cos(c + d*x) + a)) + 8*a*(6*A + 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) + 16*a*(6*A + 7*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_494():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(4*A + 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(a*cos(c + d*x) + a)) + 4*a*(4*A + 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_495():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) + 2*a*(2*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_496():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*cos(c + d*x) + a)) + 2*B*sqrt(a)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_497():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))
    F = B*a*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + sqrt(a)*(2*A + B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_498():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)/sqrt(sec(c + d*x))
    F = B*a*sin(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + sqrt(a)*(4*A + 3*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + a*(4*A + 3*B)*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_499():
    f = (A + B*cos(c + d*x))*sqrt(a*cos(c + d*x) + a)/sec(c + d*x)**(sympy.S(3)/2)
    F = B*a*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + sqrt(a)*(6*A + 5*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + a*(6*A + 5*B)*sin(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a*(6*A + 5*B)*sin(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_500():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(13)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(11)/2)/(11*d) + 2*a**2*(12*A + 11*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(99*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(168*A + 187*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(693*d*sqrt(a*cos(c + d*x) + a)) + 4*a**2*(168*A + 187*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(1155*d*sqrt(a*cos(c + d*x) + a)) + 16*a**2*(168*A + 187*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3465*d*sqrt(a*cos(c + d*x) + a)) + 32*a**2*(168*A + 187*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3465*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_501():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d) + 2*a**2*(10*A + 9*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(34*A + 39*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) + 8*a**2*(34*A + 39*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*d*sqrt(a*cos(c + d*x) + a)) + 16*a**2*(34*A + 39*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_502():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 2*a**2*(8*A + 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(52*A + 63*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) + 4*a**2*(52*A + 63*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_503():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*a**2*(6*A + 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(18*A + 25*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_504():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*B*a**(sympy.S(3)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d + 2*a**2*(4*A + 3*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_505():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sqrt(sec(c + d*x))/d + a**(sympy.S(3)/2)*(2*A + 3*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d - a**2*(2*A - B)*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_506():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(2*d*sqrt(sec(c + d*x))) + a**(sympy.S(3)/2)*(12*A + 7*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) + a**2*(4*A + 5*B)*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_507():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(3*d*sec(c + d*x)**(sympy.S(3)/2)) + a**(sympy.S(3)/2)*(14*A + 11*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + a**2*(6*A + 7*B)*sin(c + d*x)/(12*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + a**2*(14*A + 11*B)*sin(c + d*x)/(8*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_508():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = B*a*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(4*d*sec(c + d*x)**(sympy.S(5)/2)) + a**(sympy.S(3)/2)*(88*A + 75*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(64*d) + a**2*(8*A + 9*B)*sin(c + d*x)/(24*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + a**2*(88*A + 75*B)*sin(c + d*x)/(64*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a**2*(88*A + 75*B)*sin(c + d*x)/(96*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_509():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(15)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(13)/2)/(13*d) + 2*a**3*(280*A + 299*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(1287*d*sqrt(a*cos(c + d*x) + a)) + 2*a**3*(4184*A + 4615*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(9009*d*sqrt(a*cos(c + d*x) + a)) + 4*a**3*(4184*A + 4615*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(15015*d*sqrt(a*cos(c + d*x) + a)) + 16*a**3*(4184*A + 4615*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(45045*d*sqrt(a*cos(c + d*x) + a)) + 32*a**3*(4184*A + 4615*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(45045*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(16*A + 13*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(11)/2)/(143*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_510():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(13)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(11)/2)/(11*d) + 2*a**3*(194*A + 209*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(693*d*sqrt(a*cos(c + d*x) + a)) + 2*a**3*(710*A + 803*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(1155*d*sqrt(a*cos(c + d*x) + a)) + 8*a**3*(710*A + 803*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3465*d*sqrt(a*cos(c + d*x) + a)) + 16*a**3*(710*A + 803*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3465*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(14*A + 11*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(99*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_511():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d) + 2*a**3*(124*A + 135*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(315*d*sqrt(a*cos(c + d*x) + a)) + 2*a**3*(292*A + 345*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*d*sqrt(a*cos(c + d*x) + a)) + 4*a**3*(292*A + 345*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(4*A + 3*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_512():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 2*a**3*(10*A + 11*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(a*cos(c + d*x) + a)) + 2*a**3*(230*A + 301*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(10*A + 7*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_513():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + 2*B*a**(sympy.S(5)/2)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d + 2*a**3*(32*A + 35*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*cos(c + d*x) + a)) + 2*a**2*(8*A + 5*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_514():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + a**(sympy.S(5)/2)*(2*A + 5*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/d - a**3*(14*A + 3*B)*sin(c + d*x)/(3*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + 2*a**2*(2*A + B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_515():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)*sqrt(sec(c + d*x))/d + a**(sympy.S(5)/2)*(20*A + 19*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(4*d) - a**3*(4*A - 9*B)*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) - a**2*(4*A - B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_516():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + a**(sympy.S(5)/2)*(38*A + 25*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(8*d) + a**3*(54*A + 49*B)*sin(c + d*x)/(24*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a**2*(2*A + 3*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(4*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_517():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(4*d*sec(c + d*x)**(sympy.S(3)/2)) + a**(sympy.S(5)/2)*(200*A + 163*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(64*d) + a**3*(104*A + 95*B)*sin(c + d*x)/(96*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + a**3*(200*A + 163*B)*sin(c + d*x)/(64*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a**2*(8*A + 11*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(24*d*sec(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_518():
    f = (A + B*cos(c + d*x))*(a*cos(c + d*x) + a)**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = B*a*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(5)/2)) + a**(sympy.S(5)/2)*(326*A + 283*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(128*d) + a**3*(170*A + 157*B)*sin(c + d*x)/(240*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(5)/2)) + a**3*(326*A + 283*B)*sin(c + d*x)/(128*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + a**3*(326*A + 283*B)*sin(c + d*x)/(192*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + a**2*(10*A + 13*B)*sqrt(a*cos(c + d*x) + a)*sin(c + d*x)/(40*d*sec(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_519():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(11)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*A*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d*sqrt(a*cos(c + d*x) + a)) - (2*A - 18*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d*sqrt(a*cos(c + d*x) + a)) + (38*A - 6*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) - (58*A - 186*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*d*sqrt(a*cos(c + d*x) + a)) + (514*A - 258*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(315*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_520():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(9)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*A*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d*sqrt(a*cos(c + d*x) + a)) - (2*A - 14*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d*sqrt(a*cos(c + d*x) + a)) + (62*A - 14*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*d*sqrt(a*cos(c + d*x) + a)) - (86*A - 182*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(105*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_521():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*A*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d*sqrt(a*cos(c + d*x) + a)) - (2*A - 10*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d*sqrt(a*cos(c + d*x) + a)) + (26*A - 10*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(15*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_522():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*A*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d*sqrt(a*cos(c + d*x) + a)) - (2*A - 6*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_523():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/sqrt(a*cos(c + d*x) + a)
    F = 2*A*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_524():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/sqrt(a*cos(c + d*x) + a)
    F = 2*B*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_525():
    f = (A + B*cos(c + d*x))/(sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x)))
    F = B*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) - sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d) + (2*A - B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_526():
    f = (A + B*cos(c + d*x))/(sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2))
    F = B*sin(c + d*x)/(2*d*sqrt(a*cos(c + d*x) + a)*sec(c + d*x)**(sympy.S(3)/2)) + (4*A - B)*sin(c + d*x)/(4*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + sqrt(2)*(A - B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d) - (4*A - 7*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(4*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_527():
    f = (A*a + B*b*cos(c + d*x)**2 + (A*b + B*a)*cos(c + d*x))*sqrt(sec(c + d*x))/sqrt(a*cos(c + d*x) + a)
    F = B*b*sin(c + d*x)/(d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + sqrt(2)*(A - B)*(a - b)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(sqrt(a)*d) + (2*A*b + 2*B*a - B*b)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_528():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(9)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (11*A - 7*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(14*a*d*sqrt(a*cos(c + d*x) + a)) - (67*A - 63*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(70*a*d*sqrt(a*cos(c + d*x) + a)) + (397*A - 273*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(210*a*d*sqrt(a*cos(c + d*x) + a)) - (1201*A - 1029*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(210*a*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(19*A - 15*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_529():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (9*A - 5*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(10*a*d*sqrt(a*cos(c + d*x) + a)) - (39*A - 35*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(30*a*d*sqrt(a*cos(c + d*x) + a)) + (147*A - 95*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(30*a*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(15*A - 11*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_530():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (7*A - 3*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*a*d*sqrt(a*cos(c + d*x) + a)) - (19*A - 15*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*a*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(11*A - 7*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_531():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (5*A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(2*a*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(7*A - 3*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_532():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(3)/2)
    F = -(A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(3*A + B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_533():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = 2*B*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) + (A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(A - 5*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_534():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = (A - B)*sin(c + d*x)/(2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) - (A - 3*B)*sin(c + d*x)/(2*a*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + (2*A - 3*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(3)/2)*d) - sqrt(2)*(5*A - 9*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_535():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (21*A - 13*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (157*A - 85*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(80*a**2*d*sqrt(a*cos(c + d*x) + a)) - (787*A - 475*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(240*a**2*d*sqrt(a*cos(c + d*x) + a)) + (2671*A - 1495*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(240*a**2*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(283*A - 163*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_536():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (17*A - 9*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (95*A - 39*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(48*a**2*d*sqrt(a*cos(c + d*x) + a)) - (299*A - 147*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(48*a**2*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(163*A - 75*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_537():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (13*A - 5*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (49*A - 9*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(16*a**2*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(75*A - 19*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_538():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(5)/2)
    F = -(A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - (9*A - B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(19*A + 5*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_539():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = (A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + (A + 7*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(5*A + 3*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_540():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*B*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) + (A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)) + (3*A - 11*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(3*A - 43*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_541():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = (A - B)*sin(c + d*x)/(4*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)) + (7*A - 15*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) - (11*A - 35*B)*sin(c + d*x)/(16*a**2*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + (2*A - 5*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(5)/2)*d) - sqrt(2)*(43*A - 115*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_542():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -(A - B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - (23*A - 11*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (109*A - 41*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (579*A - 199*B)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(192*a**3*d*sqrt(a*cos(c + d*x) + a)) - (1887*A - 691*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(192*a**3*d*sqrt(a*cos(c + d*x) + a)) + sqrt(2)*(1015*A - 363*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_543():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -(A - B)*sin(c + d*x)*sqrt(sec(c + d*x))/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)) - (19*A - 7*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)) - (199*A - 43*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)) + (691*A - 103*B)*sin(c + d*x)*sqrt(sec(c + d*x))/(192*a**3*d*sqrt(a*cos(c + d*x) + a)) - sqrt(2)*(363*A - 63*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_544():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a*cos(c + d*x) + a)**(sympy.S(7)/2)
    F = -(A - B)*sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x))) - (5*A - B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - (103*A + 5*B)*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(63*A + 13*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_545():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x)))
    F = (A - B)*sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sqrt(sec(c + d*x))) + (A + 3*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) - (5*A - 17*B)*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(13*A + 7*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_546():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = (A - B)*sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(3)/2)) + (A - 13*B)*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sqrt(sec(c + d*x))) + (17*A + 67*B)*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(7*A + 5*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_547():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*B*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(7)/2)*d) + (A - B)*sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(5)/2)) + (5*A - 17*B)*sin(c + d*x)/(48*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)) + (5*A - 49*B)*sin(c + d*x)/(64*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) + sqrt(2)*(5*A - 177*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_548():
    f = (A + B*cos(c + d*x))/((a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(7)/2))
    F = (A - B)*sin(c + d*x)/(6*d*(a*cos(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**(sympy.S(7)/2)) + (3*A - 7*B)*sin(c + d*x)/(16*a*d*(a*cos(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)) + (79*A - 259*B)*sin(c + d*x)/(192*a**2*d*(a*cos(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)) - (49*A - 189*B)*sin(c + d*x)/(64*a**3*d*sqrt(a*cos(c + d*x) + a)*sqrt(sec(c + d*x))) + (2*A - 7*B)*sqrt(cos(c + d*x))*asin(sqrt(a)*sin(c + d*x)/sqrt(a*cos(c + d*x) + a))*sqrt(sec(c + d*x))/(a**(sympy.S(7)/2)*d) - sqrt(2)*(177*A - 637*B)*sqrt(cos(c + d*x))*atan(sqrt(2)*sqrt(a)*sin(c + d*x)/(2*sqrt(a*cos(c + d*x) + a)*sqrt(cos(c + d*x))))*sqrt(sec(c + d*x))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_549():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + (6*A*a + 10*B*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (6*A*a + 10*B*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (2*A*b + 2*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_550():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + (2*A*a + 6*B*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (-2*A*b - 2*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (2*A*b + 2*B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_551():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a*sin(c + d*x)*sqrt(sec(c + d*x))/d + (-2*A*a + 2*B*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_552():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))*sqrt(sec(c + d*x))
    F = 2*B*b*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + (6*A*a + 2*B*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_553():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))/sqrt(sec(c + d*x))
    F = 2*B*b*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + (10*A*a + 6*B*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (2*A*b + 2*B*a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + (2*A*b + 2*B*a)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_554():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + (14*A*a + 10*B*b)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (14*A*a + 10*B*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (2*A*b + 2*B*a)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + (6*A*b + 6*B*a)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_555():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*sec(c + d*x) + b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(5*d) + 2*a*(7*A*b + 5*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d) + (6*A*a**2 + 10*b*(A*b + 2*B*a))*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (6*A*a**2 + 10*b*(A*b + 2*B*a))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (4*A*a*b + 2*B*a**2 + 6*B*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_556():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*a*(a*sec(c + d*x) + b)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + 2*a*(5*A*b + 3*B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + (2*A*a**2 + 6*A*b**2 + 12*B*a*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - (4*A*a*b + 2*B*a**2 - 2*B*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_557():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*A*a**2*sin(c + d*x)*sqrt(sec(c + d*x))/d + 2*B*b**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) - (2*A*a**2 - 2*b*(A*b + 2*B*a))*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d + (12*A*a*b + 6*B*a**2 + 2*B*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_558():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2*sqrt(sec(c + d*x))
    F = 2*B*b**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*b*(A*b + 2*B*a)*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + (6*A*a**2 + 2*A*b**2 + 4*B*a*b)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) + (20*A*a*b + 10*B*a**2 + 6*B*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_559():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**2/sqrt(sec(c + d*x))
    F = 2*B*b**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*b*(A*b + 2*B*a)*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + (10*A*a**2 + 6*A*b**2 + 12*B*a*b)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (28*A*a*b + 14*B*a**2 + 10*B*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (28*A*a*b + 14*B*a**2 + 10*B*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_560():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a*sec(c + d*x) + b)**2*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(7*d) + 2*a**2*(11*A*b + 7*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + 2*a*(5*A*a**2 + 18*A*b**2 + 21*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(21*d) + (10*A*a**3 + 42*A*a*b**2 + 42*B*a**2*b + 42*B*b**3)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (18*A*a**2*b + 10*A*b**3 + 6*B*a**3 + 30*B*a*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (18*A*a**2*b + 10*A*b**3 + 6*B*a**3 + 30*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_561():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*(a*sec(c + d*x) + b)**2*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + 2*a**2*(9*A*b + 5*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d) + 2*a*(3*A*a**2 + 14*A*b**2 + 15*B*a*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) - (6*A*a**3 + 30*A*a*b**2 + 30*B*a**2*b - 10*B*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (6*A*a**2*b + 6*A*b**3 + 2*B*a**3 + 18*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_562():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*B*b*(a*sec(c + d*x) + b)**2*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*a**2*(A*a - B*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*a*(9*A*a*b + 3*B*a**2 - 2*B*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d) + (2*A*a**3 + 18*A*a*b**2 + 18*B*a**2*b + 2*B*b**3)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d) - (6*A*a**2*b - 2*A*b**3 + 2*B*a**3 - 6*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_563():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2)
    F = 2*B*b*(a*sec(c + d*x) + b)**2*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*a**2*(5*A*a - B*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(5*d) + 2*b**2*(5*A*b + 9*B*a)*sin(c + d*x)/(15*d*sqrt(sec(c + d*x))) - (10*A*a**3 - 30*A*a*b**2 - 30*B*a**2*b - 6*B*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d) + (18*A*a**2*b + 2*A*b**3 + 6*B*a**3 + 6*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_564():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3*sqrt(sec(c + d*x))
    F = 2*B*b*(a*sec(c + d*x) + b)**2*sin(c + d*x)/(7*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*b**2*(7*A*b + 11*B*a)*sin(c + d*x)/(35*d*sec(c + d*x)**(sympy.S(3)/2)) + 2*b*(21*A*a*b + 18*B*a**2 + 5*B*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (42*A*a**3 + 42*A*a*b**2 + 42*B*a**2*b + 10*B*b**3)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d) + (30*A*a**2*b + 6*A*b**3 + 10*B*a**3 + 18*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_565():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**3/sqrt(sec(c + d*x))
    F = 2*B*b*(a*sec(c + d*x) + b)**2*sin(c + d*x)/(9*d*sec(c + d*x)**(sympy.S(7)/2)) + 2*b**2*(9*A*b + 13*B*a)*sin(c + d*x)/(63*d*sec(c + d*x)**(sympy.S(5)/2)) + 2*b*(27*A*a*b + 22*B*a**2 + 7*B*b**2)*sin(c + d*x)/(45*d*sec(c + d*x)**(sympy.S(3)/2)) + (30*A*a**3 + 54*A*a*b**2 + 54*B*a**2*b + 14*B*b**3)*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(15*d) + (42*A*a**2*b + 10*A*b**3 + 14*B*a**3 + 30*B*a*b**2)*sin(c + d*x)/(21*d*sqrt(sec(c + d*x))) + (42*A*a**2*b + 10*A*b**3 + 14*B*a**3 + 30*B*a*b**2)*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_566():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('A') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_567():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))
    F = ((Integer(-2) * Symbol('A') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('A') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_568():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a + b*cos(c + d*x))
    F = ((Integer(2) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_569():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))*sqrt(sec(c + d*x)))
    F = ((Integer(2) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * Symbol('d')))**(Integer(-1))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_570():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * ((Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * (Symbol('a'))**(Integer(2)) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Symbol('c') + (Symbol('d') * x)) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * (Symbol('a') + Symbol('b')) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_571():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**2
    F = ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(3)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(4) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(3) * Symbol('a') * Symbol('b') * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_572():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a'))**(Integer(2)) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('a'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_573():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**2
    F = (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_574():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**2*sqrt(sec(c + d*x)))
    F = ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + ((((Symbol('a') * Symbol('A') * Symbol('b')) + ((Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_575():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * ((((Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_576():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**2*sec(c + d*x)**(sympy.S(5)/2))
    F = ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(2) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('b'))**(Integer(3)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(12) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(16) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(2) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(3) * (Symbol('b'))**(Integer(4)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(2) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_577():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(8) * (Symbol('a'))**(Integer(4)) * Symbol('A')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(4))) + (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(3)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(11) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_578():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**3
    F = (Integer(-1) * ((((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * Symbol('b') * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('b') * ((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('a'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_579():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**3*sqrt(sec(c + d*x)))
    F = ((((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(5) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(5)))) + ((Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * Symbol('a') * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + (Integer(-1) * ((((Integer(7) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(3) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('a') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_580():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * (((((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((((Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('B')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(10) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(5)))) + (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(15) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + (Integer(-1) * ((((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3))) + ((Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_581():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(5)/2))
    F = (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(29) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + ((((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(33) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(38) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(35) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * (((Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(7) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(11) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_582():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**3*sec(c + d*x)**(sympy.S(7)/2))
    F = ((((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(29) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(65) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * (Symbol('b'))**(Integer(4)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(45) * (Symbol('a'))**(Integer(5)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(99) * (Symbol('a'))**(Integer(3)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(72) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(105) * (Symbol('a'))**(Integer(6)) * Symbol('B'))) + (Integer(223) * (Symbol('a'))**(Integer(4)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(128) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(6)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(12) * (Symbol('b'))**(Integer(5)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * ((Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(38) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(35) * Symbol('A') * (Symbol('b'))**(Integer(5))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(86) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(63) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.Function('EllipticPi')(((Integer(2) * Symbol('b')) * ((Symbol('a') + Symbol('b')))**(Integer(-1))), ((Integer(2))**(Integer(-1)) * (Symbol('c') + (Symbol('d') * x))), Integer(2)) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Integer(4) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(2)) * (Symbol('b'))**(Integer(5)) * ((Symbol('a') + Symbol('b')))**(Integer(3)) * Symbol('d')))**(Integer(-1))) + (Integer(-1) * ((((Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(33) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(35) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(61) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(9) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(7) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(13) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * (Symbol('b') + (Symbol('a') * sympy.sec((Symbol('c') + (Symbol('d') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_583():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))
    F = 2*B*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_584():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))
    F = 2*B*sin(c + d*x)*sqrt(sec(c + d*x))/d - 2*B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_585():
    f = (B*a + B*b*cos(c + d*x))*sqrt(sec(c + d*x))/(a + b*cos(c + d*x))
    F = 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_586():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))*sqrt(sec(c + d*x)))
    F = 2*B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_587():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = 2*B*sin(c + d*x)/(3*d*sqrt(sec(c + d*x))) + 2*B*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_588():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2))
    F = 2*B*sin(c + d*x)/(5*d*sec(c + d*x)**(sympy.S(3)/2)) + 6*B*sqrt(cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)*sqrt(sec(c + d*x))/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_589():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + sqrt(a + b*cos(c + d*x))*(2*A*b + 14*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*a*d) + sqrt(a + b*cos(c + d*x))*(50*A*a**2 - 8*A*b**2 + 14*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(8*A*b**2 + a**2*(25*A - 63*B) + 2*a*b*(3*A - 7*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**3*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(19*A*a**2*b + 8*A*b**3 + 63*B*a**3 - 14*B*a*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**4*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_590():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + sqrt(a + b*cos(c + d*x))*(2*A*b + 10*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*a*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a + 2*A*b - 5*B*a)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**2*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a**2 - 2*A*b**2 + 5*B*a*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_591():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(A - 3*B)*sqrt(a + b)*(2*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(A*b + 3*B*a)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_592():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('A') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * (Symbol('A') + (Integer(-1) * Symbol('B')))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_593():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))*sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A')) + Symbol('B')) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_594():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))/sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + ((Symbol('a') + (Integer(2) * Symbol('b'))) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_595():
    f = (A + B*cos(c + d*x))*sqrt(a + b*cos(c + d*x))/sec(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(6) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * (Symbol('a') + (Integer(2) * Symbol('b'))) * ((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B'))) + (Integer(8) * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(-1) * (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(6) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_596():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d) + sqrt(a + b*cos(c + d*x))*(20*A*b + 18*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(63*d) + sqrt(a + b*cos(c + d*x))*(98*A*a**2 + 6*A*b**2 + 144*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(315*a*d) + sqrt(a + b*cos(c + d*x))*(176*A*a**2*b - 8*A*b**3 + 150*B*a**3 + 18*B*a*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(8*A*b**3 - a**3*(147*A - 75*B) + 3*a**2*b*(13*A - 57*B) + 6*a*b**2*(A - 3*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**3*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*A*a**4 + 33*A*a**2*b**2 + 8*A*b**4 + 246*B*a**3*b - 18*B*a*b**3)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**4*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_597():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + sqrt(a + b*cos(c + d*x))*(16*A*b + 14*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + sqrt(a + b*cos(c + d*x))*(50*A*a**2 + 6*A*b**2 + 84*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*a*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(6*A*b**2 - a**2*(25*A - 63*B) + 3*a*b*(19*A - 7*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**2*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(82*A*a**2*b - 6*A*b**3 + 63*B*a**3 + 21*B*a*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_598():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = 2*A*a*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*d) + sqrt(a + b*cos(c + d*x))*(12*A*b + 10*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a - 3*A*b - 5*B*a + 15*B*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a**2 + 3*A*b**2 + 20*B*a*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_599():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(2))) + ((Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(-1) * (Integer(3) * Symbol('B'))))) + (Integer(-1) * (Symbol('a') * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * Symbol('b') * Symbol('B'))))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_600():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * (Symbol('A') + (Integer(-1) * Symbol('B')))) + (Integer(-1) * (Symbol('b') * ((Integer(4) * Symbol('A')) + Symbol('B'))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_601():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * Symbol('a') * Symbol('A')) + (Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B')) + (Integer(2) * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(12) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_602():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(12) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(14) * Symbol('a') * Symbol('b') * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(12) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(7) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(12) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(30) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_603():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(3)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(156) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(4) * Symbol('A')) + Symbol('B')))) + (Integer(-1) * (Integer(8) * (Symbol('b'))**(Integer(3)) * ((Integer(16) * Symbol('A')) + (Integer(9) * Symbol('B'))))) + (Integer(-1) * (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(28) * Symbol('A')) + (Integer(39) * Symbol('B')))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(96) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(-1) * (Integer(24) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(48) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(8) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(12) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(32) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(24) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(9) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(156) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_604():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(13)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(11)/2)/(11*d) + 2*a*sqrt(a + b*cos(c + d*x))*(14*A*b + 11*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(99*d) + sqrt(a + b*cos(c + d*x))*(162*A*a**2 + 226*A*b**2 + 418*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(693*d) + sqrt(a + b*cos(c + d*x))*(2290*A*a**2*b + 30*A*b**3 + 1078*B*a**3 + 1650*B*a*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(3465*a*d) + sqrt(a + b*cos(c + d*x))*(1350*A*a**4 + 2050*A*a**2*b**2 - 40*A*b**4 + 3586*B*a**3*b + 110*B*a*b**3)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3465*a**2*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(40*A*b**4 + 3*a**4*(225*A - 539*B) - 6*a**3*b*(505*A - 209*B) + 15*a**2*b**2*(19*A - 121*B) + 10*a*b**3*(3*A - 11*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3465*a**3*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(3705*A*a**4*b + 255*A*a**2*b**3 + 40*A*b**5 + 1617*B*a**5 + 3069*B*a**3*b**2 - 110*B*a*b**4)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3465*a**4*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_605():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(11)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(9)/2)/(9*d) + 2*a*sqrt(a + b*cos(c + d*x))*(4*A*b + 3*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(21*d) + sqrt(a + b*cos(c + d*x))*(98*A*a**2 + 150*A*b**2 + 270*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(315*d) + sqrt(a + b*cos(c + d*x))*(326*A*a**2*b + 10*A*b**3 + 150*B*a**3 + 270*B*a*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(315*a*d) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(10*A*b**3 + 3*a**3*(49*A - 25*B) - 6*a**2*b*(19*A - 60*B) + 15*a*b**2*(11*A - 3*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**2*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(147*A*a**4 + 279*A*a**2*b**2 - 10*A*b**4 + 435*B*a**3*b + 45*B*a*b**3)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(315*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_606():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(9)/2)
    F = 2*A*a*(a + b*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(7)/2)/(7*d) + 2*a*sqrt(a + b*cos(c + d*x))*(10*A*b + 7*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(35*d) + sqrt(a + b*cos(c + d*x))*(50*A*a**2 + 90*A*b**2 + 154*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(105*d) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(a**2*(25*A - 63*B) - 8*a*b*(15*A - 7*B) + 15*b**2*(A - 7*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(145*A*a**2*b + 15*A*b**3 + 63*B*a**3 + 161*B*a*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(105*a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_607():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(7)/2)
    F = ((Integer(2) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(23) * Symbol('A') * (Symbol('b'))**(Integer(2))) + (Integer(35) * Symbol('a') * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(15) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(23) * Symbol('A')) + (Integer(-1) * (Integer(45) * Symbol('B')))))) + ((Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(17) * Symbol('A')) + (Integer(-1) * (Integer(35) * Symbol('B'))))) + (Integer(-1) * ((Symbol('a'))**(Integer(3)) * ((Integer(9) * Symbol('A')) + (Integer(-1) * (Integer(5) * Symbol('B'))))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(15) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(15) * Symbol('d')))**(Integer(-1))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_608():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(14) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('a') * Symbol('b') * ((Integer(7) * Symbol('A')) + (Integer(-1) * (Integer(9) * Symbol('B'))))) + (Integer(-1) * (Integer(2) * (Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(-1) * (Integer(3) * Symbol('B')))))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * ((Integer(6) * Symbol('A')) + Symbol('B'))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(5) * Symbol('a') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1))) + (Integer(-1) * ((((Integer(14) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(6) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_609():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2)
    F = (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(8) * (Symbol('a'))**(Integer(2)) * (Symbol('A') + (Integer(-1) * Symbol('B')))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('A')) + Symbol('B')))) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('b') * ((Integer(8) * Symbol('A')) + (Integer(3) * Symbol('B')))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(20) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('b') * ((Integer(4) * Symbol('a') * Symbol('A')) + (Integer(-1) * (Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((((Integer(8) * (Symbol('a'))**(Integer(2)) * Symbol('A')) + (Integer(-1) * (Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(2)))) + (Integer(-1) * (Integer(9) * Symbol('a') * Symbol('b') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d')))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * Symbol('A') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (Symbol('d'))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_610():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Integer(3) * Symbol('A')) + (Integer(4) * Symbol('B')))) + ((Symbol('a'))**(Integer(2)) * ((Integer(48) * Symbol('A')) + (Integer(33) * Symbol('B')))) + (Symbol('a') * ((Integer(54) * Symbol('A') * Symbol('b')) + (Integer(26) * Symbol('b') * Symbol('B'))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(30) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(8) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(20) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(8) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(3) * Symbol('a') * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(54) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(33) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(16) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_611():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/sqrt(sec(c + d*x))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(8) * (Symbol('b'))**(Integer(3)) * ((Integer(16) * Symbol('A')) + (Integer(9) * Symbol('B')))) + (Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('b') * ((Integer(132) * Symbol('A')) + (Integer(59) * Symbol('B')))) + (Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(52) * Symbol('A')) + (Integer(71) * Symbol('B'))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(192) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(40) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(160) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(120) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(48) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(64) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('b') * Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * Symbol('d') * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((((Integer(24) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('B')) + (Integer(12) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(32) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(8) * Symbol('A') * Symbol('b')) + (Integer(11) * Symbol('a') * Symbol('B'))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(24) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(264) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(128) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(284) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(192) * Symbol('b') * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_612():
    f = (A + B*cos(c + d*x))*(a + b*cos(c + d*x))**(sympy.S(5)/2)/sec(c + d*x)**(sympy.S(3)/2)
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(150) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(2840) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(45) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(1692) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(1024) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(1920) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * (((Integer(1920) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)) * (sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(45) * (Symbol('a'))**(Integer(4)) * Symbol('B')) + (Integer(-1) * (Integer(30) * (Symbol('a'))**(Integer(3)) * Symbol('b') * ((Integer(5) * Symbol('A')) + Symbol('B')))) + (Integer(-1) * (Integer(16) * (Symbol('b'))**(Integer(4)) * ((Integer(45) * Symbol('A')) + (Integer(64) * Symbol('B'))))) + (Integer(-1) * (Integer(8) * Symbol('a') * (Symbol('b'))**(Integer(3)) * ((Integer(355) * Symbol('A')) + (Integer(193) * Symbol('B'))))) + (Integer(-1) * (Integer(4) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * ((Integer(295) * Symbol('A')) + (Integer(423) * Symbol('B')))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(10) * (Symbol('a'))**(Integer(4)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(240) * (Symbol('a'))**(Integer(2)) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(96) * Symbol('A') * (Symbol('b'))**(Integer(5)))) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(5)) * Symbol('B'))) + (Integer(-1) * (Integer(40) * (Symbol('a'))**(Integer(3)) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(240) * Symbol('a') * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(128) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(50) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(120) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(172) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(320) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(50) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(64) * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(240) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(10) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(40) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('B') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(5) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(150) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(2840) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(-1) * (Integer(45) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(1692) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(1024) * (Symbol('b'))**(Integer(4)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(1920) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_613():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(7)/2)/sqrt(a + b*cos(c + d*x))
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(5)/2)/(5*a*d) - sqrt(a + b*cos(c + d*x))*(8*A*b - 10*B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(15*a**2*d) - 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(8*A*b**2 + a**2*(9*A - 5*B) - 2*a*b*(A + 5*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**3*d*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(9*A*a**2 + 8*A*b**2 - 10*B*a*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(15*a**4*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_614():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/sqrt(a + b*cos(c + d*x))
    F = 2*A*sqrt(a + b*cos(c + d*x))*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d) + 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*A*b + a*(A - 3*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*sqrt(sec(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*(2*A*b - 3*B*a)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_615():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/sqrt(a + b*cos(c + d*x))
    F = 2*A*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(a - b)*sqrt(a + b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(sec(c + d*x))) - 2*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(A - B)*sqrt(a + b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_616():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/sqrt(a + b*cos(c + d*x))
    F = (((Symbol('a') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)) * (Integer(2) * Symbol('A') * sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))) + (Integer(-1) * (((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)) * (Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1)))))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_617():
    f = (A + B*cos(c + d*x))/(sqrt(a + b*cos(c + d*x))*sqrt(sec(c + d*x)))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * Symbol('B') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_618():
    f = (A + B*cos(c + d*x))/(sqrt(a + b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B'))) + (Integer(2) * Symbol('b') * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(4) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + (Integer(-1) * (Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('B') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(2) * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(4) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * Symbol('d')))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_619():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a + b*cos(c + d*x))*(2*A*a**2 - 8*A*b**2 + 6*B*a*b)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*(a**2 - b**2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a + 4*b)*(4*A*b + a*(A - 3*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(a + b)*sqrt(sec(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(10*A*a**2*b - 16*A*b**3 - 6*B*a**3 + 12*B*a*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*sqrt(a + b)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_620():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(4*A*b + 2*a*(A - B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*A*a**2 - 4*A*b**2 + 2*B*a*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**3*d*sqrt(a + b)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_621():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = -(2*A*b - 2*B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*A + 2*B)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(a + b)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*A*b - 2*B*a)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(a + b)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_622():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = (Integer(-1) * ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_623():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (((Integer(3) * Symbol('a')) + Symbol('b')) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * Symbol('a') * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(3) * (Symbol('a'))**(Integer(2)) * Symbol('B'))) + ((Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * (((Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_624():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(5)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*b*(10*A*a**2*b - 6*A*b**3 - 7*B*a**3 + 3*B*a*b**2)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt(a + b*cos(c + d*x))*(2*A*a**4 - 26*A*a**2*b**2 + 16*A*b**4 + 16*B*a**3*b - 8*B*a*b**3)*sin(c + d*x)*sec(c + d*x)**(sympy.S(3)/2)/(3*a**3*d*(a**2 - b**2)**2) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(32*A*b**4 - 2*a**4*(A - 3*B) - 18*a**3*b*(A - B) - 4*a**2*b**2*(8*A + 3*B) + 8*a*b**3*(3*A - 2*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*sqrt(a + b)*(a**2 - b**2)*sqrt(sec(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(16*A*a**4*b - 56*A*a**2*b**3 + 32*A*b**5 - 6*B*a**5 + 30*B*a**3*b**2 - 16*B*a*b**4)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**5*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_625():
    f = (A + B*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)) + 2*b*(8*A*a**2*b - 4*A*b**3 - 5*B*a**3 + B*a*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a**2*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(16*A*b**3 - 6*a**3*(A - B) - 6*a**2*b*(3*A + B) + 4*a*b**2*(3*A - B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*sqrt(a + b)*(a**2 - b**2)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*A*a**4 - 30*A*a**2*b**2 + 16*A*b**4 + 12*B*a**3*b - 4*B*a*b**3)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**4*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_626():
    f = (A + B*cos(c + d*x))*sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*b*(A*b - B*a)*sin(c + d*x)/(3*a*d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(a**2 - b**2)*sqrt(sec(c + d*x))) - (12*A*a**2*b - 4*A*b**3 - 6*B*a**3 - 2*B*a*b**2)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*a*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(4*A*b**2 - 6*a**2*(A + B) + 2*a*b*(3*A + B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*sqrt(a + b)*(a**2 - b**2)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(12*A*a**2*b - 4*A*b**3 - 6*B*a**3 - 2*B*a*b**2)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**3*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_627():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sqrt(sec(c + d*x)))
    F = (6*A*a**2 + 2*A*b**2 - 8*B*a*b)*sin(c + d*x)*sqrt(sec(c + d*x))/(3*d*sqrt(a + b*cos(c + d*x))*(a**2 - b**2)**2) - (2*A*b - 2*B*a)*sin(c + d*x)/(d*(a + b*cos(c + d*x))**(sympy.S(3)/2)*(3*a**2 - 3*b**2)*sqrt(sec(c + d*x))) + sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(2*a*(3*A + B) - 2*b*(A + 3*B))*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x))) - sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*(6*A*a**2 + 2*A*b**2 - 8*B*a*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(3*a**2*d*(a - b)*(a + b)**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_628():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = ((Integer(2) * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Integer(3) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + ((Symbol('a'))**(Integer(2)) * Symbol('b') * Symbol('B')) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * (Symbol('A') + (Integer(6) * Symbol('B')))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(2)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(3)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(4) * Symbol('A') * (Symbol('b'))**(Integer(3))) + (Integer(3) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B')))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_629():
    f = (A + B*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**(sympy.S(5)/2))
    F = ((((Integer(6) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * Symbol('a') * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((((Integer(3) * (Symbol('b'))**(Integer(3)) * ((Integer(4) * Symbol('A')) + (Integer(-1) * Symbol('B')))) + (Integer(15) * (Symbol('a'))**(Integer(3)) * Symbol('B')) + (Integer(-1) * (Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(2) * Symbol('A')) + (Integer(21) * Symbol('B'))))) + (Integer(-1) * ((Symbol('a'))**(Integer(2)) * ((Integer(6) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('b') * Symbol('B'))))))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Integer(3) * (Symbol('a') + (Integer(-1) * Symbol('b'))) * (Symbol('b'))**(Integer(3)) * ((Symbol('a') + Symbol('b')))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * ((Integer(2) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(5) * Symbol('a') * Symbol('B')))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(4)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('d') * ((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (sympy.sec((Symbol('c') + (Symbol('d') * x))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * Symbol('a') * ((Integer(2) * (Symbol('a'))**(Integer(2)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(6) * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(5) * (Symbol('a'))**(Integer(3)) * Symbol('B'))) + (Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(2)) * Symbol('B'))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(6) * (Symbol('a'))**(Integer(3)) * Symbol('A') * Symbol('b')) + (Integer(-1) * (Integer(14) * Symbol('a') * Symbol('A') * (Symbol('b'))**(Integer(3)))) + (Integer(-1) * (Integer(15) * (Symbol('a'))**(Integer(4)) * Symbol('B'))) + (Integer(26) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * Symbol('B')) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(4)) * Symbol('B')))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Integer(3) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))))**(Integer(2)) * Symbol('d')))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_630():
    f = (B*a + B*b*cos(c + d*x))*sec(c + d*x)**(sympy.S(3)/2)/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = -2*B*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(sec(c + d*x))) + B*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*(2*a - 2*b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_e(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a**2*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_631():
    f = (B*a + B*b*cos(c + d*x))*sqrt(sec(c + d*x))/(a + b*cos(c + d*x))**(sympy.S(3)/2)
    F = 2*B*sqrt(a*(1 - sec(c + d*x))/(a + b))*sqrt(a*(sec(c + d*x) + 1)/(a - b))*sqrt(a + b)*sqrt(cos(c + d*x))*csc(c + d*x)*elliptic_f(asin(sqrt(a + b*cos(c + d*x))/(sqrt(a + b)*sqrt(cos(c + d*x)))), -(a + b)/(a - b))/(a*d*sqrt(sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_632():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sqrt(sec(c + d*x)))
    F = Integer(-1) * ((Integer(2) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_633():
    f = (B*a + B*b*cos(c + d*x))/((a + b*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**(sympy.S(3)/2))
    F = (Integer(-1) * (((Symbol('a') + (Integer(-1) * Symbol('b'))) * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_e(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('a') * Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))) + ((sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.elliptic_f(sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * ((Symbol('b') * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * sympy.sqrt((Symbol('a') + Symbol('b'))) * Symbol('B') * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x)))) * sympy.csc((Symbol('c') + (Symbol('d') * x))) * sympy.Function('EllipticPi')(((Symbol('a') + Symbol('b')) * (Symbol('b'))**(Integer(-1))), sympy.asin((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * ((sympy.sqrt((Symbol('a') + Symbol('b'))) * sympy.sqrt(sympy.cos((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1)))), (Integer(-1) * ((Symbol('a') + Symbol('b')) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * sympy.sqrt(((Symbol('a') * (Integer(1) + (Integer(-1) * sympy.sec((Symbol('c') + (Symbol('d') * x)))))) * ((Symbol('a') + Symbol('b')))**(Integer(-1)))) * sympy.sqrt(((Symbol('a') * (Integer(1) + sympy.sec((Symbol('c') + (Symbol('d') * x))))) * ((Symbol('a') + (Integer(-1) * Symbol('b'))))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('d') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('B') * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x)))))) * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x))))))**(Integer(-1))) + ((Symbol('a') * Symbol('B') * sympy.sqrt(sympy.sec((Symbol('c') + (Symbol('d') * x)))) * sympy.sin((Symbol('c') + (Symbol('d') * x)))) * ((Symbol('b') * Symbol('d') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('c') + (Symbol('d') * x))))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_634():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**n
    F = ((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.Function('Unintegrable')(((((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))**(Symbol('n')) * (Symbol('A') + (Symbol('B') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_635():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**4
    F = -A*a*c**5*(c*sec(e + f*x))**(m - 5)*(a*sec(e + f*x) + b)**3*tan(e + f*x)/(f*(1 - m)) - a**2*c**5*(c*sec(e + f*x))**(m - 5)*(A*a**2*(2 - m)**2 + A*b**2*(m**2 - m + 6) + 2*B*a*b*(1 - m)**2)*tan(e + f*x)*sec(e + f*x)/(f*(1 - m)*(2 - m)*(3 - m)) - a*c**5*(c*sec(e + f*x))**(m - 5)*(a*sec(e + f*x) + b)**2*(-A*b*(m + 2) + B*a*(1 - m))*tan(e + f*x)/(f*(1 - m)*(2 - m)) - a*c**5*(c*sec(e + f*x))**(m - 5)*(4*A*a**2*b*(m**2 - 4*m + 3) + 2*A*b**3*(m**2 - 2*m + 4) + B*a**3*(m**2 - 4*m + 3) + B*a*b**2*(5*m**2 - 13*m + 8))*tan(e + f*x)/(f*(1 - m)*(2 - m)*(4 - m)) - c**6*(c*sec(e + f*x))**(m - 6)*(4*A*a**3*b*(m**2 - 8*m + 15) + 4*A*a*b**3*(m**2 - 7*m + 10) + B*a**4*(m**2 - 8*m + 15) + 6*B*a**2*b**2*(m**2 - 7*m + 10) + B*b**4*(m**2 - 6*m + 8))*sin(e + f*x)*hyper((sympy.S.Half, 3 - m/2), (4 - m/2,), cos(e + f*x)**2)/(f*(2 - m)*(4 - m)*(6 - m)*sqrt(sin(e + f*x)**2)) - c**5*(c*sec(e + f*x))**(m - 5)*(A*a**4*(m**2 - 6*m + 8) + 6*A*a**2*b**2*(m**2 - 5*m + 4) + A*b**4*(m**2 - 4*m + 3) + 4*B*a**3*b*(m**2 - 5*m + 4) + 4*B*a*b**3*(m**2 - 4*m + 3))*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(5)/2 - m/2), (sympy.S(7)/2 - m/2,), cos(e + f*x)**2)/(f*(1 - m)*(3 - m)*(5 - m)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_636():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**3
    F = -A*a*c**4*(c*sec(e + f*x))**(m - 4)*(a*sec(e + f*x) + b)**2*tan(e + f*x)/(f*(1 - m)) - a**2*c**4*(c*sec(e + f*x))**(m - 4)*(-A*b*(m + 1) + B*a*(1 - m))*tan(e + f*x)*sec(e + f*x)/(f*(1 - m)*(2 - m)) - a*c**4*(c*sec(e + f*x))**(m - 4)*(A*a**2*(2 - m) - 2*A*b**2*m + 3*B*a*b*(1 - m))*tan(e + f*x)/(f*(1 - m)*(3 - m)) - c**5*(c*sec(e + f*x))**(m - 5)*(A*a**3*(m**2 - 6*m + 8) + 3*A*a*b**2*(m**2 - 5*m + 4) + 3*B*a**2*b*(m**2 - 5*m + 4) + B*b**3*(m**2 - 4*m + 3))*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(5)/2 - m/2), (sympy.S(7)/2 - m/2,), cos(e + f*x)**2)/(f*(1 - m)*(3 - m)*(5 - m)*sqrt(sin(e + f*x)**2)) - c**4*(c*sec(e + f*x))**(m - 4)*(3*A*a**2*b*(3 - m) + A*b**3*(2 - m) + B*a**3*(3 - m) + 3*B*a*b**2*(2 - m))*sin(e + f*x)*hyper((sympy.S.Half, 2 - m/2), (3 - m/2,), cos(e + f*x)**2)/(f*(2 - m)*(4 - m)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_637():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**2
    F = -A*a*c**3*(c*sec(e + f*x))**(m - 3)*(a*sec(e + f*x) + b)*tan(e + f*x)/(f*(1 - m)) - a*c**3*(c*sec(e + f*x))**(m - 3)*(-A*b*m + B*a*(1 - m))*tan(e + f*x)/(f*(1 - m)*(2 - m)) - c**4*(c*sec(e + f*x))**(m - 4)*(2*A*a*b*(3 - m) + B*a**2*(3 - m) + B*b**2*(2 - m))*sin(e + f*x)*hyper((sympy.S.Half, 2 - m/2), (3 - m/2,), cos(e + f*x)**2)/(f*(2 - m)*(4 - m)*sqrt(sin(e + f*x)**2)) - c**3*(c*sec(e + f*x))**(m - 3)*(A*a**2*(2 - m) + A*b**2*(1 - m) + 2*B*a*b*(1 - m))*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(3)/2 - m/2), (sympy.S(5)/2 - m/2,), cos(e + f*x)**2)/(f*(1 - m)*(3 - m)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_638():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))
    F = -A*a*c**2*(c*sec(e + f*x))**(m - 2)*tan(e + f*x)/(f*(1 - m)) - c**3*(c*sec(e + f*x))**(m - 3)*(A*a*(2 - m) + B*b*(1 - m))*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(3)/2 - m/2), (sympy.S(5)/2 - m/2,), cos(e + f*x)**2)/(f*(1 - m)*(3 - m)*sqrt(sin(e + f*x)**2)) - c**2*(c*sec(e + f*x))**(m - 2)*(A*b + B*a)*sin(e + f*x)*hyper((sympy.S.Half, 1 - m/2), (2 - m/2,), cos(e + f*x)**2)/(f*(2 - m)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_639():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))/(a + b*cos(e + f*x))
    F = -B*c*(c*sec(e + f*x))**(m - 1)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(e + f*x)**2)/(b*f*(1 - m)*sqrt(sin(e + f*x)**2)) + a*(c*sec(e + f*x))**(m + 1)*(A*b - B*a)*(cos(e + f*x)**2)**(m/2 + sympy.S.Half)*sin(e + f*x)*appellf1(sympy.S.Half, 1, m/2 + sympy.S.Half, sympy.S(3)/2, -b**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(b*c*f*(a**2 - b**2)) - (c*sec(e + f*x))**(m + 1)*(A*b - B*a)*(cos(e + f*x)**2)**(m/2)*sin(e + f*x)*cos(e + f*x)*appellf1(sympy.S.Half, 1, m/2, sympy.S(3)/2, -b**2*sin(e + f*x)**2/(a**2 - b**2), sin(e + f*x)**2)/(c*f*(a**2 - b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_640():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))*(a + b*cos(e + f*x))**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('b') * Symbol('B') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('f') * (Integer(5) + (Integer(-1) * (Integer(2) * Symbol('m'))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.Function('Unintegrable')(((((Integer(2))**(Integer(-1)) * Symbol('a') * Symbol('c') * ((Integer(2) * Symbol('b') * Symbol('B') * (Integer(1) + (Integer(-1) * Symbol('m')))) + (Integer(2) * Symbol('a') * Symbol('A') * ((Integer(5) * (Integer(2))**(Integer(-1))) + (Integer(-1) * Symbol('m')))))) + ((Integer(2))**(Integer(-1)) * Symbol('c') * (((Symbol('b'))**(Integer(2)) * Symbol('B') * (Integer(3) + (Integer(-1) * (Integer(2) * Symbol('m'))))) + (Symbol('a') * ((Integer(2) * Symbol('A') * Symbol('b')) + (Symbol('a') * Symbol('B'))) * (Integer(5) + (Integer(-1) * (Integer(2) * Symbol('m')))))) * sympy.cos((Symbol('e') + (Symbol('f') * x)))) + ((Integer(2))**(Integer(-1)) * Symbol('b') * Symbol('c') * ((Symbol('A') * Symbol('b') * (Integer(5) + (Integer(-1) * (Integer(2) * Symbol('m'))))) + (Integer(2) * Symbol('a') * Symbol('B') * (Integer(3) + (Integer(-1) * Symbol('m'))))) * (sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))) * ((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))), x)) * ((Symbol('c') * (Integer(5) + (Integer(-1) * (Integer(2) * Symbol('m'))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_641():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))*sqrt(a + b*cos(e + f*x))
    F = ((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.Function('Unintegrable')(((sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (Symbol('A') + (Symbol('B') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_642():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))/sqrt(a + b*cos(e + f*x))
    F = ((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.Function('Unintegrable')(((Symbol('A') + (Symbol('B') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) * ((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))), x)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_2_Cosine_4_2_3_1_a_plus_b_cos_pow_m_c_plus_d_cos_pow_n_A_plus_B_cos_643():
    f = (c*sec(e + f*x))**m*(A + B*cos(e + f*x))/(a + b*cos(e + f*x))**(sympy.S(3)/2)
    F = ((Integer(2) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('f') * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))) + ((Integer(2) * ((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * ((Symbol('c') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.Function('Unintegrable')(((((Integer(2))**(Integer(-1)) * Symbol('c') * (((Symbol('a'))**(Integer(2)) * Symbol('A')) + (Symbol('A') * (Symbol('b'))**(Integer(2)) * (Integer(1) + (Integer(-1) * (Integer(2) * Symbol('m'))))) + (Integer(-1) * (Integer(2) * Symbol('a') * Symbol('b') * Symbol('B') * (Integer(1) + (Integer(-1) * Symbol('m'))))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('a') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x))))) + (Integer(-1) * ((Integer(2))**(Integer(-1)) * Symbol('b') * ((Symbol('A') * Symbol('b')) + (Integer(-1) * (Symbol('a') * Symbol('B')))) * Symbol('c') * (Integer(3) + (Integer(-1) * (Integer(2) * Symbol('m')))) * (sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((((Symbol('c') * sympy.cos((Symbol('e') + (Symbol('f') * x)))))**(Symbol('m')) * sympy.sqrt((Symbol('a') + (Symbol('b') * sympy.cos((Symbol('e') + (Symbol('f') * x))))))))**(Integer(-1))), x)) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Symbol('b'))**(Integer(2)))) * Symbol('c')))**(Integer(-1)))
    assert integrate(f, x) == F

