"""Generated from MathematicaSyntaxTestSuite.

Source: 4 Trig functions/4.3 Tangent/4.3.1.2 (d sec)^m (a+b tan)^n.m
"""

from sympy import *
from sympy.testing.pytest import XFAIL

import sympy



x = symbols('x')

a, b, c, d, e, f, m, n = symbols('a b c d e f m n')

def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_1():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)**10
    F = a*tan(c + d*x)**9/(9*d) + 4*a*tan(c + d*x)**7/(7*d) + 6*a*tan(c + d*x)**5/(5*d) + 4*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + I*a*sec(c + d*x)**10/(10*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_2():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)**8
    F = a*tan(c + d*x)**7/(7*d) + 3*a*tan(c + d*x)**5/(5*d) + a*tan(c + d*x)**3/d + a*tan(c + d*x)/d + I*a*sec(c + d*x)**8/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_3():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)**6
    F = a*tan(c + d*x)**5/(5*d) + 2*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + I*a*sec(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_4():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + I*a*sec(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_5():
    f = I*a*tan(c + d*x) + a
    F = a*x - I*a*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_6():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)**2
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) - I*a*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_7():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)**4
    F = 3*a*x/8 + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - I*a*cos(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_8():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)**6
    F = 5*a*x/16 + a*sin(c + d*x)*cos(c + d*x)**5/(6*d) + 5*a*sin(c + d*x)*cos(c + d*x)**3/(24*d) + 5*a*sin(c + d*x)*cos(c + d*x)/(16*d) - I*a*cos(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_9():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)**8
    F = 35*a*x/128 + a*sin(c + d*x)*cos(c + d*x)**7/(8*d) + 7*a*sin(c + d*x)*cos(c + d*x)**5/(48*d) + 35*a*sin(c + d*x)*cos(c + d*x)**3/(192*d) + 35*a*sin(c + d*x)*cos(c + d*x)/(128*d) - I*a*cos(c + d*x)**8/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_10():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)**7
    F = a*tan(c + d*x)*sec(c + d*x)**5/(6*d) + 5*a*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 5*a*tan(c + d*x)*sec(c + d*x)/(16*d) + 5*a*atanh(sin(c + d*x))/(16*d) + I*a*sec(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_11():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)**5
    F = a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a*atanh(sin(c + d*x))/(8*d) + I*a*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_12():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)**3
    F = a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*atanh(sin(c + d*x))/(2*d) + I*a*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_13():
    f = (I*a*tan(c + d*x) + a)*sec(c + d*x)
    F = a*atanh(sin(c + d*x))/d + I*a*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_14():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)
    F = a*sin(c + d*x)/d - I*a*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_15():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)**3
    F = -a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d - I*a*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_16():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)**5
    F = a*sin(c + d*x)**5/(5*d) - 2*a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d - I*a*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_17():
    f = (I*a*tan(c + d*x) + a)*cos(c + d*x)**7
    F = -a*sin(c + d*x)**7/(7*d) + 3*a*sin(c + d*x)**5/(5*d) - a*sin(c + d*x)**3/d + a*sin(c + d*x)/d - I*a*cos(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_18():
    f = (I*a*tan(c + d*x) + a)**2*sec(c + d*x)**8
    F = -4*I*(I*a*tan(c + d*x) + a)**6/(3*a**4*d) + 12*I*(I*a*tan(c + d*x) + a)**7/(7*a**5*d) - 3*I*(I*a*tan(c + d*x) + a)**8/(4*a**6*d) + I*(I*a*tan(c + d*x) + a)**9/(9*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_19():
    f = (I*a*tan(c + d*x) + a)**2*sec(c + d*x)**6
    F = -4*I*(I*a*tan(c + d*x) + a)**5/(5*a**3*d) + 2*I*(I*a*tan(c + d*x) + a)**6/(3*a**4*d) - I*(I*a*tan(c + d*x) + a)**7/(7*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_20():
    f = (I*a*tan(c + d*x) + a)**2*sec(c + d*x)**4
    F = -I*(I*a*tan(c + d*x) + a)**4/(2*a**2*d) + I*(I*a*tan(c + d*x) + a)**5/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_21():
    f = (I*a*tan(c + d*x) + a)**2*sec(c + d*x)**2
    F = -I*(I*a*tan(c + d*x) + a)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_22():
    f = (I*a*tan(c + d*x) + a)**2
    F = 2*a**2*x - 2*I*a**2*log(cos(c + d*x))/d - a**2*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_23():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**2
    F = -I*a**3/(d*(-I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_24():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**4
    F = -I*a**4/(4*d*(-I*a*tan(c + d*x) + a)**2) - I*a**3/(4*d*(-I*a*tan(c + d*x) + a)) + a**2*x/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_25():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**6
    F = -I*a**5/(12*d*(-I*a*tan(c + d*x) + a)**3) - I*a**4/(8*d*(-I*a*tan(c + d*x) + a)**2) + I*a**3/(16*d*(I*a*tan(c + d*x) + a)) - 3*I*a**3/(16*d*(-I*a*tan(c + d*x) + a)) + a**2*x/4
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_26():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**8
    F = -I*a**6/(32*d*(-I*a*tan(c + d*x) + a)**4) - I*a**5/(16*d*(-I*a*tan(c + d*x) + a)**3) + I*a**4/(64*d*(I*a*tan(c + d*x) + a)**2) - 3*I*a**4/(32*d*(-I*a*tan(c + d*x) + a)**2) + 5*I*a**3/(64*d*(I*a*tan(c + d*x) + a)) - 5*I*a**3/(32*d*(-I*a*tan(c + d*x) + a)) + 15*a**2*x/64
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_27():
    f = (I*a*tan(c + d*x) + a)**2*sec(c + d*x)**5
    F = 7*a**2*tan(c + d*x)*sec(c + d*x)**3/(24*d) + 7*a**2*tan(c + d*x)*sec(c + d*x)/(16*d) + 7*a**2*atanh(sin(c + d*x))/(16*d) + 7*I*a**2*sec(c + d*x)**5/(30*d) + I*(I*a**2*tan(c + d*x) + a**2)*sec(c + d*x)**5/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_28():
    f = (I*a*tan(c + d*x) + a)**2*sec(c + d*x)**3
    F = 5*a**2*tan(c + d*x)*sec(c + d*x)/(8*d) + 5*a**2*atanh(sin(c + d*x))/(8*d) + 5*I*a**2*sec(c + d*x)**3/(12*d) + I*(I*a**2*tan(c + d*x) + a**2)*sec(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_29():
    f = (I*a*tan(c + d*x) + a)**2*sec(c + d*x)
    F = 3*a**2*atanh(sin(c + d*x))/(2*d) + 3*I*a**2*sec(c + d*x)/(2*d) + I*(I*a**2*tan(c + d*x) + a**2)*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_30():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)
    F = -a**2*atanh(sin(c + d*x))/d - 2*I*(I*a**2*tan(c + d*x) + a**2)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_31():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**3
    F = a**2*sin(c + d*x)/(3*d) - 2*I*(I*a**2*tan(c + d*x) + a**2)*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_32():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**5
    F = -a**2*sin(c + d*x)**3/(5*d) + 3*a**2*sin(c + d*x)/(5*d) - 2*I*(I*a**2*tan(c + d*x) + a**2)*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_33():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**7
    F = a**2*sin(c + d*x)**5/(7*d) - 10*a**2*sin(c + d*x)**3/(21*d) + 5*a**2*sin(c + d*x)/(7*d) - 2*I*(I*a**2*tan(c + d*x) + a**2)*cos(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_34():
    f = (I*a*tan(c + d*x) + a)**2*cos(c + d*x)**9
    F = -a**2*sin(c + d*x)**7/(9*d) + 7*a**2*sin(c + d*x)**5/(15*d) - 7*a**2*sin(c + d*x)**3/(9*d) + 7*a**2*sin(c + d*x)/(9*d) - 2*I*(I*a**2*tan(c + d*x) + a**2)*cos(c + d*x)**9/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_35():
    f = (I*a*tan(c + d*x) + a)**3*sec(c + d*x)**8
    F = -8*I*(I*a*tan(c + d*x) + a)**7/(7*a**4*d) + 3*I*(I*a*tan(c + d*x) + a)**8/(2*a**5*d) - 2*I*(I*a*tan(c + d*x) + a)**9/(3*a**6*d) + I*(I*a*tan(c + d*x) + a)**10/(10*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_36():
    f = (I*a*tan(c + d*x) + a)**3*sec(c + d*x)**6
    F = -2*I*(I*a*tan(c + d*x) + a)**6/(3*a**3*d) + 4*I*(I*a*tan(c + d*x) + a)**7/(7*a**4*d) - I*(I*a*tan(c + d*x) + a)**8/(8*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_37():
    f = (I*a*tan(c + d*x) + a)**3*sec(c + d*x)**4
    F = -2*I*(I*a*tan(c + d*x) + a)**5/(5*a**2*d) + I*(I*a*tan(c + d*x) + a)**6/(6*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_38():
    f = (I*a*tan(c + d*x) + a)**3*sec(c + d*x)**2
    F = -I*(I*a*tan(c + d*x) + a)**4/(4*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_39():
    f = (I*a*tan(c + d*x) + a)**3
    F = 4*a**3*x - 4*I*a**3*log(cos(c + d*x))/d - 2*a**3*tan(c + d*x)/d + I*a*(I*a*tan(c + d*x) + a)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_40():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**2
    F = -2*I*a**4/(d*(-I*a*tan(c + d*x) + a)) - a**3*x + I*a**3*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_41():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**4
    F = -I*a**5/(2*d*(-I*a*tan(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_42():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**6
    F = -I*a**6/(6*d*(-I*a*tan(c + d*x) + a)**3) - I*a**5/(8*d*(-I*a*tan(c + d*x) + a)**2) - I*a**4/(8*d*(-I*a*tan(c + d*x) + a)) + a**3*x/8
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_43():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**8
    F = -I*a**7/(16*d*(-I*a*tan(c + d*x) + a)**4) - I*a**6/(12*d*(-I*a*tan(c + d*x) + a)**3) - 3*I*a**5/(32*d*(-I*a*tan(c + d*x) + a)**2) + I*a**4/(32*d*(I*a*tan(c + d*x) + a)) - I*a**4/(8*d*(-I*a*tan(c + d*x) + a)) + 5*a**3*x/32
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_44():
    f = (I*a*tan(c + d*x) + a)**3*sec(c + d*x)**3
    F = 7*a**3*tan(c + d*x)*sec(c + d*x)/(8*d) + 7*a**3*atanh(sin(c + d*x))/(8*d) + 7*I*a**3*sec(c + d*x)**3/(12*d) + I*a*(I*a*tan(c + d*x) + a)**2*sec(c + d*x)**3/(5*d) + 7*I*(I*a**3*tan(c + d*x) + a**3)*sec(c + d*x)**3/(20*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_45():
    f = (I*a*tan(c + d*x) + a)**3*sec(c + d*x)
    F = 5*a**3*atanh(sin(c + d*x))/(2*d) + 5*I*a**3*sec(c + d*x)/(2*d) + I*a*(I*a*tan(c + d*x) + a)**2*sec(c + d*x)/(3*d) + 5*I*(I*a**3*tan(c + d*x) + a**3)*sec(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_46():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)
    F = -3*a**3*atanh(sin(c + d*x))/d - 3*I*a**3*sec(c + d*x)/d - 2*I*a*(I*a*tan(c + d*x) + a)**2*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_47():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**3
    F = -I*(I*a*tan(c + d*x) + a)**3*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_48():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**5
    F = -a**3*sin(c + d*x)**3/(15*d) + a**3*sin(c + d*x)/(5*d) - I*a**3*cos(c + d*x)**3/(15*d) - 2*I*a*(I*a*tan(c + d*x) + a)**2*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_49():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**7
    F = 3*a**3*sin(c + d*x)**5/(35*d) - 2*a**3*sin(c + d*x)**3/(7*d) + 3*a**3*sin(c + d*x)/(7*d) - 3*I*a**3*cos(c + d*x)**5/(35*d) - 2*I*a*(I*a*tan(c + d*x) + a)**2*cos(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_50():
    f = (I*a*tan(c + d*x) + a)**3*cos(c + d*x)**9
    F = -5*a**3*sin(c + d*x)**7/(63*d) + a**3*sin(c + d*x)**5/(3*d) - 5*a**3*sin(c + d*x)**3/(9*d) + 5*a**3*sin(c + d*x)/(9*d) - 5*I*a**3*cos(c + d*x)**7/(63*d) - 2*I*a*(I*a*tan(c + d*x) + a)**2*cos(c + d*x)**9/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_51():
    f = (I*a*tan(c + d*x) + a)**4*sec(c + d*x)**3
    F = 21*a**4*tan(c + d*x)*sec(c + d*x)/(16*d) + 21*a**4*atanh(sin(c + d*x))/(16*d) + 7*I*a**4*sec(c + d*x)**3/(8*d) + I*a*(I*a*tan(c + d*x) + a)**3*sec(c + d*x)**3/(6*d) + 3*I*(I*a**2*tan(c + d*x) + a**2)**2*sec(c + d*x)**3/(10*d) + 21*I*(I*a**4*tan(c + d*x) + a**4)*sec(c + d*x)**3/(40*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_52():
    f = (I*a*tan(c + d*x) + a)**4*sec(c + d*x)
    F = 35*a**4*atanh(sin(c + d*x))/(8*d) + 35*I*a**4*sec(c + d*x)/(8*d) + I*a*(I*a*tan(c + d*x) + a)**3*sec(c + d*x)/(4*d) + 7*I*(I*a**2*tan(c + d*x) + a**2)**2*sec(c + d*x)/(12*d) + 35*I*(I*a**4*tan(c + d*x) + a**4)*sec(c + d*x)/(24*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_53():
    f = (I*a*tan(c + d*x) + a)**4*cos(c + d*x)
    F = -15*a**4*atanh(sin(c + d*x))/(2*d) - 15*I*a**4*sec(c + d*x)/(2*d) - 2*I*a*(I*a*tan(c + d*x) + a)**3*cos(c + d*x)/d - 5*I*(I*a**4*tan(c + d*x) + a**4)*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_54():
    f = (I*a*tan(c + d*x) + a)**4*cos(c + d*x)**3
    F = a**4*atanh(sin(c + d*x))/d - 2*I*a*(I*a*tan(c + d*x) + a)**3*cos(c + d*x)**3/(3*d) + 2*I*(I*a**4*tan(c + d*x) + a**4)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_55():
    f = (I*a*tan(c + d*x) + a)**4*cos(c + d*x)**5
    F = -I*a*(I*a*tan(c + d*x) + a)**3*cos(c + d*x)**3/(15*d) - I*(I*a*tan(c + d*x) + a)**4*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_56():
    f = (I*a*tan(c + d*x) + a)**4*cos(c + d*x)**7
    F = -a**4*sin(c + d*x)**3/(35*d) + 3*a**4*sin(c + d*x)/(35*d) - 2*I*a*(I*a*tan(c + d*x) + a)**3*cos(c + d*x)**7/(7*d) - 2*I*(I*a**4*tan(c + d*x) + a**4)*cos(c + d*x)**5/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_57():
    f = (I*a*tan(c + d*x) + a)**4*cos(c + d*x)**9
    F = a**4*sin(c + d*x)**5/(21*d) - 10*a**4*sin(c + d*x)**3/(63*d) + 5*a**4*sin(c + d*x)/(21*d) - 2*I*a*(I*a*tan(c + d*x) + a)**3*cos(c + d*x)**9/(9*d) - 2*I*(I*a**4*tan(c + d*x) + a**4)*cos(c + d*x)**7/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_58():
    f = (I*a*tan(c + d*x) + a)**5*sec(c + d*x)**8
    F = -8*I*(I*a*tan(c + d*x) + a)**9/(9*a**4*d) + 6*I*(I*a*tan(c + d*x) + a)**10/(5*a**5*d) - 6*I*(I*a*tan(c + d*x) + a)**11/(11*a**6*d) + I*(I*a*tan(c + d*x) + a)**12/(12*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_59():
    f = (I*a*tan(c + d*x) + a)**5*sec(c + d*x)**6
    F = -I*(I*a*tan(c + d*x) + a)**8/(2*a**3*d) + 4*I*(I*a*tan(c + d*x) + a)**9/(9*a**4*d) - I*(I*a*tan(c + d*x) + a)**10/(10*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_60():
    f = (I*a*tan(c + d*x) + a)**5*sec(c + d*x)**4
    F = -2*I*(I*a*tan(c + d*x) + a)**7/(7*a**2*d) + I*(I*a*tan(c + d*x) + a)**8/(8*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_61():
    f = (I*a*tan(c + d*x) + a)**5*sec(c + d*x)**2
    F = -I*(I*a*tan(c + d*x) + a)**6/(6*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_62():
    f = (I*a*tan(c + d*x) + a)**5
    F = 16*a**5*x - 16*I*a**5*log(cos(c + d*x))/d - 8*a**5*tan(c + d*x)/d + 2*I*a**2*(I*a*tan(c + d*x) + a)**3/(3*d) + I*a*(I*a*tan(c + d*x) + a)**4/(4*d) + 2*I*a*(I*a**2*tan(c + d*x) + a**2)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_63():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**2
    F = -8*I*a**6/(d*(-I*a*tan(c + d*x) + a)) - 12*a**5*x + 12*I*a**5*log(cos(c + d*x))/d + I*a**5*tan(c + d*x)**2/(2*d) + 5*a**5*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_64():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**4
    F = -2*I*a**7/(d*(-I*a*tan(c + d*x) + a)**2) + 4*I*a**6/(d*(-I*a*tan(c + d*x) + a)) + a**5*x - I*a**5*log(cos(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_65():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**6
    F = -2*I*a**8/(3*d*(-I*a*tan(c + d*x) + a)**3) + I*a**7/(2*d*(-I*a*tan(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_66():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**8
    F = -I*a**9/(4*d*(-I*a*tan(c + d*x) + a)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_67():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**10
    F = -I*a**10/(10*d*(-I*a*tan(c + d*x) + a)**5) - I*a**9/(16*d*(-I*a*tan(c + d*x) + a)**4) - I*a**8/(24*d*(-I*a*tan(c + d*x) + a)**3) - I*a**7/(32*d*(-I*a*tan(c + d*x) + a)**2) - I*a**6/(32*d*(-I*a*tan(c + d*x) + a)) + a**5*x/32
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_68():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**12
    F = -I*a**11/(24*d*(-I*a*tan(c + d*x) + a)**6) - I*a**10/(20*d*(-I*a*tan(c + d*x) + a)**5) - 3*I*a**9/(64*d*(-I*a*tan(c + d*x) + a)**4) - I*a**8/(24*d*(-I*a*tan(c + d*x) + a)**3) - 5*I*a**7/(128*d*(-I*a*tan(c + d*x) + a)**2) + I*a**6/(128*d*(I*a*tan(c + d*x) + a)) - 3*I*a**6/(64*d*(-I*a*tan(c + d*x) + a)) + 7*a**5*x/128
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_69():
    f = (I*a*tan(c + d*x) + a)**5*sec(c + d*x)
    F = 63*a**5*atanh(sin(c + d*x))/(8*d) + 63*I*a**5*sec(c + d*x)/(8*d) + 9*I*a**2*(I*a*tan(c + d*x) + a)**3*sec(c + d*x)/(20*d) + I*a*(I*a*tan(c + d*x) + a)**4*sec(c + d*x)/(5*d) + 21*I*a*(I*a**2*tan(c + d*x) + a**2)**2*sec(c + d*x)/(20*d) + 21*I*(I*a**5*tan(c + d*x) + a**5)*sec(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_70():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)
    F = -35*a**5*atanh(sin(c + d*x))/(2*d) - 35*I*a**5*sec(c + d*x)/(2*d) - 7*I*a**3*(I*a*tan(c + d*x) + a)**2*sec(c + d*x)/(3*d) - 2*I*a*(I*a*tan(c + d*x) + a)**4*cos(c + d*x)/d - 35*I*(I*a**5*tan(c + d*x) + a**5)*sec(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_71():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**3
    F = 5*a**5*atanh(sin(c + d*x))/d + 5*I*a**5*sec(c + d*x)/d + 10*I*a**3*(I*a*tan(c + d*x) + a)**2*cos(c + d*x)/(3*d) - 2*I*a*(I*a*tan(c + d*x) + a)**4*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_72():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**5
    F = -I*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_73():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**7
    F = -2*I*a**2*(I*a*tan(c + d*x) + a)**3*cos(c + d*x)**3/(105*d) - 2*I*a*(I*a*tan(c + d*x) + a)**4*cos(c + d*x)**5/(35*d) - I*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_74():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**9
    F = a**5*sin(c + d*x)**5/(105*d) - 2*a**5*sin(c + d*x)**3/(63*d) + a**5*sin(c + d*x)/(21*d) - I*a**5*cos(c + d*x)**5/(105*d) - 2*I*a**3*(I*a*tan(c + d*x) + a)**2*cos(c + d*x)**7/(63*d) - 2*I*a*(I*a*tan(c + d*x) + a)**4*cos(c + d*x)**9/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_75():
    f = (I*a*tan(c + d*x) + a)**5*cos(c + d*x)**11
    F = -5*a**5*sin(c + d*x)**7/(231*d) + a**5*sin(c + d*x)**5/(11*d) - 5*a**5*sin(c + d*x)**3/(33*d) + 5*a**5*sin(c + d*x)/(33*d) - 5*I*a**5*cos(c + d*x)**7/(231*d) - 2*I*a**3*(I*a*tan(c + d*x) + a)**2*cos(c + d*x)**9/(33*d) - 2*I*a*(I*a*tan(c + d*x) + a)**4*cos(c + d*x)**11/(11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_76():
    f = (I*a*tan(c + d*x) + a)**8*sec(c + d*x)**8
    F = -2*I*(I*a*tan(c + d*x) + a)**12/(3*a**4*d) + 12*I*(I*a*tan(c + d*x) + a)**13/(13*a**5*d) - 3*I*(I*a*tan(c + d*x) + a)**14/(7*a**6*d) + I*(I*a*tan(c + d*x) + a)**15/(15*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_77():
    f = (I*a*tan(c + d*x) + a)**8*sec(c + d*x)**6
    F = -4*I*(I*a*tan(c + d*x) + a)**11/(11*a**3*d) + I*(I*a*tan(c + d*x) + a)**12/(3*a**4*d) - I*(I*a*tan(c + d*x) + a)**13/(13*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_78():
    f = (I*a*tan(c + d*x) + a)**8*sec(c + d*x)**4
    F = -I*(I*a*tan(c + d*x) + a)**10/(5*a**2*d) + I*(I*a*tan(c + d*x) + a)**11/(11*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_79():
    f = (I*a*tan(c + d*x) + a)**8*sec(c + d*x)**2
    F = -I*(I*a*tan(c + d*x) + a)**9/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_80():
    f = (I*a*tan(c + d*x) + a)**8
    F = 128*a**8*x - 128*I*a**8*log(cos(c + d*x))/d - 64*a**8*tan(c + d*x)/d + 4*I*a**3*(I*a*tan(c + d*x) + a)**5/(5*d) + I*a**2*(I*a*tan(c + d*x) + a)**6/(3*d) + 16*I*a**2*(I*a**2*tan(c + d*x) + a**2)**3/(3*d) + I*a*(I*a*tan(c + d*x) + a)**7/(7*d) + 2*I*(I*a**2*tan(c + d*x) + a**2)**4/d + 16*I*(I*a**4*tan(c + d*x) + a**4)**2/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_81():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**2
    F = -64*I*a**9/(d*(-I*a*tan(c + d*x) + a)) - 192*a**8*x + 192*I*a**8*log(cos(c + d*x))/d + a**8*tan(c + d*x)**5/(5*d) - 2*I*a**8*tan(c + d*x)**4/d - 10*a**8*tan(c + d*x)**3/d + 36*I*a**8*tan(c + d*x)**2/d + 129*a**8*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_82():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**4
    F = -16*I*a**10/(d*(-I*a*tan(c + d*x) + a)**2) + 80*I*a**9/(d*(-I*a*tan(c + d*x) + a)) + 80*a**8*x - 80*I*a**8*log(cos(c + d*x))/d + a**8*tan(c + d*x)**3/(3*d) - 4*I*a**8*tan(c + d*x)**2/d - 31*a**8*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_83():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**6
    F = -16*I*a**11/(3*d*(-I*a*tan(c + d*x) + a)**3) + 16*I*a**10/(d*(-I*a*tan(c + d*x) + a)**2) - 24*I*a**9/(d*(-I*a*tan(c + d*x) + a)) - 8*a**8*x + 8*I*a**8*log(cos(c + d*x))/d + a**8*tan(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_84():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**8
    F = -I*(I*a**3*tan(c + d*x) + a**3)**4/(8*d*(-I*a*tan(c + d*x) + a)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_85():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**10
    F = -4*I*a**13/(5*d*(-I*a*tan(c + d*x) + a)**5) + I*a**12/(d*(-I*a*tan(c + d*x) + a)**4) - I*a**11/(3*d*(-I*a*tan(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_86():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**12
    F = -I*a**14/(3*d*(-I*a*tan(c + d*x) + a)**6) + I*a**13/(5*d*(-I*a*tan(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_87():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**14
    F = -I*a**15/(7*d*(-I*a*tan(c + d*x) + a)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_88():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**16
    F = -I*a**16/(16*d*(-I*a*tan(c + d*x) + a)**8) - I*a**15/(28*d*(-I*a*tan(c + d*x) + a)**7) - I*a**14/(48*d*(-I*a*tan(c + d*x) + a)**6) - I*a**13/(80*d*(-I*a*tan(c + d*x) + a)**5) - I*a**12/(128*d*(-I*a*tan(c + d*x) + a)**4) - I*a**11/(192*d*(-I*a*tan(c + d*x) + a)**3) - I*a**10/(256*d*(-I*a*tan(c + d*x) + a)**2) - I*a**9/(256*d*(-I*a*tan(c + d*x) + a)) + a**8*x/256
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_89():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**18
    F = -I*a**17/(36*d*(-I*a*tan(c + d*x) + a)**9) - I*a**16/(32*d*(-I*a*tan(c + d*x) + a)**8) - 3*I*a**15/(112*d*(-I*a*tan(c + d*x) + a)**7) - I*a**14/(48*d*(-I*a*tan(c + d*x) + a)**6) - I*a**13/(64*d*(-I*a*tan(c + d*x) + a)**5) - 3*I*a**12/(256*d*(-I*a*tan(c + d*x) + a)**4) - 7*I*a**11/(768*d*(-I*a*tan(c + d*x) + a)**3) - I*a**10/(128*d*(-I*a*tan(c + d*x) + a)**2) + I*a**9/(1024*d*(I*a*tan(c + d*x) + a)) - 9*I*a**9/(1024*d*(-I*a*tan(c + d*x) + a)) + 5*a**8*x/512
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_90():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)
    F = -3003*a**8*atanh(sin(c + d*x))/(16*d) - 3003*I*a**8*sec(c + d*x)/(16*d) - 13*I*a**3*(I*a*tan(c + d*x) + a)**5*sec(c + d*x)/(6*d) - 429*I*a**2*(I*a**2*tan(c + d*x) + a**2)**3*sec(c + d*x)/(40*d) - 2*I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)/d - 143*I*(I*a**2*tan(c + d*x) + a**2)**4*sec(c + d*x)/(30*d) - 1001*I*(I*a**4*tan(c + d*x) + a**4)**2*sec(c + d*x)/(40*d) - 1001*I*(I*a**8*tan(c + d*x) + a**8)*sec(c + d*x)/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_91():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**3
    F = 1155*a**8*atanh(sin(c + d*x))/(8*d) + 1155*I*a**8*sec(c + d*x)/(8*d) + 22*I*a**3*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)/(3*d) + 33*I*a**2*(I*a**2*tan(c + d*x) + a**2)**3*sec(c + d*x)/(4*d) - 2*I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)**3/(3*d) + 77*I*(I*a**4*tan(c + d*x) + a**4)**2*sec(c + d*x)/(4*d) + 385*I*(I*a**8*tan(c + d*x) + a**8)*sec(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_92():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**5
    F = -63*a**8*atanh(sin(c + d*x))/(2*d) - 63*I*a**8*sec(c + d*x)/(2*d) + 6*I*a**3*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)**3/(5*d) - 42*I*a**2*(I*a**2*tan(c + d*x) + a**2)**3*cos(c + d*x)/(5*d) - 2*I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)**5/(5*d) - 21*I*(I*a**8*tan(c + d*x) + a**8)*sec(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_93():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**7
    F = a**8*atanh(sin(c + d*x))/d + 2*I*a**3*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)**5/(5*d) - 2*I*a**2*(I*a**2*tan(c + d*x) + a**2)**3*cos(c + d*x)**3/(3*d) - 2*I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)**7/(7*d) + 2*I*(I*a**8*tan(c + d*x) + a**8)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_94():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**9
    F = -I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)**7/(63*d) - I*(I*a*tan(c + d*x) + a)**8*cos(c + d*x)**9/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_95():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**11
    F = -2*I*a**3*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)**5/(1155*d) - 2*I*a**2*(I*a*tan(c + d*x) + a)**6*cos(c + d*x)**7/(231*d) - I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)**9/(33*d) - I*(I*a*tan(c + d*x) + a)**8*cos(c + d*x)**11/(11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_96():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**13
    F = -20*I*a**3*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)**7/(3003*d) - 20*I*a**2*(I*a*tan(c + d*x) + a)**6*cos(c + d*x)**9/(1287*d) - 8*I*a**2*(I*a**2*tan(c + d*x) + a**2)**3*cos(c + d*x)**3/(9009*d) - 5*I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)**11/(143*d) - I*(I*a*tan(c + d*x) + a)**8*cos(c + d*x)**13/(13*d) - 8*I*(I*a**2*tan(c + d*x) + a**2)**4*cos(c + d*x)**5/(3003*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_97():
    f = (I*a*tan(c + d*x) + a)**8*cos(c + d*x)**15
    F = -a**8*sin(c + d*x)**7/(1287*d) + 7*a**8*sin(c + d*x)**5/(2145*d) - 7*a**8*sin(c + d*x)**3/(1287*d) + 7*a**8*sin(c + d*x)/(1287*d) - 2*I*a**3*(I*a*tan(c + d*x) + a)**5*cos(c + d*x)**13/(195*d) - 2*I*a**2*(I*a**2*tan(c + d*x) + a**2)**3*cos(c + d*x)**11/(715*d) - 2*I*a*(I*a*tan(c + d*x) + a)**7*cos(c + d*x)**15/(15*d) - 2*I*(I*a**8*tan(c + d*x) + a**8)*cos(c + d*x)**9/(1287*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_98():
    f = sec(c + d*x)**10/(I*a*tan(c + d*x) + a)
    F = 8*I*(-I*a*tan(c + d*x) + a)**5/(5*a**6*d) - 2*I*(-I*a*tan(c + d*x) + a)**6/(a**7*d) + 6*I*(-I*a*tan(c + d*x) + a)**7/(7*a**8*d) - I*(-I*a*tan(c + d*x) + a)**8/(8*a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_99():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)
    F = I*(-I*a*tan(c + d*x) + a)**4/(a**5*d) - 4*I*(-I*a*tan(c + d*x) + a)**5/(5*a**6*d) + I*(-I*a*tan(c + d*x) + a)**6/(6*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_100():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)
    F = 2*I*(-I*a*tan(c + d*x) + a)**3/(3*a**4*d) - I*(-I*a*tan(c + d*x) + a)**4/(4*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_101():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)
    F = -I*tan(c + d*x)**2/(2*a*d) + tan(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_102():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)
    F = x/a + I*log(cos(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_103():
    f = 1/(I*a*tan(c + d*x) + a)
    F = I/(2*d*(I*a*tan(c + d*x) + a)) + x/(2*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_104():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)
    F = I*a/(8*d*(I*a*tan(c + d*x) + a)**2) + I/(4*d*(I*a*tan(c + d*x) + a)) - I/(8*d*(-I*a*tan(c + d*x) + a)) + 3*x/(8*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_105():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)
    F = I*a**2/(24*d*(I*a*tan(c + d*x) + a)**3) + 3*I*a/(32*d*(I*a*tan(c + d*x) + a)**2) - I*a/(32*d*(-I*a*tan(c + d*x) + a)**2) + 3*I/(16*d*(I*a*tan(c + d*x) + a)) - I/(8*d*(-I*a*tan(c + d*x) + a)) + 5*x/(16*a)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_106():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)
    F = tan(c + d*x)*sec(c + d*x)**3/(4*a*d) + 3*tan(c + d*x)*sec(c + d*x)/(8*a*d) + 3*atanh(sin(c + d*x))/(8*a*d) - I*sec(c + d*x)**5/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_107():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)
    F = tan(c + d*x)*sec(c + d*x)/(2*a*d) + atanh(sin(c + d*x))/(2*a*d) - I*sec(c + d*x)**3/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_108():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)
    F = atanh(sin(c + d*x))/(a*d) - I*sec(c + d*x)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_109():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)
    F = I*sec(c + d*x)/(d*(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_110():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)
    F = I*cos(c + d*x)/(3*d*(I*a*tan(c + d*x) + a)) + 2*sin(c + d*x)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_111():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)
    F = I*cos(c + d*x)**3/(5*d*(I*a*tan(c + d*x) + a)) - 4*sin(c + d*x)**3/(15*a*d) + 4*sin(c + d*x)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_112():
    f = cos(c + d*x)**5/(I*a*tan(c + d*x) + a)
    F = I*cos(c + d*x)**5/(7*d*(I*a*tan(c + d*x) + a)) + 6*sin(c + d*x)**5/(35*a*d) - 4*sin(c + d*x)**3/(7*a*d) + 6*sin(c + d*x)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_113():
    f = sec(c + d*x)**10/(I*a*tan(c + d*x) + a)**2
    F = 4*I*(-I*a*tan(c + d*x) + a)**5/(5*a**7*d) - 2*I*(-I*a*tan(c + d*x) + a)**6/(3*a**8*d) + I*(-I*a*tan(c + d*x) + a)**7/(7*a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_114():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)**2
    F = I*(-I*a*tan(c + d*x) + a)**4/(2*a**6*d) - I*(-I*a*tan(c + d*x) + a)**5/(5*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_115():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)**2
    F = I*(-I*a*tan(c + d*x) + a)**3/(3*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_116():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)**2
    F = 2*x/a**2 + 2*I*log(cos(c + d*x))/(a**2*d) - tan(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_117():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)**2
    F = I/(d*(I*a**2*tan(c + d*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_118():
    f = (I*a*tan(c + d*x) + a)**(-2)
    F = I/(4*d*(I*a**2*tan(c + d*x) + a**2)) + I/(4*d*(I*a*tan(c + d*x) + a)**2) + x/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_119():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)**2
    F = I*a/(12*d*(I*a*tan(c + d*x) + a)**3) + 3*I/(16*d*(I*a**2*tan(c + d*x) + a**2)) - I/(16*d*(-I*a**2*tan(c + d*x) + a**2)) + I/(8*d*(I*a*tan(c + d*x) + a)**2) + x/(4*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_120():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)**2
    F = I*a**2/(32*d*(I*a*tan(c + d*x) + a)**4) + I*a/(16*d*(I*a*tan(c + d*x) + a)**3) + 5*I/(32*d*(I*a**2*tan(c + d*x) + a**2)) - 5*I/(64*d*(-I*a**2*tan(c + d*x) + a**2)) + 3*I/(32*d*(I*a*tan(c + d*x) + a)**2) - I/(64*d*(-I*a*tan(c + d*x) + a)**2) + 15*x/(64*a**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_121():
    f = sec(c + d*x)**9/(I*a*tan(c + d*x) + a)**2
    F = -2*I*sec(c + d*x)**7/(5*d*(I*a**2*tan(c + d*x) + a**2)) + 7*tan(c + d*x)*sec(c + d*x)**5/(30*a**2*d) + 7*tan(c + d*x)*sec(c + d*x)**3/(24*a**2*d) + 7*tan(c + d*x)*sec(c + d*x)/(16*a**2*d) + 7*atanh(sin(c + d*x))/(16*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_122():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)**2
    F = -2*I*sec(c + d*x)**5/(3*d*(I*a**2*tan(c + d*x) + a**2)) + 5*tan(c + d*x)*sec(c + d*x)**3/(12*a**2*d) + 5*tan(c + d*x)*sec(c + d*x)/(8*a**2*d) + 5*atanh(sin(c + d*x))/(8*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_123():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)**2
    F = -2*I*sec(c + d*x)**3/(d*(I*a**2*tan(c + d*x) + a**2)) + 3*tan(c + d*x)*sec(c + d*x)/(2*a**2*d) + 3*atanh(sin(c + d*x))/(2*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_124():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)**2
    F = 2*I*sec(c + d*x)/(d*(I*a**2*tan(c + d*x) + a**2)) - atanh(sin(c + d*x))/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_125():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)**2
    F = I*sec(c + d*x)/(3*d*(I*a**2*tan(c + d*x) + a**2)) + I*sec(c + d*x)/(3*d*(I*a*tan(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_126():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)**2
    F = 2*I*cos(c + d*x)**3/(5*d*(I*a**2*tan(c + d*x) + a**2)) - sin(c + d*x)**3/(5*a**2*d) + 3*sin(c + d*x)/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_127():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)**2
    F = 2*I*cos(c + d*x)**5/(7*d*(I*a**2*tan(c + d*x) + a**2)) + sin(c + d*x)**5/(7*a**2*d) - 10*sin(c + d*x)**3/(21*a**2*d) + 5*sin(c + d*x)/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_128():
    f = cos(c + d*x)**5/(I*a*tan(c + d*x) + a)**2
    F = 2*I*cos(c + d*x)**7/(9*d*(I*a**2*tan(c + d*x) + a**2)) - sin(c + d*x)**7/(9*a**2*d) + 7*sin(c + d*x)**5/(15*a**2*d) - 7*sin(c + d*x)**3/(9*a**2*d) + 7*sin(c + d*x)/(9*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_129():
    f = sec(c + d*x)**14/(I*a*tan(c + d*x) + a)**3
    F = 8*I*(-I*a*tan(c + d*x) + a)**7/(7*a**10*d) - 3*I*(-I*a*tan(c + d*x) + a)**8/(2*a**11*d) + 2*I*(-I*a*tan(c + d*x) + a)**9/(3*a**12*d) - I*(-I*a*tan(c + d*x) + a)**10/(10*a**13*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_130():
    f = sec(c + d*x)**12/(I*a*tan(c + d*x) + a)**3
    F = 2*I*(-I*a*tan(c + d*x) + a)**6/(3*a**9*d) - 4*I*(-I*a*tan(c + d*x) + a)**7/(7*a**10*d) + I*(-I*a*tan(c + d*x) + a)**8/(8*a**11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_131():
    f = sec(c + d*x)**10/(I*a*tan(c + d*x) + a)**3
    F = 2*I*(-I*a*tan(c + d*x) + a)**5/(5*a**8*d) - I*(-I*a*tan(c + d*x) + a)**6/(6*a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_132():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)**3
    F = I*(-I*a*tan(c + d*x) + a)**4/(4*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_133():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)**3
    F = 4*x/a**3 + 4*I*log(cos(c + d*x))/(a**3*d) + I*tan(c + d*x)**2/(2*a**3*d) - 3*tan(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_134():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)**3
    F = 2*I/(d*(I*a**3*tan(c + d*x) + a**3)) - x/a**3 - I*log(cos(c + d*x))/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_135():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)**3
    F = I/(2*a*d*(I*a*tan(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_136():
    f = (I*a*tan(c + d*x) + a)**(-3)
    F = I/(8*d*(I*a**3*tan(c + d*x) + a**3)) + I/(6*d*(I*a*tan(c + d*x) + a)**3) + I/(8*a*d*(I*a*tan(c + d*x) + a)**2) + x/(8*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_137():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)**3
    F = I*a/(16*d*(I*a*tan(c + d*x) + a)**4) + I/(8*d*(I*a**3*tan(c + d*x) + a**3)) - I/(32*d*(-I*a**3*tan(c + d*x) + a**3)) + I/(12*d*(I*a*tan(c + d*x) + a)**3) + 3*I/(32*a*d*(I*a*tan(c + d*x) + a)**2) + 5*x/(32*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_138():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)**3
    F = I*a**2/(40*d*(I*a*tan(c + d*x) + a)**5) + 3*I*a/(64*d*(I*a*tan(c + d*x) + a)**4) + 15*I/(128*d*(I*a**3*tan(c + d*x) + a**3)) - 3*I/(64*d*(-I*a**3*tan(c + d*x) + a**3)) + I/(16*d*(I*a*tan(c + d*x) + a)**3) + 5*I/(64*a*d*(I*a*tan(c + d*x) + a)**2) - I/(128*a*d*(-I*a*tan(c + d*x) + a)**2) + 21*x/(128*a**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_139():
    f = sec(c + d*x)**9/(I*a*tan(c + d*x) + a)**3
    F = -2*I*sec(c + d*x)**7/(3*a*d*(I*a*tan(c + d*x) + a)**2) + 7*tan(c + d*x)*sec(c + d*x)**3/(12*a**3*d) + 7*tan(c + d*x)*sec(c + d*x)/(8*a**3*d) + 7*atanh(sin(c + d*x))/(8*a**3*d) - 7*I*sec(c + d*x)**5/(15*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_140():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)**3
    F = -2*I*sec(c + d*x)**5/(a*d*(I*a*tan(c + d*x) + a)**2) + 5*tan(c + d*x)*sec(c + d*x)/(2*a**3*d) + 5*atanh(sin(c + d*x))/(2*a**3*d) - 5*I*sec(c + d*x)**3/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_141():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)**3
    F = 2*I*sec(c + d*x)**3/(a*d*(I*a*tan(c + d*x) + a)**2) - 3*atanh(sin(c + d*x))/(a**3*d) + 3*I*sec(c + d*x)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_142():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)**3
    F = I*sec(c + d*x)**3/(3*d*(I*a*tan(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_143():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)**3
    F = 2*I*sec(c + d*x)/(15*d*(I*a**3*tan(c + d*x) + a**3)) + I*sec(c + d*x)/(5*d*(I*a*tan(c + d*x) + a)**3) + 2*I*sec(c + d*x)/(15*a*d*(I*a*tan(c + d*x) + a)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_144():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)**3
    F = 8*I*cos(c + d*x)**3/(35*d*(I*a**3*tan(c + d*x) + a**3)) + I*cos(c + d*x)/(7*d*(I*a*tan(c + d*x) + a)**3) - 4*sin(c + d*x)**3/(35*a**3*d) + 12*sin(c + d*x)/(35*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_145():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)**3
    F = 4*I*cos(c + d*x)**5/(21*d*(I*a**3*tan(c + d*x) + a**3)) + I*cos(c + d*x)**3/(9*d*(I*a*tan(c + d*x) + a)**3) + 2*sin(c + d*x)**5/(21*a**3*d) - 20*sin(c + d*x)**3/(63*a**3*d) + 10*sin(c + d*x)/(21*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_146():
    f = cos(c + d*x)**5/(I*a*tan(c + d*x) + a)**3
    F = 16*I*cos(c + d*x)**7/(99*d*(I*a**3*tan(c + d*x) + a**3)) + I*cos(c + d*x)**5/(11*d*(I*a*tan(c + d*x) + a)**3) - 8*sin(c + d*x)**7/(99*a**3*d) + 56*sin(c + d*x)**5/(165*a**3*d) - 56*sin(c + d*x)**3/(99*a**3*d) + 56*sin(c + d*x)/(99*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_147():
    f = sec(c + d*x)**14/(I*a*tan(c + d*x) + a)**4
    F = 4*I*(-I*a*tan(c + d*x) + a)**7/(7*a**11*d) - I*(-I*a*tan(c + d*x) + a)**8/(2*a**12*d) + I*(-I*a*tan(c + d*x) + a)**9/(9*a**13*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_148():
    f = sec(c + d*x)**12/(I*a*tan(c + d*x) + a)**4
    F = I*(-I*a*tan(c + d*x) + a)**6/(3*a**10*d) - I*(-I*a*tan(c + d*x) + a)**7/(7*a**11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_149():
    f = sec(c + d*x)**10/(I*a*tan(c + d*x) + a)**4
    F = I*(-I*a*tan(c + d*x) + a)**5/(5*a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_150():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)**4
    F = 8*x/a**4 + 8*I*log(cos(c + d*x))/(a**4*d) - 4*tan(c + d*x)/(a**4*d) - I*(-I*a*tan(c + d*x) + a)**2/(a**6*d) - I*(-I*a*tan(c + d*x) + a)**3/(3*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_151():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)**4
    F = 4*I/(d*(I*a**4*tan(c + d*x) + a**4)) - 4*x/a**4 - 4*I*log(cos(c + d*x))/(a**4*d) + tan(c + d*x)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_152():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)**4
    F = tan(c + d*x)/(d*(I*a**2*tan(c + d*x) + a**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_153():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)**4
    F = I/(3*a*d*(I*a*tan(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_154():
    f = (I*a*tan(c + d*x) + a)**(-4)
    F = I/(16*d*(I*a**4*tan(c + d*x) + a**4)) + I/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) + I/(8*d*(I*a*tan(c + d*x) + a)**4) + I/(12*a*d*(I*a*tan(c + d*x) + a)**3) + x/(16*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_155():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)**4
    F = I*a/(20*d*(I*a*tan(c + d*x) + a)**5) + 5*I/(64*d*(I*a**4*tan(c + d*x) + a**4)) - I/(64*d*(-I*a**4*tan(c + d*x) + a**4)) + I/(16*d*(I*a**2*tan(c + d*x) + a**2)**2) + I/(16*d*(I*a*tan(c + d*x) + a)**4) + I/(16*a*d*(I*a*tan(c + d*x) + a)**3) + 3*x/(32*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_156():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)**4
    F = I*a**2/(48*d*(I*a*tan(c + d*x) + a)**6) + 3*I*a/(80*d*(I*a*tan(c + d*x) + a)**5) + 21*I/(256*d*(I*a**4*tan(c + d*x) + a**4)) - 7*I/(256*d*(-I*a**4*tan(c + d*x) + a**4)) + 15*I/(256*d*(I*a**2*tan(c + d*x) + a**2)**2) - I/(256*d*(-I*a**2*tan(c + d*x) + a**2)**2) + 3*I/(64*d*(I*a*tan(c + d*x) + a)**4) + 5*I/(96*a*d*(I*a*tan(c + d*x) + a)**3) + 7*x/(64*a**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_157():
    f = sec(c + d*x)**9/(I*a*tan(c + d*x) + a)**4
    F = -14*I*sec(c + d*x)**5/(3*d*(I*a**4*tan(c + d*x) + a**4)) - 2*I*sec(c + d*x)**7/(a*d*(I*a*tan(c + d*x) + a)**3) + 35*tan(c + d*x)*sec(c + d*x)**3/(12*a**4*d) + 35*tan(c + d*x)*sec(c + d*x)/(8*a**4*d) + 35*atanh(sin(c + d*x))/(8*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_158():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)**4
    F = 10*I*sec(c + d*x)**3/(d*(I*a**4*tan(c + d*x) + a**4)) + 2*I*sec(c + d*x)**5/(a*d*(I*a*tan(c + d*x) + a)**3) - 15*tan(c + d*x)*sec(c + d*x)/(2*a**4*d) - 15*atanh(sin(c + d*x))/(2*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_159():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)**4
    F = -2*I*sec(c + d*x)/(d*(I*a**4*tan(c + d*x) + a**4)) + 2*I*sec(c + d*x)**3/(3*a*d*(I*a*tan(c + d*x) + a)**3) + atanh(sin(c + d*x))/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_160():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)**4
    F = I*sec(c + d*x)**3/(5*d*(I*a*tan(c + d*x) + a)**4) + I*sec(c + d*x)**3/(15*a*d*(I*a*tan(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_161():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)**4
    F = 2*I*sec(c + d*x)/(35*d*(I*a**4*tan(c + d*x) + a**4)) + 2*I*sec(c + d*x)/(35*d*(I*a**2*tan(c + d*x) + a**2)**2) + I*sec(c + d*x)/(7*d*(I*a*tan(c + d*x) + a)**4) + 3*I*sec(c + d*x)/(35*a*d*(I*a*tan(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_162():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)**4
    F = 8*I*cos(c + d*x)**3/(63*d*(I*a**4*tan(c + d*x) + a**4)) + I*cos(c + d*x)/(9*d*(I*a*tan(c + d*x) + a)**4) + 5*I*cos(c + d*x)/(63*a*d*(I*a*tan(c + d*x) + a)**3) - 4*sin(c + d*x)**3/(63*a**4*d) + 4*sin(c + d*x)/(21*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_163():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)**4
    F = 4*I*cos(c + d*x)**5/(33*d*(I*a**4*tan(c + d*x) + a**4)) + I*cos(c + d*x)**3/(11*d*(I*a*tan(c + d*x) + a)**4) + 7*I*cos(c + d*x)**3/(99*a*d*(I*a*tan(c + d*x) + a)**3) + 2*sin(c + d*x)**5/(33*a**4*d) - 20*sin(c + d*x)**3/(99*a**4*d) + 10*sin(c + d*x)/(33*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_164():
    f = cos(c + d*x)**5/(I*a*tan(c + d*x) + a)**4
    F = 16*I*cos(c + d*x)**7/(143*d*(I*a**4*tan(c + d*x) + a**4)) + I*cos(c + d*x)**5/(13*d*(I*a*tan(c + d*x) + a)**4) + 9*I*cos(c + d*x)**5/(143*a*d*(I*a*tan(c + d*x) + a)**3) - 8*sin(c + d*x)**7/(143*a**4*d) + 168*sin(c + d*x)**5/(715*a**4*d) - 56*sin(c + d*x)**3/(143*a**4*d) + 56*sin(c + d*x)/(143*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_165():
    f = sec(c + d*x)**14/(I*a*tan(c + d*x) + a)**8
    F = 64*I/(d*(I*a**8*tan(c + d*x) + a**8)) - 192*x/a**8 - 192*I*log(cos(c + d*x))/(a**8*d) + tan(c + d*x)**5/(5*a**8*d) + 2*I*tan(c + d*x)**4/(a**8*d) - 10*tan(c + d*x)**3/(a**8*d) - 36*I*tan(c + d*x)**2/(a**8*d) + 129*tan(c + d*x)/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_166():
    f = sec(c + d*x)**12/(I*a*tan(c + d*x) + a)**8
    F = -80*I/(d*(I*a**8*tan(c + d*x) + a**8)) + 16*I/(d*(I*a**4*tan(c + d*x) + a**4)**2) + 80*x/a**8 + 80*I*log(cos(c + d*x))/(a**8*d) + tan(c + d*x)**3/(3*a**8*d) + 4*I*tan(c + d*x)**2/(a**8*d) - 31*tan(c + d*x)/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_167():
    f = sec(c + d*x)**10/(I*a*tan(c + d*x) + a)**8
    F = 24*I/(d*(I*a**8*tan(c + d*x) + a**8)) - 16*I/(d*(I*a**4*tan(c + d*x) + a**4)**2) + 16*I/(3*a**5*d*(I*a*tan(c + d*x) + a)**3) - 8*x/a**8 - 8*I*log(cos(c + d*x))/(a**8*d) + tan(c + d*x)/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_168():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)**8
    F = I*(-I*a*tan(c + d*x) + a)**4/(8*d*(I*a**3*tan(c + d*x) + a**3)**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_169():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)**8
    F = -I/(d*(I*a**2*tan(c + d*x) + a**2)**4) + 4*I/(5*a**3*d*(I*a*tan(c + d*x) + a)**5) + I/(3*a**5*d*(I*a*tan(c + d*x) + a)**3)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_170():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)**8
    F = I/(3*a**2*d*(I*a*tan(c + d*x) + a)**6) - I/(5*a**3*d*(I*a*tan(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_171():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)**8
    F = I/(7*a*d*(I*a*tan(c + d*x) + a)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_172():
    f = (I*a*tan(c + d*x) + a)**(-8)
    F = I/(256*d*(I*a**8*tan(c + d*x) + a**8)) + I/(256*d*(I*a**4*tan(c + d*x) + a**4)**2) + I/(128*d*(I*a**2*tan(c + d*x) + a**2)**4) + I/(16*d*(I*a*tan(c + d*x) + a)**8) + I/(28*a*d*(I*a*tan(c + d*x) + a)**7) + I/(192*a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) + I/(48*a**2*d*(I*a*tan(c + d*x) + a)**6) + I/(80*a**3*d*(I*a*tan(c + d*x) + a)**5) + x/(256*a**8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_173():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)**8
    F = I*a/(36*d*(I*a*tan(c + d*x) + a)**9) + 9*I/(1024*d*(I*a**8*tan(c + d*x) + a**8)) - I/(1024*d*(-I*a**8*tan(c + d*x) + a**8)) + I/(128*d*(I*a**4*tan(c + d*x) + a**4)**2) + 3*I/(256*d*(I*a**2*tan(c + d*x) + a**2)**4) + I/(32*d*(I*a*tan(c + d*x) + a)**8) + 3*I/(112*a*d*(I*a*tan(c + d*x) + a)**7) + I/(48*a**2*d*(I*a*tan(c + d*x) + a)**6) + I/(64*a**3*d*(I*a*tan(c + d*x) + a)**5) + 7*I/(768*a**5*d*(I*a*tan(c + d*x) + a)**3) + 5*x/(512*a**8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_174():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)**8
    F = I*a**2/(80*d*(I*a*tan(c + d*x) + a)**10) + I*a/(48*d*(I*a*tan(c + d*x) + a)**9) + 55*I/(4096*d*(I*a**8*tan(c + d*x) + a**8)) - 11*I/(4096*d*(-I*a**8*tan(c + d*x) + a**8)) + 45*I/(4096*d*(I*a**4*tan(c + d*x) + a**4)**2) - I/(4096*d*(-I*a**4*tan(c + d*x) + a**4)**2) + 7*I/(512*d*(I*a**2*tan(c + d*x) + a**2)**4) + 3*I/(128*d*(I*a*tan(c + d*x) + a)**8) + 5*I/(224*a*d*(I*a*tan(c + d*x) + a)**7) + 5*I/(256*a**2*d*(I*a*tan(c + d*x) + a)**6) + 21*I/(1280*a**3*d*(I*a*tan(c + d*x) + a)**5) + 3*I/(256*a**5*d*(I*a*tan(c + d*x) + a)**3) + 33*x/(2048*a**8)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_175():
    f = sec(c + d*x)**13/(I*a*tan(c + d*x) + a)**8
    F = -154*I*sec(c + d*x)**5/(d*(I*a**8*tan(c + d*x) + a**8)) + 2*I*sec(c + d*x)**11/(3*a*d*(I*a*tan(c + d*x) + a)**7) - 66*I*sec(c + d*x)**7/(a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) - 22*I*sec(c + d*x)**9/(3*a**3*d*(I*a*tan(c + d*x) + a)**5) + 385*tan(c + d*x)*sec(c + d*x)**3/(4*a**8*d) + 1155*tan(c + d*x)*sec(c + d*x)/(8*a**8*d) + 1155*atanh(sin(c + d*x))/(8*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_176():
    f = sec(c + d*x)**11/(I*a*tan(c + d*x) + a)**8
    F = 42*I*sec(c + d*x)**3/(d*(I*a**8*tan(c + d*x) + a**8)) + 2*I*sec(c + d*x)**9/(5*a*d*(I*a*tan(c + d*x) + a)**7) + 42*I*sec(c + d*x)**5/(5*a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) - 6*I*sec(c + d*x)**7/(5*a**3*d*(I*a*tan(c + d*x) + a)**5) - 63*tan(c + d*x)*sec(c + d*x)/(2*a**8*d) - 63*atanh(sin(c + d*x))/(2*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_177():
    f = sec(c + d*x)**9/(I*a*tan(c + d*x) + a)**8
    F = -2*I*sec(c + d*x)/(d*(I*a**8*tan(c + d*x) + a**8)) + 2*I*sec(c + d*x)**7/(7*a*d*(I*a*tan(c + d*x) + a)**7) + 2*I*sec(c + d*x)**3/(3*a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) - 2*I*sec(c + d*x)**5/(5*a**3*d*(I*a*tan(c + d*x) + a)**5) + atanh(sin(c + d*x))/(a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_178():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)**8
    F = I*sec(c + d*x)**7/(9*d*(I*a*tan(c + d*x) + a)**8) + I*sec(c + d*x)**7/(63*a*d*(I*a*tan(c + d*x) + a)**7)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_179():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)**8
    F = I*sec(c + d*x)**5/(11*d*(I*a*tan(c + d*x) + a)**8) + I*sec(c + d*x)**5/(33*a*d*(I*a*tan(c + d*x) + a)**7) + 2*I*sec(c + d*x)**5/(231*a**2*d*(I*a*tan(c + d*x) + a)**6) + 2*I*sec(c + d*x)**5/(1155*a**3*d*(I*a*tan(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_180():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)**8
    F = 8*I*sec(c + d*x)**3/(3003*d*(I*a**2*tan(c + d*x) + a**2)**4) + I*sec(c + d*x)**3/(13*d*(I*a*tan(c + d*x) + a)**8) + 5*I*sec(c + d*x)**3/(143*a*d*(I*a*tan(c + d*x) + a)**7) + 8*I*sec(c + d*x)**3/(9009*a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) + 20*I*sec(c + d*x)**3/(1287*a**2*d*(I*a*tan(c + d*x) + a)**6) + 20*I*sec(c + d*x)**3/(3003*a**3*d*(I*a*tan(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_181():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)**8
    F = 16*I*sec(c + d*x)/(6435*d*(I*a**8*tan(c + d*x) + a**8)) + 16*I*sec(c + d*x)/(6435*d*(I*a**4*tan(c + d*x) + a**4)**2) + 8*I*sec(c + d*x)/(1287*d*(I*a**2*tan(c + d*x) + a**2)**4) + I*sec(c + d*x)/(15*d*(I*a*tan(c + d*x) + a)**8) + 7*I*sec(c + d*x)/(195*a*d*(I*a*tan(c + d*x) + a)**7) + 8*I*sec(c + d*x)/(2145*a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) + 14*I*sec(c + d*x)/(715*a**2*d*(I*a*tan(c + d*x) + a)**6) + 14*I*sec(c + d*x)/(1287*a**3*d*(I*a*tan(c + d*x) + a)**5)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_182():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)**8
    F = 128*I*cos(c + d*x)**3/(12155*d*(I*a**8*tan(c + d*x) + a**8)) + 112*I*cos(c + d*x)/(12155*d*(I*a**2*tan(c + d*x) + a**2)**4) + I*cos(c + d*x)/(17*d*(I*a*tan(c + d*x) + a)**8) + 3*I*cos(c + d*x)/(85*a*d*(I*a*tan(c + d*x) + a)**7) + 16*I*cos(c + d*x)/(2431*a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) + 24*I*cos(c + d*x)/(1105*a**2*d*(I*a*tan(c + d*x) + a)**6) + 168*I*cos(c + d*x)/(12155*a**3*d*(I*a*tan(c + d*x) + a)**5) - 64*sin(c + d*x)**3/(12155*a**8*d) + 192*sin(c + d*x)/(12155*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_183():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)**8
    F = 64*I*cos(c + d*x)**5/(4199*d*(I*a**8*tan(c + d*x) + a**8)) + 48*I*cos(c + d*x)**3/(4199*d*(I*a**2*tan(c + d*x) + a**2)**4) + I*cos(c + d*x)**3/(19*d*(I*a*tan(c + d*x) + a)**8) + 11*I*cos(c + d*x)**3/(323*a*d*(I*a*tan(c + d*x) + a)**7) + 112*I*cos(c + d*x)**3/(12597*a**2*d*(I*a**2*tan(c + d*x) + a**2)**3) + 22*I*cos(c + d*x)**3/(969*a**2*d*(I*a*tan(c + d*x) + a)**6) + 66*I*cos(c + d*x)**3/(4199*a**3*d*(I*a*tan(c + d*x) + a)**5) + 32*sin(c + d*x)**5/(4199*a**8*d) - 320*sin(c + d*x)**3/(12597*a**8*d) + 160*sin(c + d*x)/(4199*a**8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_184():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a)
    F = -6*a*e**4*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*a*e**3*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*d) + 2*a*e*(e*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*d) + 2*I*a*(e*sec(c + d*x))**(sympy.S(7)/2)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_185():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)
    F = 2*a*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 2*a*e*(e*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*d) + 2*I*a*(e*sec(c + d*x))**(sympy.S(5)/2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_186():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)
    F = -2*a*e**2*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*a*e*sqrt(e*sec(c + d*x))*sin(c + d*x)/d + 2*I*a*(e*sec(c + d*x))**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_187():
    f = sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)
    F = 2*a*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/d + 2*I*a*sqrt(e*sec(c + d*x))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_188():
    f = (I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x))
    F = -2*I*a/(d*sqrt(e*sec(c + d*x))) + 2*a*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_189():
    f = (I*a*tan(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*I*a/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2)) + 2*a*sin(c + d*x)/(3*d*e*sqrt(e*sec(c + d*x))) + 2*a*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_190():
    f = (I*a*tan(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*I*a/(5*d*(e*sec(c + d*x))**(sympy.S(5)/2)) + 2*a*sin(c + d*x)/(5*d*e*(e*sec(c + d*x))**(sympy.S(3)/2)) + 6*a*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_191():
    f = (I*a*tan(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(7)/2)
    F = -2*I*a/(7*d*(e*sec(c + d*x))**(sympy.S(7)/2)) + 2*a*sin(c + d*x)/(7*d*e*(e*sec(c + d*x))**(sympy.S(5)/2)) + 10*a*sin(c + d*x)/(21*d*e**3*sqrt(e*sec(c + d*x))) + 10*a*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*e**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_192():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**2
    F = -14*a**2*e**2*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 14*a**2*e*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*d) + 14*I*a**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(15*d) + 2*I*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**2*tan(c + d*x) + a**2)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_193():
    f = sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**2
    F = 10*a**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d) + 10*I*a**2*sqrt(e*sec(c + d*x))/(3*d) + 2*I*sqrt(e*sec(c + d*x))*(I*a**2*tan(c + d*x) + a**2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_194():
    f = (I*a*tan(c + d*x) + a)**2/sqrt(e*sec(c + d*x))
    F = 6*a**2*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 6*a**2*sqrt(e*sec(c + d*x))*sin(c + d*x)/(d*e) - 4*I*(I*a**2*tan(c + d*x) + a**2)/(d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_195():
    f = (I*a*tan(c + d*x) + a)**2/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*a**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2) - 4*I*(I*a**2*tan(c + d*x) + a**2)/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_196():
    f = (I*a*tan(c + d*x) + a)**2/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = 2*a**2*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 4*I*(I*a**2*tan(c + d*x) + a**2)/(5*d*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_197():
    f = (I*a*tan(c + d*x) + a)**2/(e*sec(c + d*x))**(sympy.S(7)/2)
    F = 2*a**2*sin(c + d*x)/(7*d*e**3*sqrt(e*sec(c + d*x))) + 2*a**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*d*e**4) - 4*I*(I*a**2*tan(c + d*x) + a**2)/(7*d*(e*sec(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_198():
    f = (I*a*tan(c + d*x) + a)**2/(e*sec(c + d*x))**(sympy.S(9)/2)
    F = 2*a**2*sin(c + d*x)/(9*d*e**3*(e*sec(c + d*x))**(sympy.S(3)/2)) + 2*a**2*elliptic_e(c/2 + d*x/2, 2)/(3*d*e**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 4*I*(I*a**2*tan(c + d*x) + a**2)/(9*d*(e*sec(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_199():
    f = (I*a*tan(c + d*x) + a)**2/(e*sec(c + d*x))**(sympy.S(11)/2)
    F = 2*a**2*sin(c + d*x)/(11*d*e**3*(e*sec(c + d*x))**(sympy.S(5)/2)) + 10*a**2*sin(c + d*x)/(33*d*e**5*sqrt(e*sec(c + d*x))) + 10*a**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(33*d*e**6) - 4*I*(I*a**2*tan(c + d*x) + a**2)/(11*d*(e*sec(c + d*x))**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_200():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a)**3
    F = -2*a**3*e**4*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*a**3*e**3*sqrt(e*sec(c + d*x))*sin(c + d*x)/d + 2*a**3*e*(e*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(3*d) + 10*I*a**3*(e*sec(c + d*x))**(sympy.S(7)/2)/(21*d) + 2*I*a*(e*sec(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a)**2/(11*d) + 10*I*(e*sec(c + d*x))**(sympy.S(7)/2)*(I*a**3*tan(c + d*x) + a**3)/(33*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_201():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**3
    F = 26*a**3*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d) + 26*a**3*e*(e*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(21*d) + 26*I*a**3*(e*sec(c + d*x))**(sympy.S(5)/2)/(35*d) + 2*I*a*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**2/(9*d) + 26*I*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*a**3*tan(c + d*x) + a**3)/(63*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_202():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**3
    F = -22*a**3*e**2*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 22*a**3*e*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*d) + 22*I*a**3*(e*sec(c + d*x))**(sympy.S(3)/2)/(15*d) + 2*I*a*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**2/(7*d) + 22*I*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**3*tan(c + d*x) + a**3)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_203():
    f = sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**3
    F = 6*a**3*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/d + 6*I*a**3*sqrt(e*sec(c + d*x))/d + 2*I*a*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**2/(5*d) + 6*I*sqrt(e*sec(c + d*x))*(I*a**3*tan(c + d*x) + a**3)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_204():
    f = (I*a*tan(c + d*x) + a)**3/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = -10*a**3*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*d*e**2) - 10*I*a**3*sqrt(e*sec(c + d*x))/(3*d*e**2) - 4*I*a*(I*a*tan(c + d*x) + a)**2/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_205():
    f = (I*a*tan(c + d*x) + a)**3/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = 6*I*a**3/(5*d*e**2*sqrt(e*sec(c + d*x))) - 6*a**3*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 4*I*a*(I*a*tan(c + d*x) + a)**2/(5*d*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_206():
    f = (I*a*tan(c + d*x) + a)**3/(e*sec(c + d*x))**(sympy.S(7)/2)
    F = -2*a**3*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*e**4) - 2*I*(I*a*tan(c + d*x) + a)**3/(7*d*(e*sec(c + d*x))**(sympy.S(7)/2)) - 4*I*(I*a**3*tan(c + d*x) + a**3)/(21*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_207():
    f = (I*a*tan(c + d*x) + a)**3/(e*sec(c + d*x))**(sympy.S(9)/2)
    F = 2*a**3*elliptic_e(c/2 + d*x/2, 2)/(15*d*e**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 2*I*(I*a*tan(c + d*x) + a)**3/(9*d*(e*sec(c + d*x))**(sympy.S(9)/2)) - 4*I*(I*a**3*tan(c + d*x) + a**3)/(15*d*e**2*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_208():
    f = (I*a*tan(c + d*x) + a)**3/(e*sec(c + d*x))**(sympy.S(11)/2)
    F = 10*a**3*sin(c + d*x)/(77*d*e**5*sqrt(e*sec(c + d*x))) + 10*a**3*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*d*e**6) - 2*I*(I*a*tan(c + d*x) + a)**3/(11*d*(e*sec(c + d*x))**(sympy.S(11)/2)) - 20*I*(I*a**3*tan(c + d*x) + a**3)/(77*d*e**2*(e*sec(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_209():
    f = (I*a*tan(c + d*x) + a)**3/(e*sec(c + d*x))**(sympy.S(13)/2)
    F = 14*a**3*sin(c + d*x)/(117*d*e**5*(e*sec(c + d*x))**(sympy.S(3)/2)) + 14*a**3*elliptic_e(c/2 + d*x/2, 2)/(39*d*e**6*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 2*I*(I*a*tan(c + d*x) + a)**3/(13*d*(e*sec(c + d*x))**(sympy.S(13)/2)) - 28*I*(I*a**3*tan(c + d*x) + a**3)/(117*d*e**2*(e*sec(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_210():
    f = (I*a*tan(c + d*x) + a)**3/(e*sec(c + d*x))**(sympy.S(15)/2)
    F = 6*a**3*sin(c + d*x)/(55*d*e**5*(e*sec(c + d*x))**(sympy.S(5)/2)) + 2*a**3*sin(c + d*x)/(11*d*e**7*sqrt(e*sec(c + d*x))) + 2*a**3*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(11*d*e**8) - 2*I*(I*a*tan(c + d*x) + a)**3/(15*d*(e*sec(c + d*x))**(sympy.S(15)/2)) - 12*I*(I*a**3*tan(c + d*x) + a**3)/(55*d*e**2*(e*sec(c + d*x))**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_211():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**4
    F = -22*a**4*e**2*elliptic_e(c/2 + d*x/2, 2)/(3*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 22*a**4*e*sqrt(e*sec(c + d*x))*sin(c + d*x)/(3*d) + 22*I*a**4*(e*sec(c + d*x))**(sympy.S(3)/2)/(9*d) + 2*I*a*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**3/(9*d) + 10*I*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**2*tan(c + d*x) + a**2)**2/(21*d) + 22*I*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**4*tan(c + d*x) + a**4)/(21*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_212():
    f = sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**4
    F = 78*a**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*d) + 78*I*a**4*sqrt(e*sec(c + d*x))/(7*d) + 2*I*a*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**3/(7*d) + 26*I*sqrt(e*sec(c + d*x))*(I*a**2*tan(c + d*x) + a**2)**2/(35*d) + 78*I*sqrt(e*sec(c + d*x))*(I*a**4*tan(c + d*x) + a**4)/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_213():
    f = (I*a*tan(c + d*x) + a)**4/sqrt(e*sec(c + d*x))
    F = 154*a**4*elliptic_e(c/2 + d*x/2, 2)/(5*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 154*a**4*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*d*e) - 154*I*a**4*(e*sec(c + d*x))**(sympy.S(3)/2)/(15*d*e**2) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(d*sqrt(e*sec(c + d*x))) - 22*I*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**4*tan(c + d*x) + a**4)/(5*d*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_214():
    f = (I*a*tan(c + d*x) + a)**4/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = -10*a**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*e**2) - 10*I*a**4*sqrt(e*sec(c + d*x))/(d*e**2) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2)) - 2*I*sqrt(e*sec(c + d*x))*(I*a**4*tan(c + d*x) + a**4)/(d*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_215():
    f = (I*a*tan(c + d*x) + a)**4/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = -42*a**4*elliptic_e(c/2 + d*x/2, 2)/(5*d*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 42*a**4*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*d*e**3) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(5*d*(e*sec(c + d*x))**(sympy.S(5)/2)) + 28*I*(I*a**4*tan(c + d*x) + a**4)/(5*d*e**2*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_216():
    f = (I*a*tan(c + d*x) + a)**4/(e*sec(c + d*x))**(sympy.S(7)/2)
    F = 10*a**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*d*e**4) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(7*d*(e*sec(c + d*x))**(sympy.S(7)/2)) + 20*I*(I*a**4*tan(c + d*x) + a**4)/(21*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_217():
    f = (I*a*tan(c + d*x) + a)**4/(e*sec(c + d*x))**(sympy.S(9)/2)
    F = -2*a**4*elliptic_e(c/2 + d*x/2, 2)/(15*d*e**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(9*d*(e*sec(c + d*x))**(sympy.S(9)/2)) + 4*I*(I*a**4*tan(c + d*x) + a**4)/(15*d*e**2*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_218():
    f = (I*a*tan(c + d*x) + a)**4/(e*sec(c + d*x))**(sympy.S(11)/2)
    F = -2*a**4*sin(c + d*x)/(77*d*e**5*sqrt(e*sec(c + d*x))) - 2*a**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*d*e**6) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(11*d*(e*sec(c + d*x))**(sympy.S(11)/2)) + 4*I*(I*a**4*tan(c + d*x) + a**4)/(77*d*e**2*(e*sec(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_219():
    f = (I*a*tan(c + d*x) + a)**4/(e*sec(c + d*x))**(sympy.S(13)/2)
    F = 2*a**4*sin(c + d*x)/(117*d*e**5*(e*sec(c + d*x))**(sympy.S(3)/2)) + 2*a**4*elliptic_e(c/2 + d*x/2, 2)/(39*d*e**6*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(13*d*(e*sec(c + d*x))**(sympy.S(13)/2)) - 4*I*(I*a**4*tan(c + d*x) + a**4)/(117*d*e**2*(e*sec(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_220():
    f = (I*a*tan(c + d*x) + a)**4/(e*sec(c + d*x))**(sympy.S(15)/2)
    F = 2*a**4*sin(c + d*x)/(55*d*e**5*(e*sec(c + d*x))**(sympy.S(5)/2)) + 2*a**4*sin(c + d*x)/(33*d*e**7*sqrt(e*sec(c + d*x))) + 2*a**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(33*d*e**8) - 4*I*a*(I*a*tan(c + d*x) + a)**3/(15*d*(e*sec(c + d*x))**(sympy.S(15)/2)) - 4*I*(I*a**4*tan(c + d*x) + a**4)/(55*d*e**2*(e*sec(c + d*x))**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_221():
    f = (e*sec(c + d*x))**(sympy.S(11)/2)/(I*a*tan(c + d*x) + a)
    F = -6*e**6*elliptic_e(c/2 + d*x/2, 2)/(5*a*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 6*e**5*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*a*d) + 2*e**3*(e*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(5*a*d) - 2*I*e**2*(e*sec(c + d*x))**(sympy.S(7)/2)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_222():
    f = (e*sec(c + d*x))**(sympy.S(9)/2)/(I*a*tan(c + d*x) + a)
    F = 2*e**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a*d) + 2*e**3*(e*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*a*d) - 2*I*e**2*(e*sec(c + d*x))**(sympy.S(5)/2)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_223():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)
    F = -2*e**4*elliptic_e(c/2 + d*x/2, 2)/(a*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 2*e**3*sqrt(e*sec(c + d*x))*sin(c + d*x)/(a*d) - 2*I*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_224():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)
    F = 2*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(a*d) - 2*I*e**2*sqrt(e*sec(c + d*x))/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_225():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)
    F = 2*I*e**2/(a*d*sqrt(e*sec(c + d*x))) + 2*e**2*elliptic_e(c/2 + d*x/2, 2)/(a*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_226():
    f = sqrt(e*sec(c + d*x))/(I*a*tan(c + d*x) + a)
    F = 2*I*sqrt(e*sec(c + d*x))/(3*d*(I*a*tan(c + d*x) + a)) + 2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_227():
    f = 1/(sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a))
    F = 2*I/(5*d*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)) + 6*elliptic_e(c/2 + d*x/2, 2)/(5*a*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_228():
    f = 1/((e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a))
    F = 2*I/(7*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)) + 10*sin(c + d*x)/(21*a*d*e*sqrt(e*sec(c + d*x))) + 10*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a*d*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_229():
    f = 1/((e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a))
    F = 2*I/(9*d*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)) + 14*sin(c + d*x)/(45*a*d*e*(e*sec(c + d*x))**(sympy.S(3)/2)) + 14*elliptic_e(c/2 + d*x/2, 2)/(15*a*d*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_230():
    f = 1/((e*sec(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a))
    F = 2*I/(11*d*(e*sec(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a)) + 18*sin(c + d*x)/(77*a*d*e*(e*sec(c + d*x))**(sympy.S(5)/2)) + 30*sin(c + d*x)/(77*a*d*e**3*sqrt(e*sec(c + d*x))) + 30*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*a*d*e**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_231():
    f = (e*sec(c + d*x))**(sympy.S(15)/2)/(I*a*tan(c + d*x) + a)**2
    F = -4*I*e**2*(e*sec(c + d*x))**(sympy.S(11)/2)/(7*d*(I*a**2*tan(c + d*x) + a**2)) - 22*e**8*elliptic_e(c/2 + d*x/2, 2)/(15*a**2*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 22*e**7*sqrt(e*sec(c + d*x))*sin(c + d*x)/(15*a**2*d) + 22*e**5*(e*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(45*a**2*d) + 22*e**3*(e*sec(c + d*x))**(sympy.S(9)/2)*sin(c + d*x)/(63*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_232():
    f = (e*sec(c + d*x))**(sympy.S(13)/2)/(I*a*tan(c + d*x) + a)**2
    F = -4*I*e**2*(e*sec(c + d*x))**(sympy.S(9)/2)/(5*d*(I*a**2*tan(c + d*x) + a**2)) + 6*e**6*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*a**2*d) + 6*e**5*(e*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(7*a**2*d) + 18*e**3*(e*sec(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)/(35*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_233():
    f = (e*sec(c + d*x))**(sympy.S(11)/2)/(I*a*tan(c + d*x) + a)**2
    F = -4*I*e**2*(e*sec(c + d*x))**(sympy.S(7)/2)/(3*d*(I*a**2*tan(c + d*x) + a**2)) - 14*e**6*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 14*e**5*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*a**2*d) + 14*e**3*(e*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(15*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_234():
    f = (e*sec(c + d*x))**(sympy.S(9)/2)/(I*a*tan(c + d*x) + a)**2
    F = -4*I*e**2*(e*sec(c + d*x))**(sympy.S(5)/2)/(d*(I*a**2*tan(c + d*x) + a**2)) + 10*e**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d) + 10*e**3*(e*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_235():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**2
    F = 4*I*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(d*(I*a**2*tan(c + d*x) + a**2)) + 6*e**4*elliptic_e(c/2 + d*x/2, 2)/(a**2*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 6*e**3*sqrt(e*sec(c + d*x))*sin(c + d*x)/(a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_236():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**2
    F = 4*I*e**2*sqrt(e*sec(c + d*x))/(3*d*(I*a**2*tan(c + d*x) + a**2)) - 2*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_237():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**2
    F = 4*I*e**2/(5*d*sqrt(e*sec(c + d*x))*(I*a**2*tan(c + d*x) + a**2)) + 2*e**2*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_238():
    f = sqrt(e*sec(c + d*x))/(I*a*tan(c + d*x) + a)**2
    F = 4*I*e**2/(7*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**2*tan(c + d*x) + a**2)) + 2*e*sin(c + d*x)/(7*a**2*d*sqrt(e*sec(c + d*x))) + 2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_239():
    f = 1/(sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**2)
    F = 4*I*e**2/(9*d*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*a**2*tan(c + d*x) + a**2)) + 2*e*sin(c + d*x)/(9*a**2*d*(e*sec(c + d*x))**(sympy.S(3)/2)) + 2*elliptic_e(c/2 + d*x/2, 2)/(3*a**2*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_240():
    f = 1/((e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**2)
    F = 4*I*e**2/(11*d*(e*sec(c + d*x))**(sympy.S(7)/2)*(I*a**2*tan(c + d*x) + a**2)) + 2*e*sin(c + d*x)/(11*a**2*d*(e*sec(c + d*x))**(sympy.S(5)/2)) + 10*sin(c + d*x)/(33*a**2*d*e*sqrt(e*sec(c + d*x))) + 10*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(33*a**2*d*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_241():
    f = 1/((e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**2)
    F = 4*I*e**2/(13*d*(e*sec(c + d*x))**(sympy.S(9)/2)*(I*a**2*tan(c + d*x) + a**2)) + 2*e*sin(c + d*x)/(13*a**2*d*(e*sec(c + d*x))**(sympy.S(7)/2)) + 14*sin(c + d*x)/(65*a**2*d*e*(e*sec(c + d*x))**(sympy.S(3)/2)) + 42*elliptic_e(c/2 + d*x/2, 2)/(65*a**2*d*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_242():
    f = 1/((e*sec(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a)**2)
    F = 4*I*e**2/(15*d*(e*sec(c + d*x))**(sympy.S(11)/2)*(I*a**2*tan(c + d*x) + a**2)) + 2*e*sin(c + d*x)/(15*a**2*d*(e*sec(c + d*x))**(sympy.S(9)/2)) + 6*sin(c + d*x)/(35*a**2*d*e*(e*sec(c + d*x))**(sympy.S(5)/2)) + 2*sin(c + d*x)/(7*a**2*d*e**3*sqrt(e*sec(c + d*x))) + 2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(7*a**2*d*e**4)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_243():
    f = (e*sec(c + d*x))**(sympy.S(15)/2)/(I*a*tan(c + d*x) + a)**3
    F = -4*I*e**2*(e*sec(c + d*x))**(sympy.S(11)/2)/(3*a*d*(I*a*tan(c + d*x) + a)**2) - 22*e**8*elliptic_e(c/2 + d*x/2, 2)/(5*a**3*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 22*e**7*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*a**3*d) + 22*e**5*(e*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(15*a**3*d) - 22*I*e**4*(e*sec(c + d*x))**(sympy.S(7)/2)/(21*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_244():
    f = (e*sec(c + d*x))**(sympy.S(13)/2)/(I*a*tan(c + d*x) + a)**3
    F = -4*I*e**2*(e*sec(c + d*x))**(sympy.S(9)/2)/(a*d*(I*a*tan(c + d*x) + a)**2) + 6*e**6*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(a**3*d) + 6*e**5*(e*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(a**3*d) - 18*I*e**4*(e*sec(c + d*x))**(sympy.S(5)/2)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_245():
    f = (e*sec(c + d*x))**(sympy.S(11)/2)/(I*a*tan(c + d*x) + a)**3
    F = 4*I*e**2*(e*sec(c + d*x))**(sympy.S(7)/2)/(a*d*(I*a*tan(c + d*x) + a)**2) + 14*e**6*elliptic_e(c/2 + d*x/2, 2)/(a**3*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 14*e**5*sqrt(e*sec(c + d*x))*sin(c + d*x)/(a**3*d) + 14*I*e**4*(e*sec(c + d*x))**(sympy.S(3)/2)/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_246():
    f = (e*sec(c + d*x))**(sympy.S(9)/2)/(I*a*tan(c + d*x) + a)**3
    F = 4*I*e**2*(e*sec(c + d*x))**(sympy.S(5)/2)/(3*a*d*(I*a*tan(c + d*x) + a)**2) - 10*e**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(3*a**3*d) + 10*I*e**4*sqrt(e*sec(c + d*x))/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_247():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**3
    F = 4*I*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(5*a*d*(I*a*tan(c + d*x) + a)**2) - 6*I*e**4/(5*a**3*d*sqrt(e*sec(c + d*x))) - 6*e**4*elliptic_e(c/2 + d*x/2, 2)/(5*a**3*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_248():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**3
    F = -2*I*e**2*sqrt(e*sec(c + d*x))/(21*d*(I*a**3*tan(c + d*x) + a**3)) + 4*I*e**2*sqrt(e*sec(c + d*x))/(7*a*d*(I*a*tan(c + d*x) + a)**2) - 2*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_249():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**3
    F = 2*I*e**2/(45*d*sqrt(e*sec(c + d*x))*(I*a**3*tan(c + d*x) + a**3)) + 4*I*e**2/(9*a*d*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**2) + 2*e**2*elliptic_e(c/2 + d*x/2, 2)/(15*a**3*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_250():
    f = sqrt(e*sec(c + d*x))/(I*a*tan(c + d*x) + a)**3
    F = 20*I*e**2/(77*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**3*tan(c + d*x) + a**3)) + 2*I*sqrt(e*sec(c + d*x))/(11*d*(I*a*tan(c + d*x) + a)**3) + 10*e*sin(c + d*x)/(77*a**3*d*sqrt(e*sec(c + d*x))) + 10*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_251():
    f = 1/(sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**3)
    F = 28*I*e**2/(117*d*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*a**3*tan(c + d*x) + a**3)) + 2*I/(13*d*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**3) + 14*e*sin(c + d*x)/(117*a**3*d*(e*sec(c + d*x))**(sympy.S(3)/2)) + 14*elliptic_e(c/2 + d*x/2, 2)/(39*a**3*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_252():
    f = 1/((e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**3)
    F = 12*I*e**2/(55*d*(e*sec(c + d*x))**(sympy.S(7)/2)*(I*a**3*tan(c + d*x) + a**3)) + 2*I/(15*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**3) + 6*e*sin(c + d*x)/(55*a**3*d*(e*sec(c + d*x))**(sympy.S(5)/2)) + 2*sin(c + d*x)/(11*a**3*d*e*sqrt(e*sec(c + d*x))) + 2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(11*a**3*d*e**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_253():
    f = (e*sec(c + d*x))**(sympy.S(15)/2)/(I*a*tan(c + d*x) + a)**4
    F = 44*I*e**4*(e*sec(c + d*x))**(sympy.S(7)/2)/(3*d*(I*a**4*tan(c + d*x) + a**4)) + 4*I*e**2*(e*sec(c + d*x))**(sympy.S(11)/2)/(a*d*(I*a*tan(c + d*x) + a)**3) + 154*e**8*elliptic_e(c/2 + d*x/2, 2)/(5*a**4*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) - 154*e**7*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*a**4*d) - 154*e**5*(e*sec(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)/(15*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_254():
    f = (e*sec(c + d*x))**(sympy.S(13)/2)/(I*a*tan(c + d*x) + a)**4
    F = 12*I*e**4*(e*sec(c + d*x))**(sympy.S(5)/2)/(d*(I*a**4*tan(c + d*x) + a**4)) + 4*I*e**2*(e*sec(c + d*x))**(sympy.S(9)/2)/(3*a*d*(I*a*tan(c + d*x) + a)**3) - 10*e**6*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(a**4*d) - 10*e**5*(e*sec(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)/(a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_255():
    f = (e*sec(c + d*x))**(sympy.S(11)/2)/(I*a*tan(c + d*x) + a)**4
    F = -28*I*e**4*(e*sec(c + d*x))**(sympy.S(3)/2)/(5*d*(I*a**4*tan(c + d*x) + a**4)) + 4*I*e**2*(e*sec(c + d*x))**(sympy.S(7)/2)/(5*a*d*(I*a*tan(c + d*x) + a)**3) - 42*e**6*elliptic_e(c/2 + d*x/2, 2)/(5*a**4*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))) + 42*e**5*sqrt(e*sec(c + d*x))*sin(c + d*x)/(5*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_256():
    f = (e*sec(c + d*x))**(sympy.S(9)/2)/(I*a*tan(c + d*x) + a)**4
    F = -20*I*e**4*sqrt(e*sec(c + d*x))/(21*d*(I*a**4*tan(c + d*x) + a**4)) + 4*I*e**2*(e*sec(c + d*x))**(sympy.S(5)/2)/(7*a*d*(I*a*tan(c + d*x) + a)**3) + 10*e**4*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(21*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_257():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**4
    F = -4*I*e**4/(15*d*sqrt(e*sec(c + d*x))*(I*a**4*tan(c + d*x) + a**4)) + 4*I*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(9*a*d*(I*a*tan(c + d*x) + a)**3) - 2*e**4*elliptic_e(c/2 + d*x/2, 2)/(15*a**4*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_258():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**4
    F = -4*I*e**4/(77*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**4*tan(c + d*x) + a**4)) + 4*I*e**2*sqrt(e*sec(c + d*x))/(11*a*d*(I*a*tan(c + d*x) + a)**3) - 2*e**3*sin(c + d*x)/(77*a**4*d*sqrt(e*sec(c + d*x))) - 2*e**2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(77*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_259():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**4
    F = 4*I*e**4/(117*d*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*a**4*tan(c + d*x) + a**4)) + 4*I*e**2/(13*a*d*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**3) + 2*e**3*sin(c + d*x)/(117*a**4*d*(e*sec(c + d*x))**(sympy.S(3)/2)) + 2*e**2*elliptic_e(c/2 + d*x/2, 2)/(39*a**4*d*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_260():
    f = sqrt(e*sec(c + d*x))/(I*a*tan(c + d*x) + a)**4
    F = 4*I*e**2/(33*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a**4*tan(c + d*x) + a**4)) + 2*I*sqrt(e*sec(c + d*x))/(15*d*(I*a*tan(c + d*x) + a)**4) + 14*I*sqrt(e*sec(c + d*x))/(165*a*d*(I*a*tan(c + d*x) + a)**3) + 2*e*sin(c + d*x)/(33*a**4*d*sqrt(e*sec(c + d*x))) + 2*sqrt(e*sec(c + d*x))*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(33*a**4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_261():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)*(I*a*tan(e + f*x) + a)
    F = 6*2**(sympy.S(5)/6)*I*a*(d*sec(e + f*x))**(sympy.S(5)/3)*hyper((sympy.S(-5)/6, sympy.S(5)/6), (sympy.S(11)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(5*f*(I*tan(e + f*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_262():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)
    F = 6*2**(sympy.S(1)/6)*I*a*(d*sec(e + f*x))**(sympy.S(1)/3)*hyper((sympy.S(-1)/6, sympy.S(1)/6), (sympy.S(7)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(f*(I*tan(e + f*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_263():
    f = (I*a*tan(e + f*x) + a)/(d*sec(e + f*x))**(sympy.S(1)/3)
    F = -3*2**(sympy.S(5)/6)*I*a*(I*tan(e + f*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(1)/6), (sympy.S(5)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_264():
    f = (I*a*tan(e + f*x) + a)/(d*sec(e + f*x))**(sympy.S(5)/3)
    F = -3*2**(sympy.S(1)/6)*I*a*(I*tan(e + f*x) + 1)**(sympy.S(5)/6)*hyper((sympy.S(-5)/6, sympy.S(5)/6), (sympy.S(1)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(5*f*(d*sec(e + f*x))**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_265():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)*(I*a*tan(e + f*x) + a)**2
    F = 12*2**(sympy.S(5)/6)*I*a**2*(d*sec(e + f*x))**(sympy.S(5)/3)*hyper((sympy.S(-11)/6, sympy.S(5)/6), (sympy.S(11)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(5*f*(I*tan(e + f*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_266():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**2
    F = 12*2**(sympy.S(1)/6)*I*a**2*(d*sec(e + f*x))**(sympy.S(1)/3)*hyper((sympy.S(-7)/6, sympy.S(1)/6), (sympy.S(7)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(f*(I*tan(e + f*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_267():
    f = (I*a*tan(e + f*x) + a)**2/(d*sec(e + f*x))**(sympy.S(1)/3)
    F = -6*2**(sympy.S(5)/6)*I*(I*a**2*tan(e + f*x) + a**2)*hyper((sympy.S(-5)/6, sympy.S(-1)/6), (sympy.S(5)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(f*(d*sec(e + f*x))**(sympy.S(1)/3)*(I*tan(e + f*x) + 1)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_268():
    f = (I*a*tan(e + f*x) + a)**2/(d*sec(e + f*x))**(sympy.S(5)/3)
    F = -6*2**(sympy.S(1)/6)*I*(I*a**2*tan(e + f*x) + a**2)*hyper((sympy.S(-5)/6, sympy.S(-1)/6), (sympy.S(1)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(5*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(I*tan(e + f*x) + 1)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_269():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)/(I*a*tan(e + f*x) + a)
    F = 3*2**(sympy.S(5)/6)*I*(d*sec(e + f*x))**(sympy.S(5)/3)*(I*tan(e + f*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(5)/6, sympy.S(7)/6), (sympy.S(11)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(10*f*(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_270():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)/(I*a*tan(e + f*x) + a)
    F = 3*2**(sympy.S(1)/6)*I*(d*sec(e + f*x))**(sympy.S(1)/3)*(I*tan(e + f*x) + 1)**(sympy.S(5)/6)*hyper((sympy.S(1)/6, sympy.S(11)/6), (sympy.S(7)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*(I*a*tan(e + f*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_271():
    f = 1/((d*sec(e + f*x))**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a))
    F = -3*2**(sympy.S(5)/6)*I*(I*tan(e + f*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(13)/6), (sympy.S(5)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(4*a*f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_272():
    f = 1/((d*sec(e + f*x))**(sympy.S(5)/3)*(I*a*tan(e + f*x) + a))
    F = -3*2**(sympy.S(1)/6)*I*(I*tan(e + f*x) + 1)**(sympy.S(5)/6)*hyper((sympy.S(-5)/6, sympy.S(17)/6), (sympy.S(1)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(20*a*f*(d*sec(e + f*x))**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_273():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)/(I*a*tan(e + f*x) + a)**2
    F = 3*2**(sympy.S(5)/6)*I*(d*sec(e + f*x))**(sympy.S(5)/3)*(I*tan(e + f*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(5)/6, sympy.S(13)/6), (sympy.S(11)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(20*f*(I*a**2*tan(e + f*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_274():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)/(I*a*tan(e + f*x) + a)**2
    F = 3*2**(sympy.S(1)/6)*I*(d*sec(e + f*x))**(sympy.S(1)/3)*(I*tan(e + f*x) + 1)**(sympy.S(5)/6)*hyper((sympy.S(1)/6, sympy.S(17)/6), (sympy.S(7)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(4*f*(I*a**2*tan(e + f*x) + a**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_275():
    f = 1/((d*sec(e + f*x))**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**2)
    F = -3*2**(sympy.S(5)/6)*I*(I*tan(e + f*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(-1)/6, sympy.S(19)/6), (sympy.S(5)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(8*a**2*f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_276():
    f = 1/((d*sec(e + f*x))**(sympy.S(5)/3)*(I*a*tan(e + f*x) + a)**2)
    F = -3*2**(sympy.S(1)/6)*I*(I*tan(e + f*x) + 1)**(sympy.S(5)/6)*hyper((sympy.S(-5)/6, sympy.S(23)/6), (sympy.S(1)/6,), -I*tan(e + f*x)/2 + sympy.S.Half)/(40*a**2*f*(d*sec(e + f*x))**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_277():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**8
    F = -16*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**4*d) + 24*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**5*d) - 12*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(15)/2)/(15*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_278():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**6
    F = -8*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**3*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_279():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**4
    F = -4*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**2*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_280():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**2
    F = -2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_281():
    f = sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**2
    F = -3*sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(8*d) - I*a**2/(2*d*(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 3*I*a/(4*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_282():
    f = sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**4
    F = -35*sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(128*d) - I*a**4/(4*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 7*I*a**3/(16*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 35*I*a**2/(96*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 35*I*a/(64*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_283():
    f = sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**6
    F = -231*sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(1024*d) - I*a**6/(6*d*(-I*a*tan(c + d*x) + a)**3*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) - 11*I*a**5/(48*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) - 33*I*a**4/(64*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 231*I*a**3/(640*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 77*I*a**2/(256*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 231*I*a/(512*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_284():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**7
    F = 256*I*a**4*sec(c + d*x)**7/(3003*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 64*I*a**3*sec(c + d*x)**7/(429*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 24*I*a**2*sec(c + d*x)**7/(143*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 2*I*a*sec(c + d*x)**7/(13*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_285():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**5
    F = 64*I*a**3*sec(c + d*x)**5/(315*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 16*I*a**2*sec(c + d*x)**5/(63*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 2*I*a*sec(c + d*x)**5/(9*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_286():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**3
    F = 8*I*a**2*sec(c + d*x)**3/(15*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 2*I*a*sec(c + d*x)**3/(5*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_287():
    f = sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)
    F = 2*I*a*sec(c + d*x)/(d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_288():
    f = sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)
    F = sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(2*d) - I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_289():
    f = sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3
    F = 5*sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(16*d) + 5*I*a*cos(c + d*x)/(12*d*sqrt(I*a*tan(c + d*x) + a)) - I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(3*d) - 5*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_290():
    f = sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**5
    F = 63*sqrt(2)*I*sqrt(a)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(256*d) + 9*I*a*cos(c + d*x)**3/(40*d*sqrt(I*a*tan(c + d*x) + a)) + 21*I*a*cos(c + d*x)/(64*d*sqrt(I*a*tan(c + d*x) + a)) - I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**5/(5*d) - 21*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(80*d) - 63*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_291():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**8
    F = -16*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**4*d) + 24*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**5*d) - 4*I*(I*a*tan(c + d*x) + a)**(sympy.S(15)/2)/(5*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(17)/2)/(17*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_292():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**6
    F = -8*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**3*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_293():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**4
    F = -4*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**2*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_294():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**2
    F = -2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_295():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**2
    F = -sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4*d) - I*a**2*sqrt(I*a*tan(c + d*x) + a)/(2*d*(-I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_296():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**4
    F = -15*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(64*d) - I*a**4/(4*d*(-I*a*tan(c + d*x) + a)**2*sqrt(I*a*tan(c + d*x) + a)) - 5*I*a**3/(16*d*(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 15*I*a**2/(32*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_297():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**6
    F = -105*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(512*d) - I*a**6/(6*d*(-I*a*tan(c + d*x) + a)**3*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 3*I*a**5/(16*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 21*I*a**4/(64*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 35*I*a**3/(128*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 105*I*a**2/(256*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_298():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**5
    F = 256*I*a**4*sec(c + d*x)**5/(1155*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 64*I*a**3*sec(c + d*x)**5/(231*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 8*I*a**2*sec(c + d*x)**5/(33*d*sqrt(I*a*tan(c + d*x) + a)) + 2*I*a*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**5/(11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_299():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3
    F = 64*I*a**3*sec(c + d*x)**3/(105*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 16*I*a**2*sec(c + d*x)**3/(35*d*sqrt(I*a*tan(c + d*x) + a)) + 2*I*a*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**3/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_300():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)
    F = 8*I*a**2*sec(c + d*x)/(3*d*sqrt(I*a*tan(c + d*x) + a)) + 2*I*a*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_301():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)
    F = -2*I*a*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_302():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3
    F = sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(4*d) - I*a*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(2*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_303():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**5
    F = 7*sqrt(2)*I*a**(sympy.S(3)/2)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(32*d) + 7*I*a**2*cos(c + d*x)/(24*d*sqrt(I*a*tan(c + d*x) + a)) - 7*I*a*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(30*d) - 7*I*a*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(16*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_304():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**8
    F = -16*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**4*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(15)/2)/(5*a**5*d) - 12*I*(I*a*tan(c + d*x) + a)**(sympy.S(17)/2)/(17*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(19)/2)/(19*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_305():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**6
    F = -8*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**3*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(15)/2)/(15*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_306():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**4
    F = -4*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**2*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_307():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**2
    F = -2*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_308():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**2
    F = sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2*d) - I*a**3*sqrt(I*a*tan(c + d*x) + a)/(d*(-I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_309():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**4
    F = -3*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(32*d) - I*a**4*sqrt(I*a*tan(c + d*x) + a)/(4*d*(-I*a*tan(c + d*x) + a)**2) - 3*I*a**3*sqrt(I*a*tan(c + d*x) + a)/(16*d*(-I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_310():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**6
    F = -35*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(256*d) - I*a**6/(6*d*(-I*a*tan(c + d*x) + a)**3*sqrt(I*a*tan(c + d*x) + a)) - 7*I*a**5/(48*d*(-I*a*tan(c + d*x) + a)**2*sqrt(I*a*tan(c + d*x) + a)) - 35*I*a**4/(192*d*(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 35*I*a**3/(128*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_311():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)**3
    F = 256*I*a**4*sec(c + d*x)**3/(315*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 64*I*a**3*sec(c + d*x)**3/(105*d*sqrt(I*a*tan(c + d*x) + a)) + 8*I*a**2*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**3/(21*d) + 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)**3/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_312():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)
    F = 64*I*a**3*sec(c + d*x)/(15*d*sqrt(I*a*tan(c + d*x) + a)) + 16*I*a**2*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)/(15*d) + 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_313():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)
    F = -8*I*a**2*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/d + 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_314():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3
    F = -2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_315():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5
    F = sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(8*d) - I*a**2*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(4*d) - I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3/(6*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_316():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**7
    F = 9*sqrt(2)*I*a**(sympy.S(5)/2)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(64*d) + 3*I*a**3*cos(c + d*x)/(16*d*sqrt(I*a*tan(c + d*x) + a)) - 3*I*a**2*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(20*d) - 9*I*a**2*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(32*d) - 9*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**5/(70*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_317():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**8
    F = -16*I*(I*a*tan(c + d*x) + a)**(sympy.S(15)/2)/(15*a**4*d) + 24*I*(I*a*tan(c + d*x) + a)**(sympy.S(17)/2)/(17*a**5*d) - 12*I*(I*a*tan(c + d*x) + a)**(sympy.S(19)/2)/(19*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(21)/2)/(21*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_318():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**6
    F = -8*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**3*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(15)/2)/(15*a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(17)/2)/(17*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_319():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**4
    F = -4*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**2*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_320():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)**2
    F = -2*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_321():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**2
    F = 3*sqrt(2)*I*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/d - 3*I*a**3*sqrt(I*a*tan(c + d*x) + a)/d - I*a**3*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(d*(-I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_322():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**4
    F = sqrt(2)*I*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(16*d) - I*a**5*sqrt(I*a*tan(c + d*x) + a)/(2*d*(-I*a*tan(c + d*x) + a)**2) + I*a**4*sqrt(I*a*tan(c + d*x) + a)/(8*d*(-I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_323():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**6
    F = -5*sqrt(2)*I*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(128*d) - I*a**6*sqrt(I*a*tan(c + d*x) + a)/(6*d*(-I*a*tan(c + d*x) + a)**3) - 5*I*a**5*sqrt(I*a*tan(c + d*x) + a)/(48*d*(-I*a*tan(c + d*x) + a)**2) - 5*I*a**4*sqrt(I*a*tan(c + d*x) + a)/(64*d*(-I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_324():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*sec(c + d*x)
    F = 256*I*a**4*sec(c + d*x)/(35*d*sqrt(I*a*tan(c + d*x) + a)) + 64*I*a**3*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)/(35*d) + 24*I*a**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*sec(c + d*x)/(35*d) + 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*sec(c + d*x)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_325():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)
    F = -64*I*a**3*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(3*d) + 16*I*a**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)/(3*d) + 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_326():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**3
    F = 8*I*a**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3/(3*d) - 2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**3/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_327():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**5
    F = -2*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_328():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**7
    F = sqrt(2)*I*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(16*d) - I*a**3*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(8*d) - I*a**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**3/(12*d) - I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**5/(10*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**7/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_329():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**9
    F = 11*sqrt(2)*I*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(128*d) + 11*I*a**4*cos(c + d*x)/(96*d*sqrt(I*a*tan(c + d*x) + a)) - 11*I*a**3*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(120*d) - 11*I*a**3*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(64*d) - 11*I*a**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**5/(140*d) - 11*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**7/(126*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**9/(9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_330():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**11
    F = 195*sqrt(2)*I*a**(sympy.S(7)/2)*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(2048*d) + 39*I*a**4*cos(c + d*x)**3/(448*d*sqrt(I*a*tan(c + d*x) + a)) + 65*I*a**4*cos(c + d*x)/(512*d*sqrt(I*a*tan(c + d*x) + a)) - 13*I*a**3*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**5/(168*d) - 13*I*a**3*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(128*d) - 195*I*a**3*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(1024*d) - 65*I*a**2*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)*cos(c + d*x)**7/(924*d) - 5*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)*cos(c + d*x)**9/(66*d) - I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)*cos(c + d*x)**11/(11*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_331():
    f = sec(c + d*x)**8/sqrt(I*a*tan(c + d*x) + a)
    F = -16*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**4*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(3*a**5*d) - 12*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_332():
    f = sec(c + d*x)**6/sqrt(I*a*tan(c + d*x) + a)
    F = -8*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**3*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_333():
    f = sec(c + d*x)**4/sqrt(I*a*tan(c + d*x) + a)
    F = -4*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**2*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_334():
    f = sec(c + d*x)**2/sqrt(I*a*tan(c + d*x) + a)
    F = -2*I*sqrt(I*a*tan(c + d*x) + a)/(a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_335():
    f = cos(c + d*x)**2/sqrt(I*a*tan(c + d*x) + a)
    F = -I*a**2/(2*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 5*I*a/(12*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 5*I/(8*d*sqrt(I*a*tan(c + d*x) + a)) - 5*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(16*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_336():
    f = cos(c + d*x)**4/sqrt(I*a*tan(c + d*x) + a)
    F = -I*a**4/(4*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) - 9*I*a**3/(16*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 63*I*a**2/(160*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 21*I*a/(64*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 63*I/(128*d*sqrt(I*a*tan(c + d*x) + a)) - 63*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(256*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_337():
    f = cos(c + d*x)**6/sqrt(I*a*tan(c + d*x) + a)
    F = -I*a**6/(6*d*(-I*a*tan(c + d*x) + a)**3*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) - 13*I*a**5/(48*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) - 143*I*a**4/(192*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 429*I*a**3/(896*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 429*I*a**2/(1280*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 143*I*a/(512*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 429*I/(1024*d*sqrt(I*a*tan(c + d*x) + a)) - 429*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2048*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_338():
    f = sec(c + d*x)**9/sqrt(I*a*tan(c + d*x) + a)
    F = 256*I*a**4*sec(c + d*x)**9/(6435*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 64*I*a**3*sec(c + d*x)**9/(715*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 8*I*a**2*sec(c + d*x)**9/(65*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 2*I*a*sec(c + d*x)**9/(15*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_339():
    f = sec(c + d*x)**7/sqrt(I*a*tan(c + d*x) + a)
    F = 64*I*a**3*sec(c + d*x)**7/(693*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 16*I*a**2*sec(c + d*x)**7/(99*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 2*I*a*sec(c + d*x)**7/(11*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_340():
    f = sec(c + d*x)**5/sqrt(I*a*tan(c + d*x) + a)
    F = 8*I*a**2*sec(c + d*x)**5/(35*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 2*I*a*sec(c + d*x)**5/(7*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_341():
    f = sec(c + d*x)**3/sqrt(I*a*tan(c + d*x) + a)
    F = 2*I*a*sec(c + d*x)**3/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_342():
    f = sec(c + d*x)/sqrt(I*a*tan(c + d*x) + a)
    F = sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_343():
    f = cos(c + d*x)/sqrt(I*a*tan(c + d*x) + a)
    F = I*cos(c + d*x)/(2*d*sqrt(I*a*tan(c + d*x) + a)) - 3*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(4*a*d) + 3*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(8*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_344():
    f = cos(c + d*x)**3/sqrt(I*a*tan(c + d*x) + a)
    F = I*cos(c + d*x)**3/(4*d*sqrt(I*a*tan(c + d*x) + a)) + 35*I*cos(c + d*x)/(96*d*sqrt(I*a*tan(c + d*x) + a)) - 7*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(24*a*d) - 35*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(64*a*d) + 35*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(128*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_345():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -16*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**4*d) + 24*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**5*d) - 4*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(3*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_346():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -8*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**3*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_347():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -4*I*sqrt(I*a*tan(c + d*x) + a)/(a**2*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_348():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*I/(a*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_349():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -I*a**2/(2*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 7*I*a/(20*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 7*I/(24*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 7*I/(16*a*d*sqrt(I*a*tan(c + d*x) + a)) - 7*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(32*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_350():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -I*a**4/(4*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) - 11*I*a**3/(16*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 99*I*a**2/(224*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 99*I*a/(320*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 33*I/(128*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 99*I/(256*a*d*sqrt(I*a*tan(c + d*x) + a)) - 99*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(512*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_351():
    f = cos(c + d*x)**6/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -I*a**6/(6*d*(-I*a*tan(c + d*x) + a)**3*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) - 5*I*a**5/(16*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) - 65*I*a**4/(64*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 715*I*a**3/(1152*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 715*I*a**2/(1792*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 143*I*a/(512*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 715*I/(3072*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 715*I/(2048*a*d*sqrt(I*a*tan(c + d*x) + a)) - 715*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(4096*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_352():
    f = sec(c + d*x)**11/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 256*I*a**4*sec(c + d*x)**11/(12155*d*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) + 64*I*a**3*sec(c + d*x)**11/(1105*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 8*I*a**2*sec(c + d*x)**11/(85*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 2*I*a*sec(c + d*x)**11/(17*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_353():
    f = sec(c + d*x)**9/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 64*I*a**3*sec(c + d*x)**9/(1287*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 16*I*a**2*sec(c + d*x)**9/(143*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 2*I*a*sec(c + d*x)**9/(13*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_354():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 8*I*a**2*sec(c + d*x)**7/(63*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 2*I*a*sec(c + d*x)**7/(9*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_355():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*I*a*sec(c + d*x)**5/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_356():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -2*I*sec(c + d*x)/(a*d*sqrt(I*a*tan(c + d*x) + a)) + 2*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_357():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = I*sec(c + d*x)/(2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(4*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_358():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = I*cos(c + d*x)/(4*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 5*I*cos(c + d*x)/(16*a*d*sqrt(I*a*tan(c + d*x) + a)) - 15*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(32*a**2*d) + 15*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(64*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_359():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = I*cos(c + d*x)**3/(6*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 3*I*cos(c + d*x)**3/(16*a*d*sqrt(I*a*tan(c + d*x) + a)) + 35*I*cos(c + d*x)/(128*a*d*sqrt(I*a*tan(c + d*x) + a)) - 7*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(32*a**2*d) - 105*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(256*a**2*d) + 105*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(512*a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_360():
    f = sec(c + d*x)**10/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -32*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**5*d) + 64*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**6*d) - 16*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(3*a**7*d) + 16*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**8*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)/(13*a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_361():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -16*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**4*d) + 24*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**5*d) - 12*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_362():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -8*I*sqrt(I*a*tan(c + d*x) + a)/(a**3*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_363():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 4*I/(a**2*d*sqrt(I*a*tan(c + d*x) + a)) + 2*I*sqrt(I*a*tan(c + d*x) + a)/(a**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_364():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*I/(3*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_365():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -I*a**2/(2*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 9*I*a/(28*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 9*I/(40*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 3*I/(16*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 9*I/(32*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 9*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(64*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_366():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -I*a**4/(4*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) - 13*I*a**3/(16*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 143*I*a**2/(288*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 143*I*a/(448*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 143*I/(640*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 143*I/(768*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 143*I/(512*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 143*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(1024*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_367():
    f = sec(c + d*x)**13/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 256*I*a**4*sec(c + d*x)**13/(20995*d*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)) + 64*I*a**3*sec(c + d*x)**13/(1615*d*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) + 24*I*a**2*sec(c + d*x)**13/(323*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 2*I*a*sec(c + d*x)**13/(19*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_368():
    f = sec(c + d*x)**11/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 64*I*a**3*sec(c + d*x)**11/(2145*d*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) + 16*I*a**2*sec(c + d*x)**11/(195*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 2*I*a*sec(c + d*x)**11/(15*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_369():
    f = sec(c + d*x)**9/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 8*I*a**2*sec(c + d*x)**9/(99*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 2*I*a*sec(c + d*x)**9/(11*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_370():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*I*a*sec(c + d*x)**7/(7*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_371():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -2*I*sec(c + d*x)**3/(3*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 4*I*sec(c + d*x)/(a**2*d*sqrt(I*a*tan(c + d*x) + a)) + 4*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_372():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = I*sec(c + d*x)/(a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_373():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = I*sec(c + d*x)/(4*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 3*I*sec(c + d*x)/(16*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 3*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(32*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_374():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = I*cos(c + d*x)/(6*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 7*I*cos(c + d*x)/(48*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 35*I*cos(c + d*x)/(192*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 35*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(128*a**3*d) + 35*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(256*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_375():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = I*cos(c + d*x)**3/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 11*I*cos(c + d*x)**3/(96*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 33*I*cos(c + d*x)**3/(256*a**2*d*sqrt(I*a*tan(c + d*x) + a)) + 385*I*cos(c + d*x)/(2048*a**2*d*sqrt(I*a*tan(c + d*x) + a)) - 77*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(512*a**3*d) - 1155*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(4096*a**3*d) + 1155*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(8192*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_376():
    f = sec(c + d*x)**10/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = -32*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**5*d) + 64*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**6*d) - 48*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**7*d) + 16*I*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)/(9*a**8*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)/(11*a**9*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_377():
    f = sec(c + d*x)**8/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = -16*I*sqrt(I*a*tan(c + d*x) + a)/(a**4*d) + 8*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(a**5*d) - 12*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*a**6*d) + 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)/(7*a**7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_378():
    f = sec(c + d*x)**6/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = 8*I/(a**3*d*sqrt(I*a*tan(c + d*x) + a)) + 8*I*sqrt(I*a*tan(c + d*x) + a)/(a**4*d) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*a**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_379():
    f = sec(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = 4*I/(3*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 2*I/(a**3*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_380():
    f = sec(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = 2*I/(5*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_381():
    f = cos(c + d*x)**2/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = -I*a**2/(2*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 11*I*a/(36*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 11*I/(56*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 11*I/(80*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 11*I/(96*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 11*I/(64*a**3*d*sqrt(I*a*tan(c + d*x) + a)) - 11*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_382():
    f = cos(c + d*x)**4/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = -I*a**4/(4*d*(-I*a*tan(c + d*x) + a)**2*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) - 15*I*a**3/(16*d*(-I*a*tan(c + d*x) + a)*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) + 195*I*a**2/(352*d*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) + 65*I*a/(192*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2)) + 195*I/(896*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 39*I/(256*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 65*I/(512*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 195*I/(1024*a**3*d*sqrt(I*a*tan(c + d*x) + a)) - 195*sqrt(2)*I*atanh(sqrt(2)*sqrt(I*a*tan(c + d*x) + a)/(2*sqrt(a)))/(2048*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_383():
    f = sec(c + d*x)**13/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = 64*I*a**3*sec(c + d*x)**13/(3315*d*(I*a*tan(c + d*x) + a)**(sympy.S(13)/2)) + 16*I*a**2*sec(c + d*x)**13/(255*d*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) + 2*I*a*sec(c + d*x)**13/(17*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_384():
    f = sec(c + d*x)**11/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = 8*I*a**2*sec(c + d*x)**11/(143*d*(I*a*tan(c + d*x) + a)**(sympy.S(11)/2)) + 2*I*a*sec(c + d*x)**11/(13*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_385():
    f = sec(c + d*x)**9/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = 2*I*a*sec(c + d*x)**9/(9*d*(I*a*tan(c + d*x) + a)**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_386():
    f = sec(c + d*x)**7/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = -2*I*sec(c + d*x)**5/(5*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) - 4*I*sec(c + d*x)**3/(3*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 8*I*sec(c + d*x)/(a**3*d*sqrt(I*a*tan(c + d*x) + a)) + 8*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_387():
    f = sec(c + d*x)**5/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = -2*I*sec(c + d*x)**3/(a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 6*I*sec(c + d*x)/(a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - 3*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_388():
    f = sec(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = I*sec(c + d*x)/(2*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) - I*sec(c + d*x)/(8*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(16*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_389():
    f = sec(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = I*sec(c + d*x)/(6*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 5*I*sec(c + d*x)/(48*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 5*I*sec(c + d*x)/(64*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 5*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(128*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_390():
    f = cos(c + d*x)/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = I*cos(c + d*x)/(8*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 3*I*cos(c + d*x)/(32*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 21*I*cos(c + d*x)/(256*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 105*I*cos(c + d*x)/(1024*a**3*d*sqrt(I*a*tan(c + d*x) + a)) - 315*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(2048*a**4*d) + 315*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(4096*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_391():
    f = cos(c + d*x)**3/(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = I*cos(c + d*x)**3/(10*d*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)) + 13*I*cos(c + d*x)**3/(160*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 143*I*cos(c + d*x)**3/(1920*a**2*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 429*I*cos(c + d*x)**3/(5120*a**3*d*sqrt(I*a*tan(c + d*x) + a)) + 1001*I*cos(c + d*x)/(8192*a**3*d*sqrt(I*a*tan(c + d*x) + a)) - 1001*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**3/(10240*a**4*d) - 3003*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)/(16384*a**4*d) + 3003*sqrt(2)*I*atanh(sqrt(2)*sqrt(a)*sec(c + d*x)/(2*sqrt(I*a*tan(c + d*x) + a)))/(32768*a**(sympy.S(7)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_392():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)
    F = sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(4*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(3)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(4*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(2*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(2*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + I*a*(e*sec(c + d*x))**(sympy.S(3)/2)/(d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_393():
    f = sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)
    F = -sqrt(2)*I*sqrt(a)*sqrt(e)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*d) + sqrt(2)*I*sqrt(a)*sqrt(e)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*d) + sqrt(2)*I*sqrt(a)*sqrt(e)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/d - sqrt(2)*I*sqrt(a)*sqrt(e)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_394():
    f = sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x))
    F = -2*I*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_395():
    f = sqrt(I*a*tan(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = 4*I*a*sqrt(e*sec(c + d*x))/(3*d*e**2*sqrt(I*a*tan(c + d*x) + a)) - 2*I*sqrt(I*a*tan(c + d*x) + a)/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_396():
    f = sqrt(I*a*tan(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = 8*I*a/(15*d*e**2*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)) - 2*I*sqrt(I*a*tan(c + d*x) + a)/(5*d*(e*sec(c + d*x))**(sympy.S(5)/2)) - 16*I*sqrt(I*a*tan(c + d*x) + a)/(15*d*e**2*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_397():
    f = sqrt(I*a*tan(c + d*x) + a)/(e*sec(c + d*x))**(sympy.S(7)/2)
    F = 12*I*a/(35*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)) + 32*I*a*sqrt(e*sec(c + d*x))/(35*d*e**4*sqrt(I*a*tan(c + d*x) + a)) - 2*I*sqrt(I*a*tan(c + d*x) + a)/(7*d*(e*sec(c + d*x))**(sympy.S(7)/2)) - 16*I*sqrt(I*a*tan(c + d*x) + a)/(35*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_398():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -7*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(32*d) + 7*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(5)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(32*d) + 7*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(16*d) - 7*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(16*d) + 7*I*a**2*(e*sec(c + d*x))**(sympy.S(5)/2)/(12*d*sqrt(I*a*tan(c + d*x) + a)) - 7*I*a*e**2*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/(8*d) + I*a*(e*sec(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_399():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 5*sqrt(2)*I*a**(sympy.S(5)/2)*e**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(16*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 5*sqrt(2)*I*a**(sympy.S(5)/2)*e**(sympy.S(3)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(16*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 5*sqrt(2)*I*a**(sympy.S(5)/2)*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(8*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 5*sqrt(2)*I*a**(sympy.S(5)/2)*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(8*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 5*I*a**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(4*d*sqrt(I*a*tan(c + d*x) + a)) + I*a*(e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_400():
    f = sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -3*sqrt(2)*I*a**(sympy.S(3)/2)*sqrt(e)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*d) + 3*sqrt(2)*I*a**(sympy.S(3)/2)*sqrt(e)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*d) + 3*sqrt(2)*I*a**(sympy.S(3)/2)*sqrt(e)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*d) - 3*sqrt(2)*I*a**(sympy.S(3)/2)*sqrt(e)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*d) + I*a*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_401():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/sqrt(e*sec(c + d*x))
    F = -sqrt(2)*I*a**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*a**(sympy.S(5)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*a**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*a**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 4*I*a*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_402():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = -2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_403():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = -4*I*a*sqrt(I*a*tan(c + d*x) + a)/(5*d*e**2*sqrt(e*sec(c + d*x))) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(5*d*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_404():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(e*sec(c + d*x))**(sympy.S(7)/2)
    F = 16*I*a**2*sqrt(e*sec(c + d*x))/(21*d*e**4*sqrt(I*a*tan(c + d*x) + a)) - 8*I*a*sqrt(I*a*tan(c + d*x) + a)/(21*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(7*d*(e*sec(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_405():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(e*sec(c + d*x))**(sympy.S(9)/2)
    F = 16*I*a**2/(45*d*e**4*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)) - 4*I*a*sqrt(I*a*tan(c + d*x) + a)/(15*d*e**2*(e*sec(c + d*x))**(sympy.S(5)/2)) - 32*I*a*sqrt(I*a*tan(c + d*x) + a)/(45*d*e**4*sqrt(e*sec(c + d*x))) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(9*d*(e*sec(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_406():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 15*sqrt(2)*I*a**(sympy.S(7)/2)*e**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(32*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 15*sqrt(2)*I*a**(sympy.S(7)/2)*e**(sympy.S(3)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(32*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 15*sqrt(2)*I*a**(sympy.S(7)/2)*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(16*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 15*sqrt(2)*I*a**(sympy.S(7)/2)*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(16*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 15*I*a**3*(e*sec(c + d*x))**(sympy.S(3)/2)/(8*d*sqrt(I*a*tan(c + d*x) + a)) + 3*I*a**2*(e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)/(4*d) + I*a*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_407():
    f = sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = -21*sqrt(2)*I*a**(sympy.S(5)/2)*sqrt(e)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(16*d) + 21*sqrt(2)*I*a**(sympy.S(5)/2)*sqrt(e)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(16*d) + 21*sqrt(2)*I*a**(sympy.S(5)/2)*sqrt(e)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(8*d) - 21*sqrt(2)*I*a**(sympy.S(5)/2)*sqrt(e)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(8*d) + 7*I*a**2*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/(4*d) + I*a*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_408():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/sqrt(e*sec(c + d*x))
    F = -5*sqrt(2)*I*a**(sympy.S(7)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(4*d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 5*sqrt(2)*I*a**(sympy.S(7)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(4*d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 5*sqrt(2)*I*a**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(2*d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 5*sqrt(2)*I*a**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(2*d*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 10*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(d*sqrt(e*sec(c + d*x))) + I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_409():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = sqrt(2)*I*a**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*d*e**(sympy.S(3)/2)) - sqrt(2)*I*a**(sympy.S(5)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*d*e**(sympy.S(3)/2)) - sqrt(2)*I*a**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(d*e**(sympy.S(3)/2)) + sqrt(2)*I*a**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(d*e**(sympy.S(3)/2)) - 4*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_410():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = -2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(5*d*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_411():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(e*sec(c + d*x))**(sympy.S(7)/2)
    F = -4*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(21*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(7*d*(e*sec(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_412():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(e*sec(c + d*x))**(sympy.S(9)/2)
    F = -16*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(45*d*e**4*sqrt(e*sec(c + d*x))) - 8*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(45*d*e**2*(e*sec(c + d*x))**(sympy.S(5)/2)) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(9*d*(e*sec(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_413():
    f = (I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(e*sec(c + d*x))**(sympy.S(11)/2)
    F = 32*I*a**3*sqrt(e*sec(c + d*x))/(77*d*e**6*sqrt(I*a*tan(c + d*x) + a)) - 16*I*a**2*sqrt(I*a*tan(c + d*x) + a)/(77*d*e**4*(e*sec(c + d*x))**(sympy.S(3)/2)) - 12*I*a*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)/(77*d*e**2*(e*sec(c + d*x))**(sympy.S(7)/2)) - 2*I*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)/(11*d*(e*sec(c + d*x))**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_414():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = -I*e**2*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/(a*d) - sqrt(2)*I*e**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*sqrt(a)*d) + sqrt(2)*I*e**(sympy.S(5)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*sqrt(a)*d) + sqrt(2)*I*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*sqrt(a)*d) - sqrt(2)*I*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*sqrt(a)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_415():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = sqrt(2)*I*sqrt(a)*e**(sympy.S(3)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*sqrt(a)*e**(sympy.S(3)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*sqrt(a)*e**(sympy.S(3)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*sqrt(a)*e**(sympy.S(3)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_416():
    f = sqrt(e*sec(c + d*x))/sqrt(I*a*tan(c + d*x) + a)
    F = 2*I*sqrt(e*sec(c + d*x))/(d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_417():
    f = 1/(sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a))
    F = 2*I/(3*d*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)) - 4*I*sqrt(I*a*tan(c + d*x) + a)/(3*a*d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_418():
    f = 1/((e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a))
    F = 2*I/(5*d*(e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)) + 16*I*sqrt(e*sec(c + d*x))/(15*d*e**2*sqrt(I*a*tan(c + d*x) + a)) - 8*I*sqrt(I*a*tan(c + d*x) + a)/(15*a*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_419():
    f = 1/((e*sec(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a))
    F = 2*I/(7*d*(e*sec(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)) + 16*I/(35*d*e**2*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)) - 12*I*sqrt(I*a*tan(c + d*x) + a)/(35*a*d*(e*sec(c + d*x))**(sympy.S(5)/2)) - 32*I*sqrt(I*a*tan(c + d*x) + a)/(35*a*d*e**2*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_420():
    f = 1/((e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a))
    F = 2*I/(9*d*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a)) + 32*I/(105*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)) + 256*I*sqrt(e*sec(c + d*x))/(315*d*e**4*sqrt(I*a*tan(c + d*x) + a)) - 16*I*sqrt(I*a*tan(c + d*x) + a)/(63*a*d*(e*sec(c + d*x))**(sympy.S(7)/2)) - 128*I*sqrt(I*a*tan(c + d*x) + a)/(315*a*d*e**2*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_421():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = -I*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(a*d*sqrt(I*a*tan(c + d*x) + a)) + 3*sqrt(2)*I*e**(sympy.S(7)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(4*sqrt(a)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 3*sqrt(2)*I*e**(sympy.S(7)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(4*sqrt(a)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 3*sqrt(2)*I*e**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(2*sqrt(a)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 3*sqrt(2)*I*e**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(2*sqrt(a)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_422():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 4*I*e**2*sqrt(e*sec(c + d*x))/(a*d*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*e**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*a**(sympy.S(3)/2)*d) - sqrt(2)*I*e**(sympy.S(5)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*a**(sympy.S(3)/2)*d) - sqrt(2)*I*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(a**(sympy.S(3)/2)*d) + sqrt(2)*I*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(a**(sympy.S(3)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_423():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*I*(e*sec(c + d*x))**(sympy.S(3)/2)/(3*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_424():
    f = sqrt(e*sec(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2*I*sqrt(e*sec(c + d*x))/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 4*I*sqrt(e*sec(c + d*x))/(5*a*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_425():
    f = 1/(sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    F = 2*I/(7*d*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 8*I/(21*a*d*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)) - 16*I*sqrt(I*a*tan(c + d*x) + a)/(21*a**2*d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_426():
    f = 1/((e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    F = 2*I/(9*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 4*I/(15*a*d*(e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)) + 32*I*sqrt(e*sec(c + d*x))/(45*a*d*e**2*sqrt(I*a*tan(c + d*x) + a)) - 16*I*sqrt(I*a*tan(c + d*x) + a)/(45*a**2*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_427():
    f = 1/((e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    F = 2*I/(11*d*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 16*I/(77*a*d*(e*sec(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)) + 128*I/(385*a*d*e**2*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)) - 96*I*sqrt(I*a*tan(c + d*x) + a)/(385*a**2*d*(e*sec(c + d*x))**(sympy.S(5)/2)) - 256*I*sqrt(I*a*tan(c + d*x) + a)/(385*a**2*d*e**2*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_428():
    f = (e*sec(c + d*x))**(sympy.S(9)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 4*I*e**2*(e*sec(c + d*x))**(sympy.S(5)/2)/(a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 5*I*e**4*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/(a**3*d) + 5*sqrt(2)*I*e**(sympy.S(9)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*a**(sympy.S(5)/2)*d) - 5*sqrt(2)*I*e**(sympy.S(9)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*a**(sympy.S(5)/2)*d) - 5*sqrt(2)*I*e**(sympy.S(9)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*a**(sympy.S(5)/2)*d) + 5*sqrt(2)*I*e**(sympy.S(9)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*a**(sympy.S(5)/2)*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_429():
    f = (e*sec(c + d*x))**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 4*I*e**2*(e*sec(c + d*x))**(sympy.S(3)/2)/(3*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) - sqrt(2)*I*e**(sympy.S(7)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*a**(sympy.S(3)/2)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*e**(sympy.S(7)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*a**(sympy.S(3)/2)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*e**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(a**(sympy.S(3)/2)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*e**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(a**(sympy.S(3)/2)*d*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_430():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*I*(e*sec(c + d*x))**(sympy.S(5)/2)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_431():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*I*(e*sec(c + d*x))**(sympy.S(3)/2)/(7*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 4*I*(e*sec(c + d*x))**(sympy.S(3)/2)/(21*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_432():
    f = sqrt(e*sec(c + d*x))/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2*I*sqrt(e*sec(c + d*x))/(9*d*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 8*I*sqrt(e*sec(c + d*x))/(45*a*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 16*I*sqrt(e*sec(c + d*x))/(45*a**2*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_433():
    f = 1/(sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    F = 2*I/(11*d*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 12*I/(77*a*d*sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 16*I/(77*a**2*d*sqrt(e*sec(c + d*x))*sqrt(I*a*tan(c + d*x) + a)) - 32*I*sqrt(I*a*tan(c + d*x) + a)/(77*a**3*d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_434():
    f = 1/((e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2))
    F = 2*I/(13*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)) + 16*I/(117*a*d*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)) + 32*I/(195*a**2*d*(e*sec(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)) + 256*I*sqrt(e*sec(c + d*x))/(585*a**2*d*e**2*sqrt(I*a*tan(c + d*x) + a)) - 128*I*sqrt(I*a*tan(c + d*x) + a)/(585*a**3*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_435():
    f = (e*sec(c + d*x))**(sympy.S(7)/3)/sqrt(I*a*tan(c + d*x) + a)
    F = 3*2**(sympy.S(2)/3)*I*a*(e*sec(c + d*x))**(sympy.S(7)/3)*(I*tan(c + d*x) + 1)**(sympy.S(1)/3)*hyper((sympy.S(1)/3, sympy.S(7)/6), (sympy.S(13)/6,), -I*tan(c + d*x)/2 + sympy.S.Half)/(7*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_436():
    f = (e*sec(c + d*x))**(sympy.S(5)/3)/sqrt(I*a*tan(c + d*x) + a)
    F = 3*2**(sympy.S(1)/3)*I*a*(e*sec(c + d*x))**(sympy.S(5)/3)*(I*tan(c + d*x) + 1)**(sympy.S(2)/3)*hyper((sympy.S(2)/3, sympy.S(5)/6), (sympy.S(11)/6,), -I*tan(c + d*x)/2 + sympy.S.Half)/(5*d*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_437():
    f = (e*sec(c + d*x))**(sympy.S(2)/3)/sqrt(I*a*tan(c + d*x) + a)
    F = 3*2**(sympy.S(5)/6)*I*(e*sec(c + d*x))**(sympy.S(2)/3)*(I*tan(c + d*x) + 1)**(sympy.S(1)/6)*hyper((sympy.S(1)/3, sympy.S(7)/6), (sympy.S(4)/3,), -I*tan(c + d*x)/2 + sympy.S.Half)/(4*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_438():
    f = (e*sec(c + d*x))**(sympy.S(1)/3)/sqrt(I*a*tan(c + d*x) + a)
    F = 3*2**(sympy.S(2)/3)*I*(e*sec(c + d*x))**(sympy.S(1)/3)*(I*tan(c + d*x) + 1)**(sympy.S(1)/3)*hyper((sympy.S(1)/6, sympy.S(4)/3), (sympy.S(7)/6,), -I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_439():
    f = 1/((e*sec(c + d*x))**(sympy.S(1)/3)*sqrt(I*a*tan(c + d*x) + a))
    F = -3*2**(sympy.S(1)/3)*I*(I*tan(c + d*x) + 1)**(sympy.S(2)/3)*hyper((sympy.S(-1)/6, sympy.S(5)/3), (sympy.S(5)/6,), -I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*(e*sec(c + d*x))**(sympy.S(1)/3)*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_440():
    f = 1/((e*sec(c + d*x))**(sympy.S(4)/3)*sqrt(I*a*tan(c + d*x) + a))
    F = -3*2**(sympy.S(5)/6)*I*(I*tan(c + d*x) + 1)**(sympy.S(1)/6)*sqrt(I*a*tan(c + d*x) + a)*hyper((sympy.S(-2)/3, sympy.S(13)/6), (sympy.S(1)/3,), -I*tan(c + d*x)/2 + sympy.S.Half)/(16*a*d*(e*sec(c + d*x))**(sympy.S(4)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_441():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)/(I*a*tan(e + f*x) + a)**(sympy.S(7)/3)
    F = 5*I*(d*sec(e + f*x))**(sympy.S(2)/3)/(24*f*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a**2*tan(e + f*x) + a**2)) + I*(d*sec(e + f*x))**(sympy.S(2)/3)/(4*f*(I*a*tan(e + f*x) + a)**(sympy.S(7)/3)) - 5*2**(sympy.S(1)/3)*x*(d*sec(e + f*x))**(sympy.S(2)/3)/(144*a**(sympy.S(5)/3)*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) - 5*2**(sympy.S(1)/3)*I*(d*sec(e + f*x))**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (-I*a*tan(e + f*x) + a)**(sympy.S(1)/3))/(48*a**(sympy.S(5)/3)*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) - 5*2**(sympy.S(1)/3)*I*(d*sec(e + f*x))**(sympy.S(2)/3)*log(cos(e + f*x))/(144*a**(sympy.S(5)/3)*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) + 5*2**(sympy.S(1)/3)*sqrt(3)*I*(d*sec(e + f*x))**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(72*a**(sympy.S(5)/3)*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_442():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)/(I*a*tan(e + f*x) + a)**(sympy.S(4)/3)
    F = I*(d*sec(e + f*x))**(sympy.S(2)/3)/(2*f*(I*a*tan(e + f*x) + a)**(sympy.S(4)/3)) - 2**(sympy.S(1)/3)*x*(d*sec(e + f*x))**(sympy.S(2)/3)/(12*a**(sympy.S(2)/3)*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*I*(d*sec(e + f*x))**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (-I*a*tan(e + f*x) + a)**(sympy.S(1)/3))/(4*a**(sympy.S(2)/3)*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*I*(d*sec(e + f*x))**(sympy.S(2)/3)*log(cos(e + f*x))/(12*a**(sympy.S(2)/3)*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*I*(d*sec(e + f*x))**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(6*a**(sympy.S(2)/3)*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_443():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)/(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)
    F = -2**(sympy.S(1)/3)*a**(sympy.S(1)/3)*x*(d*sec(e + f*x))**(sympy.S(2)/3)/(4*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) - 3*2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*(d*sec(e + f*x))**(sympy.S(2)/3)*log(2**(sympy.S(1)/3)*a**(sympy.S(1)/3) - (-I*a*tan(e + f*x) + a)**(sympy.S(1)/3))/(4*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) - 2**(sympy.S(1)/3)*I*a**(sympy.S(1)/3)*(d*sec(e + f*x))**(sympy.S(2)/3)*log(cos(e + f*x))/(4*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) + 2**(sympy.S(1)/3)*sqrt(3)*I*a**(sympy.S(1)/3)*(d*sec(e + f*x))**(sympy.S(2)/3)*atan(sqrt(3)*(a**(sympy.S(1)/3) + 2**(sympy.S(2)/3)*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3))/(3*a**(sympy.S(1)/3)))/(2*f*(-I*a*tan(e + f*x) + a)**(sympy.S(1)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_444():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(2)/3)
    F = 3*I*a*(d*sec(e + f*x))**(sympy.S(2)/3)/(f*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_445():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/3)
    F = 9*I*a**2*(d*sec(e + f*x))**(sympy.S(2)/3)/(2*f*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) + 3*I*a*(d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(2)/3)/(4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_446():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(8)/3)
    F = 54*I*a**3*(d*sec(e + f*x))**(sympy.S(2)/3)/(7*f*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) + 9*I*a**2*(d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(2)/3)/(7*f) + 3*I*a*(d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/3)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_447():
    f = (d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(11)/3)
    F = 486*I*a**4*(d*sec(e + f*x))**(sympy.S(2)/3)/(35*f*(I*a*tan(e + f*x) + a)**(sympy.S(1)/3)) + 81*I*a**3*(d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(2)/3)/(35*f) + 27*I*a**2*(d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(5)/3)/(35*f) + 3*I*a*(d*sec(e + f*x))**(sympy.S(2)/3)*(I*a*tan(e + f*x) + a)**(sympy.S(8)/3)/(10*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_448():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)**5
    F = 2**(m/2 + 5)*I*a**5*(e*sec(c + d*x))**m*hyper((m/2, -m/2 - 4), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*(I*tan(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_449():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)**3
    F = 2**(m/2 + 3)*I*a**3*(e*sec(c + d*x))**m*hyper((m/2, -m/2 - 2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*(I*tan(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_450():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)**2
    F = 2**(m/2 + 2)*I*a**2*(e*sec(c + d*x))**m*hyper((m/2, -m/2 - 1), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*(I*tan(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_451():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)
    F = 2**(m/2 + 1)*I*a*(e*sec(c + d*x))**m*hyper((-m/2, m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*(I*tan(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_452():
    f = (e*sec(c + d*x))**m/(I*a*tan(c + d*x) + a)
    F = 2**(m/2 - 1)*I*(e*sec(c + d*x))**m*hyper((m/2, 2 - m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(a*d*m*(I*tan(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_453():
    f = (e*sec(c + d*x))**m/(I*a*tan(c + d*x) + a)**2
    F = 2**(m/2 - 2)*I*(e*sec(c + d*x))**m*hyper((m/2, 3 - m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(a**2*d*m*(I*tan(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_454():
    f = (e*sec(c + d*x))**m/(I*a*tan(c + d*x) + a)**3
    F = 2**(m/2 - 3)*I*(e*sec(c + d*x))**m*hyper((m/2, 4 - m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(a**3*d*m*(I*tan(c + d*x) + 1)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_455():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)**(sympy.S(7)/2)
    F = 2**(m/2 + sympy.S(7)/2)*I*a**3*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(-m/2 + sympy.S(-1)/2)*sqrt(I*a*tan(c + d*x) + a)*hyper((m/2, -m/2 + sympy.S(-5)/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_456():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2**(m/2 + sympy.S(5)/2)*I*a**2*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(-m/2 + sympy.S(-1)/2)*sqrt(I*a*tan(c + d*x) + a)*hyper((m/2, -m/2 + sympy.S(-3)/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_457():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2**(m/2 + sympy.S(3)/2)*I*a*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(-m/2 + sympy.S(-1)/2)*sqrt(I*a*tan(c + d*x) + a)*hyper((m/2, -m/2 + sympy.S(-1)/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_458():
    f = (e*sec(c + d*x))**m*sqrt(I*a*tan(c + d*x) + a)
    F = 2**(m/2 + sympy.S.Half)*I*a*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(sympy.S.Half - m/2)*hyper((m/2, sympy.S.Half - m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_459():
    f = (e*sec(c + d*x))**m/sqrt(I*a*tan(c + d*x) + a)
    F = 2**(m/2 + sympy.S(-1)/2)*I*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(sympy.S.Half - m/2)*hyper((m/2, sympy.S(3)/2 - m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_460():
    f = (e*sec(c + d*x))**m/(I*a*tan(c + d*x) + a)**(sympy.S(3)/2)
    F = 2**(m/2 + sympy.S(-3)/2)*I*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(sympy.S.Half - m/2)*hyper((m/2, sympy.S(5)/2 - m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(a*d*m*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_461():
    f = (e*sec(c + d*x))**m/(I*a*tan(c + d*x) + a)**(sympy.S(5)/2)
    F = 2**(m/2 + sympy.S(-5)/2)*I*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(sympy.S.Half - m/2)*hyper((m/2, sympy.S(7)/2 - m/2), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(a**2*d*m*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_462():
    f = (e*sec(c + d*x))**m*(I*a*tan(c + d*x) + a)**n
    F = 2**(m/2 + n)*I*(e*sec(c + d*x))**m*(I*tan(c + d*x) + 1)**(-m/2 - n)*(I*a*tan(c + d*x) + a)**n*hyper((m/2, -m/2 - n + 1), (m/2 + 1,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_463():
    f = (I*a*tan(c + d*x) + a)**n*sec(c + d*x)**6
    F = -4*I*(I*a*tan(c + d*x) + a)**(n + 3)/(a**3*d*(n + 3)) + 4*I*(I*a*tan(c + d*x) + a)**(n + 4)/(a**4*d*(n + 4)) - I*(I*a*tan(c + d*x) + a)**(n + 5)/(a**5*d*(n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_464():
    f = (I*a*tan(c + d*x) + a)**n*sec(c + d*x)**4
    F = -2*I*(I*a*tan(c + d*x) + a)**(n + 2)/(a**2*d*(n + 2)) + I*(I*a*tan(c + d*x) + a)**(n + 3)/(a**3*d*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_465():
    f = (I*a*tan(c + d*x) + a)**n*sec(c + d*x)**2
    F = -I*(I*a*tan(c + d*x) + a)**(n + 1)/(a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_466():
    f = (I*a*tan(c + d*x) + a)**n*cos(c + d*x)**2
    F = I*a*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((2, n - 1), (n,), I*tan(c + d*x)/2 + sympy.S.Half)/(4*d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_467():
    f = (I*a*tan(c + d*x) + a)**n*cos(c + d*x)**4
    F = I*a**2*(I*a*tan(c + d*x) + a)**(n - 2)*hyper((3, n - 2), (n - 1,), I*tan(c + d*x)/2 + sympy.S.Half)/(8*d*(2 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_468():
    f = (I*a*tan(c + d*x) + a)**n*cos(c + d*x)**6
    F = I*a**3*(I*a*tan(c + d*x) + a)**(n - 3)*hyper((4, n - 3), (n - 2,), I*tan(c + d*x)/2 + sympy.S.Half)/(16*d*(3 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_469():
    f = (I*a*tan(c + d*x) + a)**n*sec(c + d*x)**5
    F = 2**(n + sympy.S(5)/2)*I*a**2*(I*tan(c + d*x) + 1)**(-n + sympy.S(-1)/2)*(I*a*tan(c + d*x) + a)**(n - 2)*hyper((sympy.S(5)/2, -n + sympy.S(-3)/2), (sympy.S(7)/2,), -I*tan(c + d*x)/2 + sympy.S.Half)*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_470():
    f = (I*a*tan(c + d*x) + a)**n*sec(c + d*x)**3
    F = 2**(n + sympy.S(3)/2)*I*a*(I*tan(c + d*x) + 1)**(-n + sympy.S(-1)/2)*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((sympy.S(3)/2, -n + sympy.S(-1)/2), (sympy.S(5)/2,), -I*tan(c + d*x)/2 + sympy.S.Half)*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_471():
    f = (I*a*tan(c + d*x) + a)**n*sec(c + d*x)
    F = 2**(n + sympy.S.Half)*I*a*(I*tan(c + d*x) + 1)**(sympy.S.Half - n)*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((sympy.S.Half, sympy.S.Half - n), (sympy.S(3)/2,), -I*tan(c + d*x)/2 + sympy.S.Half)*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_472():
    f = (I*a*tan(c + d*x) + a)**n*cos(c + d*x)
    F = -2**(n + sympy.S(-1)/2)*I*(I*tan(c + d*x) + 1)**(sympy.S.Half - n)*(I*a*tan(c + d*x) + a)**n*cos(c + d*x)*hyper((sympy.S(-1)/2, sympy.S(3)/2 - n), (sympy.S.Half,), -I*tan(c + d*x)/2 + sympy.S.Half)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_473():
    f = (I*a*tan(c + d*x) + a)**n*cos(c + d*x)**3
    F = -2**(n + sympy.S(-3)/2)*I*(I*tan(c + d*x) + 1)**(sympy.S.Half - n)*(I*a*tan(c + d*x) + a)**(n + 1)*cos(c + d*x)**3*hyper((sympy.S(-3)/2, sympy.S(5)/2 - n), (sympy.S(-1)/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_474():
    f = (I*a*tan(c + d*x) + a)**n*cos(c + d*x)**5
    F = -2**(n + sympy.S(-5)/2)*I*(I*tan(c + d*x) + 1)**(sympy.S.Half - n)*(I*a*tan(c + d*x) + a)**(n + 2)*cos(c + d*x)**5*hyper((sympy.S(-5)/2, sympy.S(7)/2 - n), (sympy.S(-3)/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(5*a**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_475():
    f = (e*sec(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**n
    F = 2**(n + sympy.S(9)/4)*I*a*(e*sec(c + d*x))**(sympy.S(5)/2)*(I*tan(c + d*x) + 1)**(-n + sympy.S(-1)/4)*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((sympy.S(5)/4, -n + sympy.S(-1)/4), (sympy.S(9)/4,), -I*tan(c + d*x)/2 + sympy.S.Half)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_476():
    f = (e*sec(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**n
    F = 2**(n + sympy.S(7)/4)*I*a*(e*sec(c + d*x))**(sympy.S(3)/2)*(I*tan(c + d*x) + 1)**(sympy.S(1)/4 - n)*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((sympy.S(3)/4, sympy.S(1)/4 - n), (sympy.S(7)/4,), -I*tan(c + d*x)/2 + sympy.S.Half)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_477():
    f = sqrt(e*sec(c + d*x))*(I*a*tan(c + d*x) + a)**n
    F = 2**(n + sympy.S(5)/4)*I*a*sqrt(e*sec(c + d*x))*(I*tan(c + d*x) + 1)**(sympy.S(3)/4 - n)*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((sympy.S(1)/4, sympy.S(3)/4 - n), (sympy.S(5)/4,), -I*tan(c + d*x)/2 + sympy.S.Half)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_478():
    f = (I*a*tan(c + d*x) + a)**n/sqrt(e*sec(c + d*x))
    F = -2**(n + sympy.S(3)/4)*I*(I*tan(c + d*x) + 1)**(sympy.S(1)/4 - n)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S(-1)/4, sympy.S(5)/4 - n), (sympy.S(3)/4,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*sqrt(e*sec(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_479():
    f = (I*a*tan(c + d*x) + a)**n/(e*sec(c + d*x))**(sympy.S(3)/2)
    F = -2**(n + sympy.S(1)/4)*I*(I*tan(c + d*x) + 1)**(sympy.S(3)/4 - n)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S(-3)/4, sympy.S(7)/4 - n), (sympy.S(1)/4,), -I*tan(c + d*x)/2 + sympy.S.Half)/(3*d*(e*sec(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_480():
    f = (I*a*tan(c + d*x) + a)**n/(e*sec(c + d*x))**(sympy.S(5)/2)
    F = -2**(n + sympy.S(-1)/4)*I*(I*tan(c + d*x) + 1)**(sympy.S(1)/4 - n)*(I*a*tan(c + d*x) + a)**(n + 1)*hyper((sympy.S(-5)/4, sympy.S(9)/4 - n), (sympy.S(-1)/4,), -I*tan(c + d*x)/2 + sympy.S.Half)/(5*a*d*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_481():
    f = (e*sec(c + d*x))**(-n - 4)*(I*a*tan(c + d*x) + a)**n
    F = I*(e*sec(c + d*x))**(-n - 4)*(I*a*tan(c + d*x) + a)**n/(d*(4 - n)) + 4*I*(e*sec(c + d*x))**(-n - 4)*(I*a*tan(c + d*x) + a)**(n + 1)/(a*d*(n**2 - 6*n + 8)) - 12*I*(e*sec(c + d*x))**(-n - 4)*(I*a*tan(c + d*x) + a)**(n + 2)/(a**2*d*n*(2 - n)*(4 - n)) + 24*I*(e*sec(c + d*x))**(-n - 4)*(I*a*tan(c + d*x) + a)**(n + 3)/(a**3*d*n*(4 - n)*(4 - n**2)) - 24*I*(e*sec(c + d*x))**(-n - 4)*(I*a*tan(c + d*x) + a)**(n + 4)/(a**4*d*n*(n**4 - 20*n**2 + 64))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_482():
    f = (e*sec(c + d*x))**(-n - 3)*(I*a*tan(c + d*x) + a)**n
    F = I*(e*sec(c + d*x))**(-n - 3)*(I*a*tan(c + d*x) + a)**n/(d*(3 - n)) + 3*I*(e*sec(c + d*x))**(-n - 3)*(I*a*tan(c + d*x) + a)**(n + 1)/(a*d*(n**2 - 4*n + 3)) - 6*I*(e*sec(c + d*x))**(-n - 3)*(I*a*tan(c + d*x) + a)**(n + 2)/(a**2*d*(1 - n**2)*(3 - n)) + 6*I*(e*sec(c + d*x))**(-n - 3)*(I*a*tan(c + d*x) + a)**(n + 3)/(a**3*d*(n**4 - 10*n**2 + 9))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_483():
    f = (e*sec(c + d*x))**(-n - 2)*(I*a*tan(c + d*x) + a)**n
    F = I*(e*sec(c + d*x))**(-n - 2)*(I*a*tan(c + d*x) + a)**n/(d*(2 - n)) - 2*I*(e*sec(c + d*x))**(-n - 2)*(I*a*tan(c + d*x) + a)**(n + 1)/(a*d*n*(2 - n)) + 2*I*(e*sec(c + d*x))**(-n - 2)*(I*a*tan(c + d*x) + a)**(n + 2)/(a**2*d*n*(4 - n**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_484():
    f = (e*sec(c + d*x))**(-n - 1)*(I*a*tan(c + d*x) + a)**n
    F = I*(e*sec(c + d*x))**(-n - 1)*(I*a*tan(c + d*x) + a)**n/(d*(1 - n)) - I*(e*sec(c + d*x))**(-n - 1)*(I*a*tan(c + d*x) + a)**(n + 1)/(a*d*(1 - n**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_485():
    f = (I*a*tan(c + d*x) + a)**n/(e*sec(c + d*x))**n
    F = -I*(I*a*tan(c + d*x) + a)**n/(d*n*(e*sec(c + d*x))**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_486():
    f = (e*sec(c + d*x))**(1 - n)*(I*a*tan(c + d*x) + a)**n
    F = 2**(n/2 + sympy.S.Half)*I*(e*sec(c + d*x))**(1 - n)*(I*tan(c + d*x) + 1)**(-n/2 + sympy.S(-1)/2)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S.Half - n/2, sympy.S.Half - n/2), (sympy.S(3)/2 - n/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_487():
    f = (e*sec(c + d*x))**(2 - n)*(I*a*tan(c + d*x) + a)**n
    F = 2**(n/2 + 1)*I*a*(e*sec(c + d*x))**(2 - n)*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((-n/2, 1 - n/2), (2 - n/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*(2 - n)*(I*tan(c + d*x) + 1)**(n/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_488():
    f = (e*sec(c + d*x))**(3 - n)*(I*a*tan(c + d*x) + a)**n
    F = 2**(n/2 + sympy.S(3)/2)*I*a*(e*sec(c + d*x))**(3 - n)*(I*tan(c + d*x) + 1)**(-n/2 + sympy.S(-1)/2)*(I*a*tan(c + d*x) + a)**(n - 1)*hyper((sympy.S(3)/2 - n/2, -n/2 + sympy.S(-1)/2), (sympy.S(5)/2 - n/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*(3 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_489():
    f = (e*sec(c + d*x))**(6 - 2*n)*(I*a*tan(c + d*x) + a)**n
    F = 8*I*a**3*(e*sec(c + d*x))**(6 - 2*n)*(I*a*tan(c + d*x) + a)**(n - 3)/(d*(5 - n)*(n**2 - 7*n + 12)) + 4*I*a**2*(e*sec(c + d*x))**(6 - 2*n)*(I*a*tan(c + d*x) + a)**(n - 2)/(d*(n**2 - 9*n + 20)) + I*a*(e*sec(c + d*x))**(6 - 2*n)*(I*a*tan(c + d*x) + a)**(n - 1)/(d*(5 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_490():
    f = (e*sec(c + d*x))**(5 - 2*n)*(I*a*tan(c + d*x) + a)**n
    F = -2**(sympy.S(5)/2 - n)*I*(e*sec(c + d*x))**(5 - 2*n)*(-I*tan(c + d*x) + 1)**(n + sympy.S(-5)/2)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S(5)/2, n + sympy.S(-3)/2), (sympy.S(7)/2,), I*tan(c + d*x)/2 + sympy.S.Half)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_491():
    f = (e*sec(c + d*x))**(4 - 2*n)*(I*a*tan(c + d*x) + a)**n
    F = 2*I*a**2*(e*sec(c + d*x))**(4 - 2*n)*(I*a*tan(c + d*x) + a)**(n - 2)/(d*(n**2 - 5*n + 6)) + I*a*(e*sec(c + d*x))**(4 - 2*n)*(I*a*tan(c + d*x) + a)**(n - 1)/(d*(3 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_492():
    f = (e*sec(c + d*x))**(3 - 2*n)*(I*a*tan(c + d*x) + a)**n
    F = -2**(sympy.S(3)/2 - n)*I*(e*sec(c + d*x))**(3 - 2*n)*(-I*tan(c + d*x) + 1)**(n + sympy.S(-3)/2)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S(3)/2, n + sympy.S(-1)/2), (sympy.S(5)/2,), I*tan(c + d*x)/2 + sympy.S.Half)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_493():
    f = (e*sec(c + d*x))**(2 - 2*n)*(I*a*tan(c + d*x) + a)**n
    F = I*a*(e*sec(c + d*x))**(2 - 2*n)*(I*a*tan(c + d*x) + a)**(n - 1)/(d*(1 - n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_494():
    f = (e*sec(c + d*x))**(1 - 2*n)*(I*a*tan(c + d*x) + a)**n
    F = -2**(sympy.S.Half - n)*I*(e*sec(c + d*x))**(1 - 2*n)*(-I*tan(c + d*x) + 1)**(n + sympy.S(-1)/2)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S.Half, n + sympy.S.Half), (sympy.S(3)/2,), I*tan(c + d*x)/2 + sympy.S.Half)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_495():
    f = (I*a*tan(c + d*x) + a)**n/(e*sec(c + d*x))**(2*n)
    F = -I*(I*a*tan(c + d*x) + a)**n*hyper((1, -n), (1 - n,), -I*tan(c + d*x)/2 + sympy.S.Half)/(2*d*n*(e*sec(c + d*x))**(2*n))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_496():
    f = (e*sec(c + d*x))**(-2*n - 1)*(I*a*tan(c + d*x) + a)**n
    F = 2**(-n + sympy.S(-1)/2)*I*(e*sec(c + d*x))**(-2*n - 1)*(-I*tan(c + d*x) + 1)**(n + sympy.S.Half)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S(-1)/2, n + sympy.S(3)/2), (sympy.S.Half,), I*tan(c + d*x)/2 + sympy.S.Half)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_497():
    f = (e*sec(c + d*x))**(-2*n - 2)*(I*a*tan(c + d*x) + a)**n
    F = -I*(e*sec(c + d*x))**(-2*n - 2)*(I*a*tan(c + d*x) + a)**(n + 1)*hyper((2, -n - 1), (-n,), -I*tan(c + d*x)/2 + sympy.S.Half)/(4*a*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_498():
    f = (e*sec(c + d*x))**(-2*n - 3)*(I*a*tan(c + d*x) + a)**n
    F = 2**(-n + sympy.S(-3)/2)*I*(e*sec(c + d*x))**(-2*n - 3)*(-I*tan(c + d*x) + 1)**(n + sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**n*hyper((sympy.S(-3)/2, n + sympy.S(5)/2), (sympy.S(-1)/2,), I*tan(c + d*x)/2 + sympy.S.Half)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_499():
    f = (d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(-n - 2)
    F = I*(d*sec(e + f*x))**(2*n)*hyper((3, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(8*a**2*f*n*(I*a*tan(e + f*x) + a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_500():
    f = (d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(-n - 1)
    F = I*(d*sec(e + f*x))**(2*n)*hyper((2, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(4*a*f*n*(I*a*tan(e + f*x) + a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_501():
    f = (d*sec(e + f*x))**(2*n)/(I*a*tan(e + f*x) + a)**n
    F = I*(d*sec(e + f*x))**(2*n)*hyper((1, n), (n + 1,), -I*tan(e + f*x)/2 + sympy.S.Half)/(2*f*n*(I*a*tan(e + f*x) + a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_502():
    f = (d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(1 - n)
    F = I*a*(d*sec(e + f*x))**(2*n)/(f*n*(I*a*tan(e + f*x) + a)**n)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_503():
    f = (d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(2 - n)
    F = 2*I*a**2*(d*sec(e + f*x))**(2*n)/(f*n*(n + 1)*(I*a*tan(e + f*x) + a)**n) + I*a*(d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(1 - n)/(f*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_504():
    f = (d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(3 - n)
    F = 8*I*a**3*(d*sec(e + f*x))**(2*n)/(f*n*(I*a*tan(e + f*x) + a)**n*(n**2 + 3*n + 2)) + 4*I*a**2*(d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(1 - n)/(f*(n**2 + 3*n + 2)) + I*a*(d*sec(e + f*x))**(2*n)*(I*a*tan(e + f*x) + a)**(2 - n)/(f*(n + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_505():
    f = (a + b*tan(c + d*x))*sec(c + d*x)**6
    F = a*tan(c + d*x)**5/(5*d) + 2*a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*sec(c + d*x)**6/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_506():
    f = (a + b*tan(c + d*x))*sec(c + d*x)**5
    F = a*tan(c + d*x)*sec(c + d*x)**3/(4*d) + 3*a*tan(c + d*x)*sec(c + d*x)/(8*d) + 3*a*atanh(sin(c + d*x))/(8*d) + b*sec(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_507():
    f = (a + b*tan(c + d*x))*sec(c + d*x)**4
    F = a*tan(c + d*x)**3/(3*d) + a*tan(c + d*x)/d + b*sec(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_508():
    f = (a + b*tan(c + d*x))*sec(c + d*x)**3
    F = a*tan(c + d*x)*sec(c + d*x)/(2*d) + a*atanh(sin(c + d*x))/(2*d) + b*sec(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_509():
    f = (a + b*tan(c + d*x))*sec(c + d*x)**2
    F = a*tan(c + d*x)/d + b*sec(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_510():
    f = (a + b*tan(c + d*x))*sec(c + d*x)
    F = a*atanh(sin(c + d*x))/d + b*sec(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_511():
    f = (a + b*tan(c + d*x))*cos(c + d*x)
    F = a*sin(c + d*x)/d - b*cos(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_512():
    f = (a + b*tan(c + d*x))*cos(c + d*x)**2
    F = a*x/2 + a*sin(c + d*x)*cos(c + d*x)/(2*d) - b*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_513():
    f = (a + b*tan(c + d*x))*cos(c + d*x)**3
    F = -a*sin(c + d*x)**3/(3*d) + a*sin(c + d*x)/d - b*cos(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_514():
    f = (a + b*tan(c + d*x))*cos(c + d*x)**4
    F = 3*a*x/8 + a*sin(c + d*x)*cos(c + d*x)**3/(4*d) + 3*a*sin(c + d*x)*cos(c + d*x)/(8*d) - b*cos(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_515():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)**8
    F = a**2*tan(c + d*x)/d + a*b*sec(c + d*x)**8/(4*d) + b**2*tan(c + d*x)**9/(9*d) + (a**2 + 3*b**2)*tan(c + d*x)**7/(7*d) + (3*a**2 + b**2)*tan(c + d*x)**3/(3*d) + (3*a**2 + 3*b**2)*tan(c + d*x)**5/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_516():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)**6
    F = a**2*tan(c + d*x)/d + a*b*sec(c + d*x)**6/(3*d) + b**2*tan(c + d*x)**7/(7*d) + (a**2 + 2*b**2)*tan(c + d*x)**5/(5*d) + (2*a**2 + b**2)*tan(c + d*x)**3/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_517():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)**4
    F = -a*(a + b*tan(c + d*x))**4/(2*b**3*d) + (a + b*tan(c + d*x))**5/(5*b**3*d) + (a + b*tan(c + d*x))**3*(a**2 + b**2)/(3*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_518():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)**2
    F = (a + b*tan(c + d*x))**3/(3*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_519():
    f = (a + b*tan(c + d*x))**2*cos(c + d*x)**2
    F = x*(a**2/2 + b**2/2) - (a + b*tan(c + d*x))*(-a*tan(c + d*x) + b)*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_520():
    f = (a + b*tan(c + d*x))**2*cos(c + d*x)**4
    F = x*(3*a**2/8 + b**2/8) - (a + b*tan(c + d*x))*(-a*tan(c + d*x) + b)*cos(c + d*x)**4/(4*d) - (2*a*b - (3*a**2 + b**2)*tan(c + d*x))*cos(c + d*x)**2/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_521():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)**7
    F = 9*a*b*sec(c + d*x)**7/(56*d) + b*(a + b*tan(c + d*x))*sec(c + d*x)**7/(8*d) + (8*a**2 - b**2)*tan(c + d*x)*sec(c + d*x)**5/(48*d) + (40*a**2 - 5*b**2)*tan(c + d*x)*sec(c + d*x)**3/(192*d) + (40*a**2 - 5*b**2)*tan(c + d*x)*sec(c + d*x)/(128*d) + (40*a**2 - 5*b**2)*atanh(sin(c + d*x))/(128*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_522():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)**5
    F = 7*a*b*sec(c + d*x)**5/(30*d) + b*(a + b*tan(c + d*x))*sec(c + d*x)**5/(6*d) + (6*a**2 - b**2)*tan(c + d*x)*sec(c + d*x)**3/(24*d) + (6*a**2 - b**2)*tan(c + d*x)*sec(c + d*x)/(16*d) + (6*a**2 - b**2)*atanh(sin(c + d*x))/(16*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_523():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)**3
    F = 5*a*b*sec(c + d*x)**3/(12*d) + b*(a + b*tan(c + d*x))*sec(c + d*x)**3/(4*d) + (4*a**2 - b**2)*tan(c + d*x)*sec(c + d*x)/(8*d) + (4*a**2 - b**2)*atanh(sin(c + d*x))/(8*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_524():
    f = (a + b*tan(c + d*x))**2*sec(c + d*x)
    F = 3*a*b*sec(c + d*x)/(2*d) + b*(a + b*tan(c + d*x))*sec(c + d*x)/(2*d) + (2*a**2 - b**2)*atanh(sin(c + d*x))/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_525():
    f = (a + b*tan(c + d*x))**2*cos(c + d*x)
    F = -2*a*b*cos(c + d*x)/d + b**2*atanh(sin(c + d*x))/d + (a**2 - b**2)*sin(c + d*x)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_526():
    f = (a + b*tan(c + d*x))**2*cos(c + d*x)**3
    F = -a*b*cos(c + d*x)**3/(6*d) - b*(a + b*tan(c + d*x))*cos(c + d*x)**3/(2*d) - (2*a**2 + b**2)*sin(c + d*x)**3/(6*d) + (2*a**2 + b**2)*sin(c + d*x)/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_527():
    f = (a + b*tan(c + d*x))**2*cos(c + d*x)**5
    F = -3*a*b*cos(c + d*x)**5/(20*d) - b*(a + b*tan(c + d*x))*cos(c + d*x)**5/(4*d) + (4*a**2 + b**2)*sin(c + d*x)**5/(20*d) - (4*a**2 + b**2)*sin(c + d*x)**3/(6*d) + (4*a**2 + b**2)*sin(c + d*x)/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_528():
    f = (a + b*tan(c + d*x))**2*cos(c + d*x)**7
    F = -5*a*b*cos(c + d*x)**7/(42*d) - b*(a + b*tan(c + d*x))*cos(c + d*x)**7/(6*d) - (6*a**2 + b**2)*sin(c + d*x)**7/(42*d) + (6*a**2 + b**2)*sin(c + d*x)**5/(10*d) - (6*a**2 + b**2)*sin(c + d*x)**3/(6*d) + (6*a**2 + b**2)*sin(c + d*x)/(6*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_529():
    f = (a + b*tan(c + d*x))**3*sec(c + d*x)**8
    F = a**3*tan(c + d*x)/d + 3*a**2*b*sec(c + d*x)**8/(8*d) + a*b**2*tan(c + d*x)**9/(3*d) + a*(a**2 + b**2)*tan(c + d*x)**3/d + 3*a*(a**2 + 3*b**2)*tan(c + d*x)**5/(5*d) + a*(a**2 + 9*b**2)*tan(c + d*x)**7/(7*d) + b**3*tan(c + d*x)**10/(10*d) + 3*b**3*tan(c + d*x)**8/(8*d) + b**3*tan(c + d*x)**6/(2*d) + b**3*tan(c + d*x)**4/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_530():
    f = (a + b*tan(c + d*x))**3*sec(c + d*x)**6
    F = -4*a*(a + b*tan(c + d*x))**7/(7*b**5*d) - 4*a*(a + b*tan(c + d*x))**5*(a**2 + b**2)/(5*b**5*d) + (a + b*tan(c + d*x))**8/(8*b**5*d) + (a + b*tan(c + d*x))**6*(3*a**2 + b**2)/(3*b**5*d) + (a + b*tan(c + d*x))**4*(a**2 + b**2)**2/(4*b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_531():
    f = (a + b*tan(c + d*x))**3*sec(c + d*x)**4
    F = -2*a*(a + b*tan(c + d*x))**5/(5*b**3*d) + (a + b*tan(c + d*x))**6/(6*b**3*d) + (a + b*tan(c + d*x))**4*(a**2 + b**2)/(4*b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_532():
    f = (a + b*tan(c + d*x))**3*sec(c + d*x)**2
    F = (a + b*tan(c + d*x))**4/(4*b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_533():
    f = (a + b*tan(c + d*x))**3*cos(c + d*x)**2
    F = -a*b**2*tan(c + d*x)/(2*d) + a*x*(a**2 + 3*b**2)/2 - b**3*log(cos(c + d*x))/d - (a + b*tan(c + d*x))**2*(-a*tan(c + d*x) + b)*cos(c + d*x)**2/(2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_534():
    f = (a + b*tan(c + d*x))**3*cos(c + d*x)**4
    F = 3*a*x*(a**2 + b**2)/8 - 3*a*(a + b*tan(c + d*x))*(-a*tan(c + d*x) + b)*cos(c + d*x)**2/(8*d) + (a + b*tan(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**3/(4*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_535():
    f = (a + b*tan(c + d*x))**3*cos(c + d*x)**3
    F = -(a + b*tan(c + d*x))**2*(-a*tan(c + d*x) + b)*cos(c + d*x)**3/(3*d) - (2*a**2 + 2*b**2)*(-a*tan(c + d*x) + b)*cos(c + d*x)/(3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_536():
    f = (a + b*tan(c + d*x))**3*cos(c + d*x)**5
    F = (a + b*tan(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**4/(5*d) - (a + b*tan(c + d*x))**2*(-4*a*tan(c + d*x) + b)*cos(c + d*x)**3/(15*d) - (8*a**2 + 2*b**2)*(-a*tan(c + d*x) + b)*cos(c + d*x)/(15*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_537():
    f = (a + b*tan(c + d*x))**3*cos(c + d*x)**7
    F = 8*a*(2*a**2 + b**2)*sin(c + d*x)/(35*d) + (a + b*tan(c + d*x))**3*sin(c + d*x)*cos(c + d*x)**6/(7*d) - 3*(a + b*tan(c + d*x))**2*(-2*a*tan(c + d*x) + b)*cos(c + d*x)**5/(35*d) - 2*(-a*(4*a**2 - b**2)*tan(c + d*x) + b*(6*a**2 + b**2))*cos(c + d*x)**3/(35*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_538():
    f = sec(c + d*x)**6/(a + b*tan(c + d*x))
    F = -a*tan(c + d*x)**3/(3*b**2*d) - a*(a**2 + 2*b**2)*tan(c + d*x)/(b**4*d) + tan(c + d*x)**4/(4*b*d) + (a**2 + 2*b**2)*tan(c + d*x)**2/(2*b**3*d) + (a**2 + b**2)**2*log(a + b*tan(c + d*x))/(b**5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_539():
    f = sec(c + d*x)**4/(a + b*tan(c + d*x))
    F = -a*tan(c + d*x)/(b**2*d) + tan(c + d*x)**2/(2*b*d) + (a**2 + b**2)*log(a + b*tan(c + d*x))/(b**3*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_540():
    f = sec(c + d*x)**2/(a + b*tan(c + d*x))
    F = log(a + b*tan(c + d*x))/(b*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_541():
    f = cos(c + d*x)**2/(a + b*tan(c + d*x))
    F = a*x*(a**2 + 3*b**2)/(2*(a**2 + b**2)**2) + b**3*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**2) + (a*tan(c + d*x) + b)*cos(c + d*x)**2/(d*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_542():
    f = cos(c + d*x)**4/(a + b*tan(c + d*x))
    F = a*x*(3*a**4 + 10*a**2*b**2 + 15*b**4)/(8*(a**2 + b**2)**3) + b**5*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + (a*tan(c + d*x) + b)*cos(c + d*x)**4/(d*(4*a**2 + 4*b**2)) + (a*(3*a**2 + 7*b**2)*tan(c + d*x) + 4*b**3)*cos(c + d*x)**2/(8*d*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_543():
    f = sec(c + d*x)**3/(a + b*tan(c + d*x))
    F = -a*atanh(sin(c + d*x))/(b**2*d) + sec(c + d*x)/(b*d) - sqrt(a**2 + b**2)*atanh((-a*tan(c + d*x) + b)*cos(c + d*x)/sqrt(a**2 + b**2))/(b**2*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_544():
    f = sec(c + d*x)/(a + b*tan(c + d*x))
    F = -atanh((-a*tan(c + d*x) + b)*cos(c + d*x)/sqrt(a**2 + b**2))/(d*sqrt(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_545():
    f = cos(c + d*x)/(a + b*tan(c + d*x))
    F = a*sin(c + d*x)/(d*(a**2 + b**2)) - b**2*atanh((-a*tan(c + d*x) + b)*cos(c + d*x)/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(3)/2)) + b*cos(c + d*x)/(d*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_546():
    f = cos(c + d*x)**3/(a + b*tan(c + d*x))
    F = a*b**2*sin(c + d*x)/(d*(a**2 + b**2)**2) - a*sin(c + d*x)**3/(d*(3*a**2 + 3*b**2)) + a*sin(c + d*x)/(d*(a**2 + b**2)) - b**4*atanh((-a*tan(c + d*x) + b)*cos(c + d*x)/sqrt(a**2 + b**2))/(d*(a**2 + b**2)**(sympy.S(5)/2)) + b**3*cos(c + d*x)/(d*(a**2 + b**2)**2) + b*cos(c + d*x)**3/(d*(3*a**2 + 3*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_547():
    f = sec(c + d*x)**8/(a + b*tan(c + d*x))**2
    F = -a*tan(c + d*x)**4/(2*b**3*d) - a*(2*a**2 + 3*b**2)*tan(c + d*x)**2/(b**5*d) - 6*a*(a**2 + b**2)**2*log(a + b*tan(c + d*x))/(b**7*d) + tan(c + d*x)**5/(5*b**2*d) + (a**2 + b**2)*tan(c + d*x)**3/(b**4*d) + (5*a**4 + 9*a**2*b**2 + 3*b**4)*tan(c + d*x)/(b**6*d) - (a**2 + b**2)**3/(b**7*d*(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_548():
    f = sec(c + d*x)**6/(a + b*tan(c + d*x))**2
    F = -a*tan(c + d*x)**2/(b**3*d) - 4*a*(a**2 + b**2)*log(a + b*tan(c + d*x))/(b**5*d) + tan(c + d*x)**3/(3*b**2*d) + (3*a**2 + 2*b**2)*tan(c + d*x)/(b**4*d) - (a**2 + b**2)**2/(b**5*d*(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_549():
    f = sec(c + d*x)**4/(a + b*tan(c + d*x))**2
    F = -2*a*log(a + b*tan(c + d*x))/(b**3*d) + tan(c + d*x)/(b**2*d) - (a**2 + b**2)/(b**3*d*(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_550():
    f = sec(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = -1/(b*d*(a + b*tan(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_551():
    f = cos(c + d*x)**2/(a + b*tan(c + d*x))**2
    F = 4*a*b**3*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**3) + b*(a**2 - 3*b**2)/(2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + x*(a**4 + 6*a**2*b**2 - 3*b**4)/(2*(a**2 + b**2)**3) + (a*tan(c + d*x) + b)*cos(c + d*x)**2/(d*(a + b*tan(c + d*x))*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_552():
    f = cos(c + d*x)**4/(a + b*tan(c + d*x))**2
    F = 6*a*b**5*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + 3*b*(a**2 - b**2)*(a**2 + 5*b**2)/(8*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + x*(3*a**6 + 15*a**4*b**2 + 45*a**2*b**4 - 15*b**6)/(8*(a**2 + b**2)**4) + (a*tan(c + d*x) + b)*cos(c + d*x)**4/(d*(a + b*tan(c + d*x))*(4*a**2 + 4*b**2)) - (-3*a*(a**2 + 3*b**2)*tan(c + d*x) + b*(a**2 - 5*b**2))*cos(c + d*x)**2/(8*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_553():
    f = sec(c + d*x)**7/(a + b*tan(c + d*x))**2
    F = 5*a*(a**2 + b**2)**(sympy.S(3)/2)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))*sec(c + d*x)/(b**6*d*sqrt(sec(c + d*x)**2)) - sec(c + d*x)**5/(b*d*(a + b*tan(c + d*x))) - 5*(4*a - 3*b*tan(c + d*x))*sec(c + d*x)**3/(12*b**3*d) - 5*(8*a*(a**2 + b**2) - b*(4*a**2 + 3*b**2)*tan(c + d*x))*sec(c + d*x)/(8*b**5*d) + (40*a**4 + 60*a**2*b**2 + 15*b**4)*asinh(tan(c + d*x))*sec(c + d*x)/(8*b**6*d*sqrt(sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_554():
    f = sec(c + d*x)**5/(a + b*tan(c + d*x))**2
    F = 3*a*sqrt(a**2 + b**2)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))*sec(c + d*x)/(b**4*d*sqrt(sec(c + d*x)**2)) - sec(c + d*x)**3/(b*d*(a + b*tan(c + d*x))) - 3*(2*a - b*tan(c + d*x))*sec(c + d*x)/(2*b**3*d) + (6*a**2 + 3*b**2)*asinh(tan(c + d*x))*sec(c + d*x)/(2*b**4*d*sqrt(sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_555():
    f = cos(c + d*x)/(a + b*tan(c + d*x))**2
    F = -3*a*b**2*sqrt(sec(c + d*x)**2)*cos(c + d*x)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))/(d*(a**2 + b**2)**(sympy.S(5)/2)) + b*(a**2 - 2*b**2)*sec(c + d*x)/(d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) + (a*tan(c + d*x) + b)*cos(c + d*x)/(d*(a + b*tan(c + d*x))*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_556():
    f = cos(c + d*x)**3/(a + b*tan(c + d*x))**2
    F = -5*a*b**4*sqrt(sec(c + d*x)**2)*cos(c + d*x)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))/(d*(a**2 + b**2)**(sympy.S(7)/2)) + b*(2*a**4 + 9*a**2*b**2 - 8*b**4)*sec(c + d*x)/(3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + (a*tan(c + d*x) + b)*cos(c + d*x)**3/(d*(a + b*tan(c + d*x))*(3*a**2 + 3*b**2)) - (-a*(2*a**2 + 7*b**2)*tan(c + d*x) + b*(a**2 - 4*b**2))*cos(c + d*x)/(3*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_557():
    f = sec(c + d*x)**8/(a + b*tan(c + d*x))**3
    F = -a*tan(c + d*x)**3/(b**4*d) - a*(10*a**2 + 9*b**2)*tan(c + d*x)/(b**6*d) + 6*a*(a**2 + b**2)**2/(b**7*d*(a + b*tan(c + d*x))) + tan(c + d*x)**4/(4*b**3*d) + (6*a**2 + 3*b**2)*tan(c + d*x)**2/(2*b**5*d) + (3*a**2 + 3*b**2)*(5*a**2 + b**2)*log(a + b*tan(c + d*x))/(b**7*d) - (a**2 + b**2)**3/(2*b**7*d*(a + b*tan(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_558():
    f = sec(c + d*x)**6/(a + b*tan(c + d*x))**3
    F = -3*a*tan(c + d*x)/(b**4*d) + 4*a*(a**2 + b**2)/(b**5*d*(a + b*tan(c + d*x))) + tan(c + d*x)**2/(2*b**3*d) + (6*a**2 + 2*b**2)*log(a + b*tan(c + d*x))/(b**5*d) - (a**2 + b**2)**2/(2*b**5*d*(a + b*tan(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_559():
    f = sec(c + d*x)**4/(a + b*tan(c + d*x))**3
    F = 2*a/(b**3*d*(a + b*tan(c + d*x))) + log(a + b*tan(c + d*x))/(b**3*d) - (a**2 + b**2)/(2*b**3*d*(a + b*tan(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_560():
    f = sec(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = -1/(2*b*d*(a + b*tan(c + d*x))**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_561():
    f = cos(c + d*x)**2/(a + b*tan(c + d*x))**3
    F = a*b*(a**2 - 11*b**2)/(2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) + a*x*(a**4 + 10*a**2*b**2 - 15*b**4)/(2*(a**2 + b**2)**4) + 2*b**3*(5*a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**4) + b*(a**2 - 2*b**2)/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + (a*tan(c + d*x) + b)*cos(c + d*x)**2/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_562():
    f = cos(c + d*x)**4/(a + b*tan(c + d*x))**3
    F = 3*a*b*(a**4 + 6*a**2*b**2 - 27*b**4)/(8*d*(a + b*tan(c + d*x))*(a**2 + b**2)**4) + 3*a*x*(a**6 + 7*a**4*b**2 + 35*a**2*b**4 - 35*b**6)/(8*(a**2 + b**2)**5) + 3*b**5*(7*a**2 - b**2)*log(a*cos(c + d*x) + b*sin(c + d*x))/(d*(a**2 + b**2)**5) + 3*b*(a**4 + 5*a**2*b**2 - 4*b**4)/(8*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**3) + (a*tan(c + d*x) + b)*cos(c + d*x)**4/(d*(a + b*tan(c + d*x))**2*(4*a**2 + 4*b**2)) - (-a*(3*a**2 + 11*b**2)*tan(c + d*x) + 2*b*(a**2 - 3*b**2))*cos(c + d*x)**2/(8*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_563():
    f = sec(c + d*x)**7/(a + b*tan(c + d*x))**3
    F = -5*a*(4*a**2 + 3*b**2)*asinh(tan(c + d*x))*sec(c + d*x)/(2*b**6*d*sqrt(sec(c + d*x)**2)) - sec(c + d*x)**5/(2*b*d*(a + b*tan(c + d*x))**2) + 5*(4*a + b*tan(c + d*x))*sec(c + d*x)**3/(6*b**3*d*(a + b*tan(c + d*x))) + 5*(4*a**2 - 2*a*b*tan(c + d*x) + b**2)*sec(c + d*x)/(2*b**5*d) - 5*sqrt(a**2 + b**2)*(4*a**2 + b**2)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))*sec(c + d*x)/(2*b**6*d*sqrt(sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_564():
    f = sec(c + d*x)/(a + b*tan(c + d*x))**3
    F = -3*a*b*sec(c + d*x)/(2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**2) - b*sec(c + d*x)/(d*(a + b*tan(c + d*x))**2*(2*a**2 + 2*b**2)) - (2*a**2 - b**2)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))*sec(c + d*x)/(2*d*(a**2 + b**2)**(sympy.S(5)/2)*sqrt(sec(c + d*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_565():
    f = cos(c + d*x)/(a + b*tan(c + d*x))**3
    F = a*b*(2*a**2 - 13*b**2)*sec(c + d*x)/(2*d*(a + b*tan(c + d*x))*(a**2 + b**2)**3) - 3*b**2*(4*a**2 - b**2)*sqrt(sec(c + d*x)**2)*cos(c + d*x)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))/(2*d*(a**2 + b**2)**(sympy.S(7)/2)) + b*(2*a**2 - 3*b**2)*sec(c + d*x)/(2*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2) + (a*tan(c + d*x) + b)*cos(c + d*x)/(d*(a + b*tan(c + d*x))**2*(a**2 + b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_566():
    f = cos(c + d*x)**3/(a + b*tan(c + d*x))**3
    F = a*b*(4*a**4 + 28*a**2*b**2 - 81*b**4)*sec(c + d*x)/(6*d*(a + b*tan(c + d*x))*(a**2 + b**2)**4) - 5*b**4*(6*a**2 - b**2)*sqrt(sec(c + d*x)**2)*cos(c + d*x)*atanh((-a*tan(c + d*x) + b)/(sqrt(a**2 + b**2)*sqrt(sec(c + d*x)**2)))/(2*d*(a**2 + b**2)**(sympy.S(9)/2)) + b*(4*a**4 + 24*a**2*b**2 - 15*b**4)*sec(c + d*x)/(6*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**3) + (a*tan(c + d*x) + b)*cos(c + d*x)**3/(d*(a + b*tan(c + d*x))**2*(3*a**2 + 3*b**2)) - (-a*(2*a**2 + 9*b**2)*tan(c + d*x) + b*(2*a**2 - 5*b**2))*cos(c + d*x)/(3*d*(a + b*tan(c + d*x))**2*(a**2 + b**2)**2)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_567():
    f = (d*sec(e + f*x))**(sympy.S(7)/2)*(a + b*tan(e + f*x))
    F = -6*a*d**4*elliptic_e(e/2 + f*x/2, 2)/(5*f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))) + 6*a*d**3*sqrt(d*sec(e + f*x))*sin(e + f*x)/(5*f) + 2*a*d*(d*sec(e + f*x))**(sympy.S(5)/2)*sin(e + f*x)/(5*f) + 2*b*(d*sec(e + f*x))**(sympy.S(7)/2)/(7*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_568():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)*(a + b*tan(e + f*x))
    F = 2*a*d**2*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f) + 2*a*d*(d*sec(e + f*x))**(sympy.S(3)/2)*sin(e + f*x)/(3*f) + 2*b*(d*sec(e + f*x))**(sympy.S(5)/2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_569():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x))
    F = -2*a*d**2*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))) + 2*a*d*sqrt(d*sec(e + f*x))*sin(e + f*x)/f + 2*b*(d*sec(e + f*x))**(sympy.S(3)/2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_570():
    f = sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))
    F = 2*a*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/f + 2*b*sqrt(d*sec(e + f*x))/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_571():
    f = (a + b*tan(e + f*x))/sqrt(d*sec(e + f*x))
    F = 2*a*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))) - 2*b/(f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_572():
    f = (a + b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*a*sin(e + f*x)/(3*d*f*sqrt(d*sec(e + f*x))) + 2*a*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*d**2*f) - 2*b/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_573():
    f = (a + b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(5)/2)
    F = 2*a*sin(e + f*x)/(5*d*f*(d*sec(e + f*x))**(sympy.S(3)/2)) + 6*a*elliptic_e(e/2 + f*x/2, 2)/(5*d**2*f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))) - 2*b/(5*f*(d*sec(e + f*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_574():
    f = (a + b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(7)/2)
    F = 2*a*sin(e + f*x)/(7*d*f*(d*sec(e + f*x))**(sympy.S(5)/2)) + 10*a*sin(e + f*x)/(21*d**3*f*sqrt(d*sec(e + f*x))) + 10*a*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(21*d**4*f) - 2*b/(7*f*(d*sec(e + f*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_575():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)*(a + b*tan(e + f*x))**2
    F = 18*a*b*(d*sec(e + f*x))**(sympy.S(5)/2)/(35*f) + 2*b*(d*sec(e + f*x))**(sympy.S(5)/2)*(a + b*tan(e + f*x))/(7*f) + d**2*sqrt(d*sec(e + f*x))*(14*a**2 - 4*b**2)*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(21*f) + d*(d*sec(e + f*x))**(sympy.S(3)/2)*(14*a**2 - 4*b**2)*sin(e + f*x)/(21*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_576():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x))**2
    F = 14*a*b*(d*sec(e + f*x))**(sympy.S(3)/2)/(15*f) + 2*b*(d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x))/(5*f) - d**2*(10*a**2 - 4*b**2)*elliptic_e(e/2 + f*x/2, 2)/(5*f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x))) + d*sqrt(d*sec(e + f*x))*(10*a**2 - 4*b**2)*sin(e + f*x)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_577():
    f = sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))**2
    F = 10*a*b*sqrt(d*sec(e + f*x))/(3*f) + 2*b*sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))/(3*f) + sqrt(d*sec(e + f*x))*(6*a**2 - 4*b**2)*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_578():
    f = (a + b*tan(e + f*x))**2/sqrt(d*sec(e + f*x))
    F = -6*a*b/(f*sqrt(d*sec(e + f*x))) + 2*b*(a + b*tan(e + f*x))/(f*sqrt(d*sec(e + f*x))) + (2*a**2 - 4*b**2)*elliptic_e(e/2 + f*x/2, 2)/(f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_579():
    f = (a + b*tan(e + f*x))**2/(d*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*a*b/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2)) - 2*b*(a + b*tan(e + f*x))/(f*(d*sec(e + f*x))**(sympy.S(3)/2)) + (2*a**2 + 4*b**2)*sin(e + f*x)/(3*d*f*sqrt(d*sec(e + f*x))) + sqrt(d*sec(e + f*x))*(2*a**2 + 4*b**2)*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(3*d**2*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_580():
    f = (a + b*tan(e + f*x))**2/(d*sec(e + f*x))**(sympy.S(5)/2)
    F = -2*a*b/(15*f*(d*sec(e + f*x))**(sympy.S(5)/2)) - 2*b*(a + b*tan(e + f*x))/(3*f*(d*sec(e + f*x))**(sympy.S(5)/2)) + (6*a**2 + 4*b**2)*sin(e + f*x)/(15*d*f*(d*sec(e + f*x))**(sympy.S(3)/2)) + (6*a**2 + 4*b**2)*elliptic_e(e/2 + f*x/2, 2)/(5*d**2*f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_581():
    f = (a + b*tan(e + f*x))**2/(d*sec(e + f*x))**(sympy.S(7)/2)
    F = -6*a*b/(35*f*(d*sec(e + f*x))**(sympy.S(7)/2)) - 2*b*(a + b*tan(e + f*x))/(5*f*(d*sec(e + f*x))**(sympy.S(7)/2)) + (10*a**2 + 4*b**2)*sin(e + f*x)/(35*d*f*(d*sec(e + f*x))**(sympy.S(5)/2)) + (10*a**2 + 4*b**2)*sin(e + f*x)/(21*d**3*f*sqrt(d*sec(e + f*x))) + sqrt(d*sec(e + f*x))*(10*a**2 + 4*b**2)*sqrt(cos(e + f*x))*elliptic_f(e/2 + f*x/2, 2)/(21*d**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_582():
    f = (a + b*tan(e + f*x))**2/(d*sec(e + f*x))**(sympy.S(9)/2)
    F = -10*a*b/(63*f*(d*sec(e + f*x))**(sympy.S(9)/2)) - 2*b*(a + b*tan(e + f*x))/(7*f*(d*sec(e + f*x))**(sympy.S(9)/2)) + (14*a**2 + 4*b**2)*sin(e + f*x)/(63*d*f*(d*sec(e + f*x))**(sympy.S(7)/2)) + (14*a**2 + 4*b**2)*sin(e + f*x)/(45*d**3*f*(d*sec(e + f*x))**(sympy.S(3)/2)) + (14*a**2 + 4*b**2)*elliptic_e(e/2 + f*x/2, 2)/(15*d**4*f*sqrt(d*sec(e + f*x))*sqrt(cos(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_583():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)*(a + b*tan(e + f*x))**3
    F = 2*a*d**2*sqrt(d*sec(e + f*x))*(7*a**2 - 6*b**2)*tan(e + f*x)/(21*f) + 2*a*d**2*sqrt(d*sec(e + f*x))*(7*a**2 - 6*b**2)*elliptic_f(atan(tan(e + f*x))/2, 2)/(21*f*(sec(e + f*x)**2)**(sympy.S(1)/4)) + 2*b*d**2*sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))**2*sec(e + f*x)**2/(9*f) + 2*b*d**2*sqrt(d*sec(e + f*x))*(154*a**2 + 65*a*b*tan(e + f*x) - 28*b**2)*sec(e + f*x)**2/(315*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_584():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x))**3
    F = 2*a*(d*sec(e + f*x))**(sympy.S(3)/2)*(5*a**2 - 6*b**2)*sin(e + f*x)*cos(e + f*x)/(5*f) - 2*a*(d*sec(e + f*x))**(sympy.S(3)/2)*(5*a**2 - 6*b**2)*elliptic_e(atan(tan(e + f*x))/2, 2)/(5*f*(sec(e + f*x)**2)**(sympy.S(3)/4)) + 2*b*(d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x))**2/(7*f) + 2*b*(d*sec(e + f*x))**(sympy.S(3)/2)*(90*a**2 + 33*a*b*tan(e + f*x) - 20*b**2)/(105*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_585():
    f = sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))**3
    F = 2*a*sqrt(d*sec(e + f*x))*(a**2 - 2*b**2)*elliptic_f(atan(tan(e + f*x))/2, 2)/(f*(sec(e + f*x)**2)**(sympy.S(1)/4)) + 2*b*sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))**2/(5*f) + 2*b*sqrt(d*sec(e + f*x))*(14*a**2 + 3*a*b*tan(e + f*x) - 4*b**2)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_586():
    f = (a + b*tan(e + f*x))**3/sqrt(d*sec(e + f*x))
    F = 2*a*(a**2 - 6*b**2)*(sec(e + f*x)**2)**(sympy.S(1)/4)*elliptic_e(atan(tan(e + f*x))/2, 2)/(f*sqrt(d*sec(e + f*x))) - 2*a*(a**2 - 6*b**2)*tan(e + f*x)/(f*sqrt(d*sec(e + f*x))) - 2*b*(6*a**2 + 3*a*b*tan(e + f*x) - 4*b**2)*sec(e + f*x)**2/(3*f*sqrt(d*sec(e + f*x))) - (a + b*tan(e + f*x))**2*(-2*a*tan(e + f*x) + 2*b)/(f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_587():
    f = (a + b*tan(e + f*x))**3/(d*sec(e + f*x))**(sympy.S(3)/2)
    F = 2*a*(a**2 + 6*b**2)*(sec(e + f*x)**2)**(sympy.S(3)/4)*elliptic_f(atan(tan(e + f*x))/2, 2)/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2)) - 2*b*(2*a**2 + a*b*tan(e + f*x) - 4*b**2)*sec(e + f*x)**2/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2)) - (a + b*tan(e + f*x))**2*(-2*a*tan(e + f*x) + 2*b)/(3*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_588():
    f = (a + b*tan(e + f*x))**3/(d*sec(e + f*x))**(sympy.S(5)/2)
    F = 6*a*(a**2 + 2*b**2)*(sec(e + f*x)**2)**(sympy.S(1)/4)*elliptic_e(atan(tan(e + f*x))/2, 2)/(5*d**2*f*sqrt(d*sec(e + f*x))) - 6*a*(a**2 + 2*b**2)*tan(e + f*x)/(5*d**2*f*sqrt(d*sec(e + f*x))) - 2*(a + b*tan(e + f*x))**2*(-a*tan(e + f*x) + b)*cos(e + f*x)**2/(5*d**2*f*sqrt(d*sec(e + f*x))) - (-2*a*(3*a**2 + 5*b**2)*tan(e + f*x) + 4*b*(a**2 + 2*b**2))/(5*d**2*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_589():
    f = (a + b*tan(e + f*x))**3/(d*sec(e + f*x))**(sympy.S(7)/2)
    F = 2*a*(5*a**2 + 6*b**2)*(sec(e + f*x)**2)**(sympy.S(3)/4)*elliptic_f(atan(tan(e + f*x))/2, 2)/(21*d**2*f*(d*sec(e + f*x))**(sympy.S(3)/2)) - 2*(a + b*tan(e + f*x))**2*(-a*tan(e + f*x) + b)*cos(e + f*x)**2/(7*d**2*f*(d*sec(e + f*x))**(sympy.S(3)/2)) - (-2*a*(5*a**2 + 3*b**2)*tan(e + f*x) + 4*b*(3*a**2 + 2*b**2))/(21*d**2*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_590():
    f = (a + b*tan(e + f*x))**3/(d*sec(e + f*x))**(sympy.S(9)/2)
    F = 2*a*(7*a**2 + 6*b**2)*(sec(e + f*x)**2)**(sympy.S(1)/4)*elliptic_e(atan(tan(e + f*x))/2, 2)/(15*d**4*f*sqrt(d*sec(e + f*x))) - 2*(a + b*tan(e + f*x))**2*(-a*tan(e + f*x) + b)*cos(e + f*x)**4/(9*d**4*f*sqrt(d*sec(e + f*x))) - 2*(-a*(7*a**2 + b**2)*tan(e + f*x) + 2*b*(5*a**2 + 2*b**2))*cos(e + f*x)**2/(45*d**4*f*sqrt(d*sec(e + f*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_591():
    f = (a + b*tan(e + f*x))**3/(d*sec(e + f*x))**(sympy.S(11)/2)
    F = 10*a*(3*a**2 + 2*b**2)*(sec(e + f*x)**2)**(sympy.S(3)/4)*elliptic_f(atan(tan(e + f*x))/2, 2)/(77*d**4*f*(d*sec(e + f*x))**(sympy.S(3)/2)) + 10*a*(3*a**2 + 2*b**2)*tan(e + f*x)/(77*d**4*f*(d*sec(e + f*x))**(sympy.S(3)/2)) - 2*(a + b*tan(e + f*x))**2*(-a*tan(e + f*x) + b)*cos(e + f*x)**4/(11*d**4*f*(d*sec(e + f*x))**(sympy.S(3)/2)) - 2*(-a*(9*a**2 - b**2)*tan(e + f*x) + 2*b*(7*a**2 + 2*b**2))*cos(e + f*x)**2/(77*d**4*f*(d*sec(e + f*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_592():
    f = (d*sec(e + f*x))**(sympy.S(7)/2)/(a + b*tan(e + f*x))
    F = ((Integer(2) * (Symbol('d'))**(Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(3) * Symbol('b') * Symbol('f')))**(Integer(-1))) + (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((Symbol('b'))**(Integer(3)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((Symbol('b'))**(Integer(3)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_593():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))
    F = ((Integer(2) * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * (Symbol('d'))**(Integer(2)) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(2) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + ((Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_594():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))
    F = ((sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((sympy.sqrt(Symbol('b')) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Symbol('b') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Symbol('b') * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_595():
    f = sqrt(d*sec(e + f*x))/(a + b*tan(e + f*x))
    F = (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.sqrt(Symbol('b')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((Symbol('a') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + ((Symbol('a') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_596():
    f = 1/(sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x)))
    F = (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * sympy.elliptic_e((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * Symbol('b') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Symbol('a') * Symbol('b') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_597():
    f = 1/((d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x)))
    F = (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * sympy.elliptic_f((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('a') * (Symbol('b'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_598():
    f = 1/((d*sec(e + f*x))**(sympy.S(5)/2)*(a + b*tan(e + f*x)))
    F = (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(2) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(8) * (Symbol('b'))**(Integer(2)))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(8) * (Symbol('b'))**(Integer(2)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * (Symbol('b'))**(Integer(3)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('b'))**(Integer(3)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(2) * (sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(5) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(2) * ((Integer(5) * (Symbol('b'))**(Integer(3))) + (Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(8) * (Symbol('b'))**(Integer(2)))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_599():
    f = (d*sec(e + f*x))**(sympy.S(7)/2)/(a + b*tan(e + f*x))**2
    F = (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * (Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * (((Symbol('b'))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * (Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(3)) * sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Symbol('b') * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_600():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))**2
    F = ((Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + ((Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(2) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (((Symbol('d'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * (((Symbol('b'))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Symbol('b') * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_601():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))**2
    F = ((Symbol('a') * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * sympy.sqrt(Symbol('b')) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((sympy.elliptic_e((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * (((Symbol('a'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (((Symbol('a'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_602():
    f = sqrt(d*sec(e + f*x))/(a + b*tan(e + f*x))**2
    F = ((Integer(-3) * Symbol('a') * sympy.sqrt(Symbol('b')) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * sympy.sqrt(Symbol('b')) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((sympy.elliptic_f((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + ((Integer(3) * (Symbol('a'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_603():
    f = 1/(sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))**2)
    F = ((Integer(5) * Symbol('a') * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * Symbol('a') * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(5) * (Symbol('a'))**(Integer(2)) * Symbol('b') * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(3) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * (((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_604():
    f = 1/((d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x))**2)
    F = ((Integer(-7) * Symbol('a') * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * Symbol('a') * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(3) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(7) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(7) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(3) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(2) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_605():
    f = 1/((d*sec(e + f*x))**(sympy.S(5)/2)*(a + b*tan(e + f*x))**2)
    F = ((Integer(9) * Symbol('a') * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * Symbol('a') * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(3) * ((Integer(2) * (Symbol('a'))**(Integer(4))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(4))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Integer(2) * (Symbol('a'))**(Integer(4))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(4))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(9) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(3)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * ((Integer(2) * (Symbol('a'))**(Integer(4))) + (Integer(10) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(4))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(2) * (sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(5) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * ((Integer(2) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2)))))) + (Integer(-1) * (Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(4) * (Symbol('b'))**(Integer(2)))) * sympy.tan((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_606():
    f = (d*sec(e + f*x))**(sympy.S(7)/2)/(a + b*tan(e + f*x))**3
    F = ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('d'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('d'))**(Integer(2)) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f')))**(Integer(-1)))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(2) * (Symbol('b'))**(Integer(2)))) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (Symbol('b'))**(Integer(3)) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * Symbol('b') * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * (Symbol('d'))**(Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_607():
    f = (d*sec(e + f*x))**(sympy.S(5)/2)/(a + b*tan(e + f*x))**3
    F = ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2)) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + ((((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2)) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(8) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + ((Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.elliptic_f(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (Symbol('b'))**(Integer(2)) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Symbol('a') * ((Symbol('a'))**(Integer(2)) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * (Symbol('d'))**(Integer(2)) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (Symbol('b'))**(Integer(2)) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * (((Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(2) * Symbol('b') * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)))) + ((Symbol('a') * (Symbol('d'))**(Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * Symbol('b') * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_608():
    f = (d*sec(e + f*x))**(sympy.S(3)/2)/(a + b*tan(e + f*x))**3
    F = ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(8) * sympy.sqrt(Symbol('b')) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * sympy.elliptic_e((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Integer(5) * Symbol('a') * sympy.cos((Symbol('e') + (Symbol('f') * x))) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sin((Symbol('e') + (Symbol('f') * x)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f')))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * Symbol('b') * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(5) * (Integer(2))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * Symbol('b') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_609():
    f = sqrt(d*sec(e + f*x))/(a + b*tan(e + f*x))**3
    F = ((Integer(-3) * sympy.sqrt(Symbol('b')) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * sympy.sqrt(Symbol('b')) * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(11) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('a') * sympy.elliptic_f((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * ((Integer(5) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))))**(Integer(-1))) + (Integer(-1) * ((Symbol('b') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(2) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1)))) + (Integer(-1) * ((Integer(7) * Symbol('a') * Symbol('b') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_610():
    f = 1/(sqrt(d*sec(e + f*x))*(a + b*tan(e + f*x))**3)
    F = ((Integer(5) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(5) * (Symbol('b'))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(13) * (Integer(4))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(37) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_e((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(37) * (Symbol('b'))**(Integer(2))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(5) * Symbol('a') * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(5) * Symbol('a') * Symbol('b') * ((Integer(7) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(7) * (Integer(2))**(Integer(-1)))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(5) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(2) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * Symbol('b') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(37) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(4) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_611():
    f = 1/((d*sec(e + f*x))**(sympy.S(3)/2)*(a + b*tan(e + f*x))**3)
    F = ((Integer(-7) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(15) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + (Integer(-1) * ((Integer(7) * (Symbol('b'))**((Integer(5) * (Integer(2))**(Integer(-1)))) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(15) * (Integer(4))**(Integer(-1)))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1)))) + ((Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(69) * (Symbol('b'))**(Integer(2))))) * sympy.elliptic_f((sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x)))) * (Integer(2))**(Integer(-1))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1))))) * ((Integer(12) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(4)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Integer(7) * Symbol('a') * (Symbol('b'))**(Integer(2)) * ((Integer(9) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(3) * (Integer(4))**(Integer(-1)))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(4)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1))))))**(Integer(-1))) + ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(7) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(6) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(3) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Symbol('a') * Symbol('b') * ((Integer(8) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(69) * (Symbol('b'))**(Integer(2))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(12) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * Symbol('f') * ((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))**((Integer(3) * (Integer(2))**(Integer(-1)))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_612():
    f = 1/((d*sec(e + f*x))**(sympy.S(5)/2)*(a + b*tan(e + f*x))**3)
    F = ((Integer(9) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atan(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(17) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(9) * (Symbol('b'))**((Integer(7) * (Integer(2))**(Integer(-1)))) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.atanh(((sympy.sqrt(Symbol('b')) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(4))**(Integer(-1))))**(Integer(-1)))) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(17) * (Integer(4))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(3) * Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(64) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(139) * (Symbol('b'))**(Integer(4))))) * sympy.elliptic_e(((Integer(2))**(Integer(-1)) * sympy.atan(sympy.tan((Symbol('e') + (Symbol('f') * x))))), Integer(2)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))) * ((Integer(20) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(4)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(3) * Symbol('a') * ((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(64) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(139) * (Symbol('b'))**(Integer(4))))) * sympy.tan((Symbol('e') + (Symbol('f') * x)))) * ((Integer(20) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(4)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + (Integer(-1) * ((Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(3)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Integer(-1) * (Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1)))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1)))) + ((Integer(9) * Symbol('a') * (Symbol('b'))**(Integer(3)) * ((Integer(11) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(2) * (Symbol('b'))**(Integer(2))))) * sympy.cot((Symbol('e') + (Symbol('f') * x))) * sympy.Function('EllipticPi')((Symbol('b') * (sympy.sqrt(((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2)))))**(Integer(-1))), sympy.asin(((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1)))), Integer(-1)) * ((sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2)))**((Integer(4))**(Integer(-1))) * sympy.sqrt((Integer(-1) * (sympy.tan((Symbol('e') + (Symbol('f') * x))))**(Integer(2))))) * ((Integer(8) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**((Integer(9) * (Integer(2))**(Integer(-1)))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + ((Integer(3) * Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(4))) + (Integer(28) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(15) * (Symbol('b'))**(Integer(4))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(10) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(3)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(2) * (sympy.cos((Symbol('e') + (Symbol('f') * x))))**(Integer(2)) * (Symbol('b') + (Symbol('a') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))) * ((Integer(5) * ((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))) + ((Integer(3) * Symbol('a') * Symbol('b') * ((Integer(8) * (Symbol('a'))**(Integer(4))) + (Integer(64) * (Symbol('a'))**(Integer(2)) * (Symbol('b'))**(Integer(2))) + (Integer(-1) * (Integer(139) * (Symbol('b'))**(Integer(4))))) * (sympy.sec((Symbol('e') + (Symbol('f') * x))))**(Integer(2))) * ((Integer(20) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(4)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * (Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x)))))))**(Integer(-1))) + (Integer(-1) * ((Integer(2) * ((Symbol('b') * ((Integer(4) * (Symbol('a'))**(Integer(2))) + (Integer(-1) * (Integer(9) * (Symbol('b'))**(Integer(2)))))) + (Integer(-1) * (Symbol('a') * ((Integer(3) * (Symbol('a'))**(Integer(2))) + (Integer(16) * (Symbol('b'))**(Integer(2)))) * sympy.tan((Symbol('e') + (Symbol('f') * x))))))) * ((Integer(5) * (((Symbol('a'))**(Integer(2)) + (Symbol('b'))**(Integer(2))))**(Integer(2)) * (Symbol('d'))**(Integer(2)) * Symbol('f') * sympy.sqrt((Symbol('d') * sympy.sec((Symbol('e') + (Symbol('f') * x))))) * ((Symbol('a') + (Symbol('b') * sympy.tan((Symbol('e') + (Symbol('f') * x))))))**(Integer(2))))**(Integer(-1))))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_613():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)*(a + b*tan(e + f*x))
    F = 3*a*d*(d*sec(e + f*x))**(sympy.S(2)/3)*sin(e + f*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(e + f*x)**2)/(2*f*sqrt(sin(e + f*x)**2)) + 3*b*(d*sec(e + f*x))**(sympy.S(5)/3)/(5*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_614():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)*(a + b*tan(e + f*x))
    F = -3*a*d*sin(e + f*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(e + f*x)**2)/(2*f*(d*sec(e + f*x))**(sympy.S(2)/3)*sqrt(sin(e + f*x)**2)) + 3*b*(d*sec(e + f*x))**(sympy.S(1)/3)/f
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_615():
    f = (a + b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(1)/3)
    F = -3*a*d*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(e + f*x)**2)/(4*f*(d*sec(e + f*x))**(sympy.S(4)/3)*sqrt(sin(e + f*x)**2)) - 3*b/(f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_616():
    f = (a + b*tan(e + f*x))/(d*sec(e + f*x))**(sympy.S(5)/3)
    F = -3*a*d*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(e + f*x)**2)/(8*f*(d*sec(e + f*x))**(sympy.S(8)/3)*sqrt(sin(e + f*x)**2)) - 3*b/(5*f*(d*sec(e + f*x))**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_617():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)*(a + b*tan(e + f*x))**2
    F = 33*a*b*(d*sec(e + f*x))**(sympy.S(5)/3)/(40*f) + 3*b*(d*sec(e + f*x))**(sympy.S(5)/3)*(a + b*tan(e + f*x))/(8*f) + d*(d*sec(e + f*x))**(sympy.S(2)/3)*(24*a**2 - 9*b**2)*sin(e + f*x)*hyper((sympy.S(-1)/3, sympy.S.Half), (sympy.S(2)/3,), cos(e + f*x)**2)/(16*f*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_618():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)*(a + b*tan(e + f*x))**2
    F = 21*a*b*(d*sec(e + f*x))**(sympy.S(1)/3)/(4*f) + 3*b*(d*sec(e + f*x))**(sympy.S(1)/3)*(a + b*tan(e + f*x))/(4*f) - d*(12*a**2 - 9*b**2)*sin(e + f*x)*hyper((sympy.S(1)/3, sympy.S.Half), (sympy.S(4)/3,), cos(e + f*x)**2)/(8*f*(d*sec(e + f*x))**(sympy.S(2)/3)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_619():
    f = (a + b*tan(e + f*x))**2/(d*sec(e + f*x))**(sympy.S(1)/3)
    F = -15*a*b/(2*f*(d*sec(e + f*x))**(sympy.S(1)/3)) + 3*b*(a + b*tan(e + f*x))/(2*f*(d*sec(e + f*x))**(sympy.S(1)/3)) - d*(6*a**2 - 9*b**2)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(2)/3), (sympy.S(5)/3,), cos(e + f*x)**2)/(8*f*(d*sec(e + f*x))**(sympy.S(4)/3)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_620():
    f = (a + b*tan(e + f*x))**2/(d*sec(e + f*x))**(sympy.S(5)/3)
    F = 3*a*b/(10*f*(d*sec(e + f*x))**(sympy.S(5)/3)) - 3*b*(a + b*tan(e + f*x))/(2*f*(d*sec(e + f*x))**(sympy.S(5)/3)) - d*(6*a**2 + 9*b**2)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S(4)/3), (sympy.S(7)/3,), cos(e + f*x)**2)/(16*f*(d*sec(e + f*x))**(sympy.S(8)/3)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_621():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)/(a + b*tan(e + f*x))
    F = (d*sec(e + f*x))**(sympy.S(5)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) - (d*sec(e + f*x))**(sympy.S(5)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) + sqrt(3)*(d*sec(e + f*x))**(sympy.S(5)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(2*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) + sqrt(3)*(d*sec(e + f*x))**(sympy.S(5)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(2*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) - (d*sec(e + f*x))**(sympy.S(5)/3)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) + (d*sec(e + f*x))**(sympy.S(5)/3)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(1)/6, 1, sympy.S(3)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(a*f*(sec(e + f*x)**2)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_622():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)/(a + b*tan(e + f*x))
    F = b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*f*(a**2 + b**2)**(sympy.S(5)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*f*(a**2 + b**2)**(sympy.S(5)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - sqrt(3)*b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(2*f*(a**2 + b**2)**(sympy.S(5)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - sqrt(3)*b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(2*f*(a**2 + b**2)**(sympy.S(5)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(f*(a**2 + b**2)**(sympy.S(5)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) + (d*sec(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(5)/6, 1, sympy.S(3)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(a*f*(sec(e + f*x)**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_623():
    f = 1/((d*sec(e + f*x))**(sympy.S(1)/3)*(a + b*tan(e + f*x)))
    F = b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(7)/6)) - b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(7)/6)) + sqrt(3)*b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(2*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(7)/6)) + sqrt(3)*b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(2*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(7)/6)) - b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(7)/6)) + 3*b/(f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)) + (sec(e + f*x)**2)**(sympy.S(1)/6)*tan(e + f*x)*appellf1(sympy.S.Half, 1, sympy.S(7)/6, sympy.S(3)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(a*f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_624():
    f = 1/((d*sec(e + f*x))**(sympy.S(5)/3)*(a + b*tan(e + f*x)))
    F = b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(11)/6)) - b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(4*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(11)/6)) - sqrt(3)*b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(2*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(11)/6)) - sqrt(3)*b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(2*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(11)/6)) - b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(11)/6)) + 3*b/(f*(d*sec(e + f*x))**(sympy.S(5)/3)*(5*a**2 + 5*b**2)) + (sec(e + f*x)**2)**(sympy.S(5)/6)*tan(e + f*x)*appellf1(sympy.S.Half, 1, sympy.S(11)/6, sympy.S(3)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(a*f*(d*sec(e + f*x))**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_625():
    f = (d*sec(e + f*x))**(sympy.S(5)/3)/(a + b*tan(e + f*x))**2
    F = -a*b*(d*sec(e + f*x))**(sympy.S(5)/3)/(f*(a**2 + b**2)*(a**2 - b**2*tan(e + f*x)**2)) + a*(d*sec(e + f*x))**(sympy.S(5)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(7)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) - a*(d*sec(e + f*x))**(sympy.S(5)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(7)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) + sqrt(3)*a*(d*sec(e + f*x))**(sympy.S(5)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(6*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(7)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) + sqrt(3)*a*(d*sec(e + f*x))**(sympy.S(5)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(6*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(7)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) - a*(d*sec(e + f*x))**(sympy.S(5)/3)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(3*b**(sympy.S(2)/3)*f*(a**2 + b**2)**(sympy.S(7)/6)*(sec(e + f*x)**2)**(sympy.S(5)/6)) + (d*sec(e + f*x))**(sympy.S(5)/3)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(1)/6, 2, sympy.S(3)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(a**2*f*(sec(e + f*x)**2)**(sympy.S(5)/6)) + b**2*(d*sec(e + f*x))**(sympy.S(5)/3)*tan(e + f*x)**3*appellf1(sympy.S(3)/2, sympy.S(1)/6, 2, sympy.S(5)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(3*a**4*f*(sec(e + f*x)**2)**(sympy.S(5)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_626():
    f = (d*sec(e + f*x))**(sympy.S(1)/3)/(a + b*tan(e + f*x))**2
    F = 5*a*b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*f*(a**2 + b**2)**(sympy.S(11)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - 5*a*b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*f*(a**2 + b**2)**(sympy.S(11)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - 5*sqrt(3)*a*b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(6*f*(a**2 + b**2)**(sympy.S(11)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - 5*sqrt(3)*a*b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(6*f*(a**2 + b**2)**(sympy.S(11)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - 5*a*b**(sympy.S(2)/3)*(d*sec(e + f*x))**(sympy.S(1)/3)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(3*f*(a**2 + b**2)**(sympy.S(11)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6)) - a*b*(d*sec(e + f*x))**(sympy.S(1)/3)/(f*(a**2 + b**2)*(a**2 - b**2*tan(e + f*x)**2)) + (d*sec(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(5)/6, 2, sympy.S(3)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(a**2*f*(sec(e + f*x)**2)**(sympy.S(1)/6)) + b**2*(d*sec(e + f*x))**(sympy.S(1)/3)*tan(e + f*x)**3*appellf1(sympy.S(3)/2, sympy.S(5)/6, 2, sympy.S(5)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(3*a**4*f*(sec(e + f*x)**2)**(sympy.S(1)/6))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_627():
    f = 1/((d*sec(e + f*x))**(sympy.S(1)/3)*(a + b*tan(e + f*x))**2)
    F = 7*a*b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(13)/6)) - 7*a*b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(13)/6)) + 7*sqrt(3)*a*b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(6*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(13)/6)) + 7*sqrt(3)*a*b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(6*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(13)/6)) - 7*a*b**(sympy.S(4)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(3*f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(13)/6)) - a*b/(f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)*(a**2 - b**2*tan(e + f*x)**2)) + 7*a*b/(f*(d*sec(e + f*x))**(sympy.S(1)/3)*(a**2 + b**2)**2) + (sec(e + f*x)**2)**(sympy.S(1)/6)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(7)/6, 2, sympy.S(3)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(a**2*f*(d*sec(e + f*x))**(sympy.S(1)/3)) + b**2*(sec(e + f*x)**2)**(sympy.S(1)/6)*tan(e + f*x)**3*appellf1(sympy.S(3)/2, sympy.S(7)/6, 2, sympy.S(5)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(3*a**4*f*(d*sec(e + f*x))**(sympy.S(1)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_628():
    f = 1/((d*sec(e + f*x))**(sympy.S(5)/3)*(a + b*tan(e + f*x))**2)
    F = 11*a*b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) - b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(17)/6)) - 11*a*b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*log(b**(sympy.S(2)/3)*(sec(e + f*x)**2)**(sympy.S(1)/3) + b**(sympy.S(1)/3)*(a**2 + b**2)**(sympy.S(1)/6)*(sec(e + f*x)**2)**(sympy.S(1)/6) + (a**2 + b**2)**(sympy.S(1)/3))/(12*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(17)/6)) - 11*sqrt(3)*a*b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) - sqrt(3)/3)/(6*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(17)/6)) - 11*sqrt(3)*a*b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*atan(2*sqrt(3)*b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(3*(a**2 + b**2)**(sympy.S(1)/6)) + sqrt(3)/3)/(6*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(17)/6)) - 11*a*b**(sympy.S(8)/3)*(sec(e + f*x)**2)**(sympy.S(5)/6)*atanh(b**(sympy.S(1)/3)*(sec(e + f*x)**2)**(sympy.S(1)/6)/(a**2 + b**2)**(sympy.S(1)/6))/(3*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**(sympy.S(17)/6)) - a*b/(f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)*(a**2 - b**2*tan(e + f*x)**2)) + 11*a*b/(5*f*(d*sec(e + f*x))**(sympy.S(5)/3)*(a**2 + b**2)**2) + (sec(e + f*x)**2)**(sympy.S(5)/6)*tan(e + f*x)*appellf1(sympy.S.Half, sympy.S(11)/6, 2, sympy.S(3)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(a**2*f*(d*sec(e + f*x))**(sympy.S(5)/3)) + b**2*(sec(e + f*x)**2)**(sympy.S(5)/6)*tan(e + f*x)**3*appellf1(sympy.S(3)/2, sympy.S(11)/6, 2, sympy.S(5)/2, -tan(e + f*x)**2, b**2*tan(e + f*x)**2/a**2)/(3*a**4*f*(d*sec(e + f*x))**(sympy.S(5)/3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_629():
    f = (d*sec(e + f*x))**m*(a + b*tan(e + f*x))**3
    F = -a*(d*sec(e + f*x))**m*(-a**2*(m + 1) + 3*b**2)*tan(e + f*x)*hyper((sympy.S.Half, 1 - m/2), (sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(m + 1)*(sec(e + f*x)**2)**(m/2)) + b*(d*sec(e + f*x))**m*(a + b*tan(e + f*x))**2/(f*(m + 2)) - b*(d*sec(e + f*x))**m*(-a*b*m*(m + 4)*tan(e + f*x) + (2*m + 2)*(-a**2*(m + 3) + b**2))/(f*m*(m**2 + 3*m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_630():
    f = (d*sec(e + f*x))**m*(a + b*tan(e + f*x))**2
    F = a*b*(d*sec(e + f*x))**m*(m + 2)/(f*m*(m + 1)) + b*(d*sec(e + f*x))**m*(a + b*tan(e + f*x))/(f*(m + 1)) + d*(d*sec(e + f*x))**(m - 1)*(-a**2*(m + 1) + b**2)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(e + f*x)**2)/(f*(1 - m)*(m + 1)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_631():
    f = (d*sec(e + f*x))**m*(a + b*tan(e + f*x))
    F = -a*d*(d*sec(e + f*x))**(m - 1)*sin(e + f*x)*hyper((sympy.S.Half, sympy.S.Half - m/2), (sympy.S(3)/2 - m/2,), cos(e + f*x)**2)/(f*(1 - m)*sqrt(sin(e + f*x)**2)) + b*(d*sec(e + f*x))**m/(f*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_632():
    f = (d*sec(e + f*x))**m/(a + b*tan(e + f*x))
    F = -b*(d*sec(e + f*x))**m*hyper((1, m/2), (m/2 + 1,), b**2*sec(e + f*x)**2/(a**2 + b**2))/(f*m*(a**2 + b**2)) + (d*sec(e + f*x))**m*tan(e + f*x)*appellf1(sympy.S.Half, 1, 1 - m/2, sympy.S(3)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(a*f*(sec(e + f*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_633():
    f = (d*sec(e + f*x))**m/(a + b*tan(e + f*x))**2
    F = -2*a*b*(d*sec(e + f*x))**m*hyper((2, m/2), (m/2 + 1,), b**2*sec(e + f*x)**2/(a**2 + b**2))/(f*m*(a**2 + b**2)**2) + (d*sec(e + f*x))**m*tan(e + f*x)*appellf1(sympy.S.Half, 2, 1 - m/2, sympy.S(3)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(a**2*f*(sec(e + f*x)**2)**(m/2)) + b**2*(d*sec(e + f*x))**m*tan(e + f*x)**3*appellf1(sympy.S(3)/2, 2, 1 - m/2, sympy.S(5)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(3*a**4*f*(sec(e + f*x)**2)**(m/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_634():
    f = (a + b*tan(c + d*x))**n*sec(c + d*x)**6
    F = -4*a*(a + b*tan(c + d*x))**(n + 2)*(a**2 + b**2)/(b**5*d*(n + 2)) - 4*a*(a + b*tan(c + d*x))**(n + 4)/(b**5*d*(n + 4)) + (a + b*tan(c + d*x))**(n + 1)*(a**2 + b**2)**2/(b**5*d*(n + 1)) + (a + b*tan(c + d*x))**(n + 3)*(6*a**2 + 2*b**2)/(b**5*d*(n + 3)) + (a + b*tan(c + d*x))**(n + 5)/(b**5*d*(n + 5))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_635():
    f = (a + b*tan(c + d*x))**n*sec(c + d*x)**4
    F = -2*a*(a + b*tan(c + d*x))**(n + 2)/(b**3*d*(n + 2)) + (a + b*tan(c + d*x))**(n + 1)*(a**2 + b**2)/(b**3*d*(n + 1)) + (a + b*tan(c + d*x))**(n + 3)/(b**3*d*(n + 3))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_636():
    f = (a + b*tan(c + d*x))**n*sec(c + d*x)**2
    F = (a + b*tan(c + d*x))**(n + 1)/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_637():
    f = (a + b*tan(c + d*x))**n*cos(c + d*x)**2
    F = b*(a + b*tan(c + d*x))**(n + 1)*(a*n + sqrt(-b**2)*(a**2/b**2 - n + 1))*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(d*(a + sqrt(-b**2))*(4*a**2 + 4*b**2)*(n + 1)) + (a + b*tan(c + d*x))**(n + 1)*(a*tan(c + d*x) + b)*cos(c + d*x)**2/(d*(2*a**2 + 2*b**2)) - (a + b*tan(c + d*x))**(n + 1)*(-a*n + sqrt(-b**2)*(a**2/b**2 - n + 1))*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(b*d*(a - sqrt(-b**2))*(n + 1)*(4*a**2/b**2 + 4))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_638():
    f = (a + b*tan(c + d*x))**n*cos(c + d*x)**4
    F = b*(a + b*tan(c + d*x))**(n + 1)*(a**2*(n + 1) + a*b*(3*a**2/b**2 - 2*n + 5)*tan(c + d*x) + b**2*(3 - n))*cos(c + d*x)**2/(8*d*(a**2 + b**2)**2) + b*(a + b*tan(c + d*x))**(n + 1)*(a*n*(3*a**2/b**2 - 2*n + 5)/b**2 + sqrt(-b**2)*(3*a**4 + a**2*b**2*(-n**2 - 2*n + 6) + b**4*(n**2 - 4*n + 3))/b**6)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(16*d*(a + sqrt(-b**2))*(n + 1)*(a**2/b**2 + 1)**2) + b*(a + b*tan(c + d*x))**(n + 1)*(a*n*(3*a**2/b**2 - 2*n + 5)/b**2 - sqrt(-b**2)*(3*a**4 + a**2*b**2*(-n**2 - 2*n + 6) + b**4*(n**2 - 4*n + 3))/b**6)*hyper((1, n + 1), (n + 2,), (a + b*tan(c + d*x))/(a - sqrt(-b**2)))/(16*d*(a - sqrt(-b**2))*(n + 1)*(a**2/b**2 + 1)**2) + (a + b*tan(c + d*x))**(n + 1)*(a*tan(c + d*x) + b)*cos(c + d*x)**4/(d*(4*a**2 + 4*b**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_639():
    f = (a + b*tan(c + d*x))**n*sec(c + d*x)**3
    F = (a + b*tan(c + d*x))**(n + 1)*appellf1(n + 1, sympy.S(-1)/2, sympy.S(-1)/2, n + 2, (a + b*tan(c + d*x))/(a - sqrt(-b**2)), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))*sec(c + d*x)/(b*d*sqrt(1 - (a + b*tan(c + d*x))/(a - sqrt(-b**2)))*sqrt(1 - (a + b*tan(c + d*x))/(a + sqrt(-b**2)))*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_640():
    f = (a + b*tan(c + d*x))**n*sec(c + d*x)
    F = sqrt(1 - (a + b*tan(c + d*x))/(a - sqrt(-b**2)))*sqrt(1 - (a + b*tan(c + d*x))/(a + sqrt(-b**2)))*(a + b*tan(c + d*x))**(n + 1)*cos(c + d*x)*appellf1(n + 1, sympy.S.Half, sympy.S.Half, n + 2, (a + b*tan(c + d*x))/(a - sqrt(-b**2)), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_641():
    f = (a + b*tan(c + d*x))**n*cos(c + d*x)
    F = (1 - (a + b*tan(c + d*x))/(a - sqrt(-b**2)))**(sympy.S(3)/2)*(1 - (a + b*tan(c + d*x))/(a + sqrt(-b**2)))**(sympy.S(3)/2)*(a + b*tan(c + d*x))**(n + 1)*cos(c + d*x)**3*appellf1(n + 1, sympy.S(3)/2, sympy.S(3)/2, n + 2, (a + b*tan(c + d*x))/(a - sqrt(-b**2)), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_642():
    f = (a + b*tan(c + d*x))**n*cos(c + d*x)**3
    F = (1 - (a + b*tan(c + d*x))/(a - sqrt(-b**2)))**(sympy.S(5)/2)*(1 - (a + b*tan(c + d*x))/(a + sqrt(-b**2)))**(sympy.S(5)/2)*(a + b*tan(c + d*x))**(n + 1)*cos(c + d*x)**5*appellf1(n + 1, sympy.S(5)/2, sympy.S(5)/2, n + 2, (a + b*tan(c + d*x))/(a - sqrt(-b**2)), (a + b*tan(c + d*x))/(a + sqrt(-b**2)))/(b*d*(n + 1))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_643():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a)
    F = 10*a*(e*cos(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)*sec(c + d*x)**2/(21*d) + 2*a*(e*cos(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)/(7*d) - 2*I*a*(e*cos(c + d*x))**(sympy.S(7)/2)/(7*d) + 10*a*(e*cos(c + d*x))**(sympy.S(7)/2)*elliptic_f(c/2 + d*x/2, 2)/(21*d*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_644():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)
    F = 2*a*(e*cos(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(5*d) - 2*I*a*(e*cos(c + d*x))**(sympy.S(5)/2)/(5*d) + 6*a*(e*cos(c + d*x))**(sympy.S(5)/2)*elliptic_e(c/2 + d*x/2, 2)/(5*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_645():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)
    F = 2*a*(e*cos(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(3*d) - 2*I*a*(e*cos(c + d*x))**(sympy.S(3)/2)/(3*d) + 2*a*(e*cos(c + d*x))**(sympy.S(3)/2)*elliptic_f(c/2 + d*x/2, 2)/(3*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_646():
    f = sqrt(e*cos(c + d*x))*(I*a*tan(c + d*x) + a)
    F = -2*I*a*sqrt(e*cos(c + d*x))/d + 2*a*sqrt(e*cos(c + d*x))*elliptic_e(c/2 + d*x/2, 2)/(d*sqrt(cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_647():
    f = (I*a*tan(c + d*x) + a)/sqrt(e*cos(c + d*x))
    F = 2*a*sqrt(cos(c + d*x))*elliptic_f(c/2 + d*x/2, 2)/(d*sqrt(e*cos(c + d*x))) + 2*I*a/(d*sqrt(e*cos(c + d*x)))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_648():
    f = (I*a*tan(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = 2*a*sin(c + d*x)*cos(c + d*x)/(3*d*(e*cos(c + d*x))**(sympy.S(5)/2)) + 2*a*cos(c + d*x)**(sympy.S(5)/2)*elliptic_f(c/2 + d*x/2, 2)/(3*d*(e*cos(c + d*x))**(sympy.S(5)/2)) + 2*I*a/(5*d*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_649():
    f = (I*a*tan(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 6*a*sin(c + d*x)*cos(c + d*x)**3/(5*d*(e*cos(c + d*x))**(sympy.S(7)/2)) + 2*a*sin(c + d*x)*cos(c + d*x)/(5*d*(e*cos(c + d*x))**(sympy.S(7)/2)) - 6*a*cos(c + d*x)**(sympy.S(7)/2)*elliptic_e(c/2 + d*x/2, 2)/(5*d*(e*cos(c + d*x))**(sympy.S(7)/2)) + 2*I*a/(7*d*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_650():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)/(I*a*tan(c + d*x) + a)**2
    F = 4*I*(e*cos(c + d*x))**(sympy.S(7)/2)*cos(c + d*x)**2/(15*d*(I*a**2*tan(c + d*x) + a**2)) + 2*(e*cos(c + d*x))**(sympy.S(7)/2)*sin(c + d*x)*cos(c + d*x)/(15*a**2*d) + 2*(e*cos(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)*sec(c + d*x)**2/(7*a**2*d) + 6*(e*cos(c + d*x))**(sympy.S(7)/2)*tan(c + d*x)/(35*a**2*d) + 2*(e*cos(c + d*x))**(sympy.S(7)/2)*elliptic_f(c/2 + d*x/2, 2)/(7*a**2*d*cos(c + d*x)**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_651():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/(I*a*tan(c + d*x) + a)**2
    F = 4*I*(e*cos(c + d*x))**(sympy.S(5)/2)*cos(c + d*x)**2/(13*d*(I*a**2*tan(c + d*x) + a**2)) + 2*(e*cos(c + d*x))**(sympy.S(5)/2)*sin(c + d*x)*cos(c + d*x)/(13*a**2*d) + 14*(e*cos(c + d*x))**(sympy.S(5)/2)*tan(c + d*x)/(65*a**2*d) + 42*(e*cos(c + d*x))**(sympy.S(5)/2)*elliptic_e(c/2 + d*x/2, 2)/(65*a**2*d*cos(c + d*x)**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_652():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/(I*a*tan(c + d*x) + a)**2
    F = 4*I*(e*cos(c + d*x))**(sympy.S(3)/2)*cos(c + d*x)**2/(11*d*(I*a**2*tan(c + d*x) + a**2)) + 2*(e*cos(c + d*x))**(sympy.S(3)/2)*sin(c + d*x)*cos(c + d*x)/(11*a**2*d) + 10*(e*cos(c + d*x))**(sympy.S(3)/2)*tan(c + d*x)/(33*a**2*d) + 10*(e*cos(c + d*x))**(sympy.S(3)/2)*elliptic_f(c/2 + d*x/2, 2)/(33*a**2*d*cos(c + d*x)**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_653():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*(I*a*tan(c + d*x) + a)**2)
    F = 4*I*cos(c + d*x)**2/(5*d*(e*cos(c + d*x))**(sympy.S(3)/2)*(I*a**2*tan(c + d*x) + a**2)) + 2*cos(c + d*x)**(sympy.S(3)/2)*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d*(e*cos(c + d*x))**(sympy.S(3)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_654():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*(I*a*tan(c + d*x) + a)**2)
    F = 4*I*cos(c + d*x)**2/(3*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(I*a**2*tan(c + d*x) + a**2)) - 2*cos(c + d*x)**(sympy.S(5)/2)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_655():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*(I*a*tan(c + d*x) + a)**2)
    F = 4*I*cos(c + d*x)**2/(d*(e*cos(c + d*x))**(sympy.S(7)/2)*(I*a**2*tan(c + d*x) + a**2)) - 6*sin(c + d*x)*cos(c + d*x)**3/(a**2*d*(e*cos(c + d*x))**(sympy.S(7)/2)) + 6*cos(c + d*x)**(sympy.S(7)/2)*elliptic_e(c/2 + d*x/2, 2)/(a**2*d*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_656():
    f = 1/((e*cos(c + d*x))**(sympy.S(9)/2)*(I*a*tan(c + d*x) + a)**2)
    F = -4*I*cos(c + d*x)**2/(d*(e*cos(c + d*x))**(sympy.S(9)/2)*(I*a**2*tan(c + d*x) + a**2)) + 10*sin(c + d*x)*cos(c + d*x)**3/(3*a**2*d*(e*cos(c + d*x))**(sympy.S(9)/2)) + 10*cos(c + d*x)**(sympy.S(9)/2)*elliptic_f(c/2 + d*x/2, 2)/(3*a**2*d*(e*cos(c + d*x))**(sympy.S(9)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_657():
    f = 1/((e*cos(c + d*x))**(sympy.S(11)/2)*(I*a*tan(c + d*x) + a)**2)
    F = -4*I*cos(c + d*x)**2/(3*d*(e*cos(c + d*x))**(sympy.S(11)/2)*(I*a**2*tan(c + d*x) + a**2)) + 14*sin(c + d*x)*cos(c + d*x)**5/(5*a**2*d*(e*cos(c + d*x))**(sympy.S(11)/2)) + 14*sin(c + d*x)*cos(c + d*x)**3/(15*a**2*d*(e*cos(c + d*x))**(sympy.S(11)/2)) - 14*cos(c + d*x)**(sympy.S(11)/2)*elliptic_e(c/2 + d*x/2, 2)/(5*a**2*d*(e*cos(c + d*x))**(sympy.S(11)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_658():
    f = (e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a)
    F = 32*I*a*(e*cos(c + d*x))**(sympy.S(7)/2)*sec(c + d*x)**4/(35*d*sqrt(I*a*tan(c + d*x) + a)) + 12*I*a*(e*cos(c + d*x))**(sympy.S(7)/2)*sec(c + d*x)**2/(35*d*sqrt(I*a*tan(c + d*x) + a)) - 16*I*(e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**2/(35*d) - 2*I*(e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a)/(7*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_659():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)
    F = 8*I*a*(e*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2/(15*d*sqrt(I*a*tan(c + d*x) + a)) - 16*I*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**2/(15*d) - 2*I*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)/(5*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_660():
    f = sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a)
    F = -2*I*sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/d
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_661():
    f = sqrt(I*a*tan(c + d*x) + a)/sqrt(e*cos(c + d*x))
    F = -sqrt(2)*I*sqrt(a)*log(-sqrt(2)*sqrt(a)*sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a) + a*sqrt(e) + sqrt(e)*(I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*d*sqrt(e)) + sqrt(2)*I*sqrt(a)*log(sqrt(2)*sqrt(a)*sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a) + a*sqrt(e) + sqrt(e)*(I*a*tan(c + d*x) + a)*cos(c + d*x))/(2*d*sqrt(e)) + sqrt(2)*I*sqrt(a)*atan(1 - sqrt(2)*sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e)))/(d*sqrt(e)) - sqrt(2)*I*sqrt(a)*atan(1 + sqrt(2)*sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e)))/(d*sqrt(e))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_662():
    f = sqrt(I*a*tan(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(5)/2)
    F = -3*sqrt(2)*I*sqrt(a)*e**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(16*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2)) + 3*sqrt(2)*I*sqrt(a)*e**(sympy.S(5)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(16*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2)) + 3*sqrt(2)*I*sqrt(a)*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(8*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2)) - 3*sqrt(2)*I*sqrt(a)*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(8*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2)) + I*a/(2*d*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)) - 3*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**2/(4*d*(e*cos(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_663():
    f = sqrt(I*a*tan(c + d*x) + a)/(e*cos(c + d*x))**(sympy.S(7)/2)
    F = 5*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(7)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(32*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 5*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(7)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(32*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 5*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(16*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 5*sqrt(2)*I*a**(sympy.S(3)/2)*e**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(16*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 5*I*a*cos(c + d*x)**2/(8*d*(e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a)) + I*a/(3*d*(e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a)) - 5*I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**2/(12*d*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_664():
    f = (e*cos(c + d*x))**(sympy.S(5)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = 16*I*(e*cos(c + d*x))**(sympy.S(5)/2)*sec(c + d*x)**2/(35*d*sqrt(I*a*tan(c + d*x) + a)) + 2*I*(e*cos(c + d*x))**(sympy.S(5)/2)/(7*d*sqrt(I*a*tan(c + d*x) + a)) - 32*I*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)*sec(c + d*x)**2/(35*a*d) - 12*I*(e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a)/(35*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_665():
    f = (e*cos(c + d*x))**(sympy.S(3)/2)/sqrt(I*a*tan(c + d*x) + a)
    F = 16*I*(e*cos(c + d*x))**(sympy.S(3)/2)*sec(c + d*x)**2/(15*d*sqrt(I*a*tan(c + d*x) + a)) + 2*I*(e*cos(c + d*x))**(sympy.S(3)/2)/(5*d*sqrt(I*a*tan(c + d*x) + a)) - 8*I*(e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a)/(15*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_666():
    f = sqrt(e*cos(c + d*x))/sqrt(I*a*tan(c + d*x) + a)
    F = 2*I*sqrt(e*cos(c + d*x))/(3*d*sqrt(I*a*tan(c + d*x) + a)) - 4*I*sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a)/(3*a*d)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_667():
    f = 1/(sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a))
    F = 2*I/(d*sqrt(e*cos(c + d*x))*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_668():
    f = 1/((e*cos(c + d*x))**(sympy.S(3)/2)*sqrt(I*a*tan(c + d*x) + a))
    F = sqrt(2)*I*sqrt(a)*log(-sqrt(2)*sqrt(a)*sqrt(e*cos(c + d*x))*sqrt(-I*a*tan(c + d*x) + a) + a*sqrt(e) + sqrt(e)*(-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*d*e**(sympy.S(3)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*sqrt(a)*log(sqrt(2)*sqrt(a)*sqrt(e*cos(c + d*x))*sqrt(-I*a*tan(c + d*x) + a) + a*sqrt(e) + sqrt(e)*(-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(2*d*e**(sympy.S(3)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - sqrt(2)*I*sqrt(a)*atan(1 - sqrt(2)*sqrt(e*cos(c + d*x))*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e)))*sec(c + d*x)/(d*e**(sympy.S(3)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + sqrt(2)*I*sqrt(a)*atan(1 + sqrt(2)*sqrt(e*cos(c + d*x))*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e)))*sec(c + d*x)/(d*e**(sympy.S(3)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_669():
    f = 1/((e*cos(c + d*x))**(sympy.S(5)/2)*sqrt(I*a*tan(c + d*x) + a))
    F = -I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**2/(a*d*(e*cos(c + d*x))**(sympy.S(5)/2)) - sqrt(2)*I*e**(sympy.S(5)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*sqrt(a)*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2)) + sqrt(2)*I*e**(sympy.S(5)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (I*a*tan(c + d*x) + a)*cos(c + d*x))/(4*sqrt(a)*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2)) + sqrt(2)*I*e**(sympy.S(5)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*sqrt(a)*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2)) - sqrt(2)*I*e**(sympy.S(5)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))/(2*sqrt(a)*d*(e*cos(c + d*x))**(sympy.S(5)/2)*(e*sec(c + d*x))**(sympy.S(5)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_670():
    f = 1/((e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a))
    F = 3*sqrt(2)*I*sqrt(a)*e**(sympy.S(7)/2)*log(-sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(16*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 3*sqrt(2)*I*sqrt(a)*e**(sympy.S(7)/2)*log(sqrt(2)*sqrt(a)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/sqrt(e*sec(c + d*x)) + a + (-I*a*tan(c + d*x) + a)*cos(c + d*x))*sec(c + d*x)/(16*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) - 3*sqrt(2)*I*sqrt(a)*e**(sympy.S(7)/2)*atan(1 - sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(8*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 3*sqrt(2)*I*sqrt(a)*e**(sympy.S(7)/2)*atan(1 + sqrt(2)*sqrt(e)*sqrt(-I*a*tan(c + d*x) + a)/(sqrt(a)*sqrt(e*sec(c + d*x))))*sec(c + d*x)/(8*d*(e*cos(c + d*x))**(sympy.S(7)/2)*(e*sec(c + d*x))**(sympy.S(7)/2)*sqrt(-I*a*tan(c + d*x) + a)*sqrt(I*a*tan(c + d*x) + a)) + 3*I*cos(c + d*x)**2/(4*d*(e*cos(c + d*x))**(sympy.S(7)/2)*sqrt(I*a*tan(c + d*x) + a)) - I*sqrt(I*a*tan(c + d*x) + a)*cos(c + d*x)**2/(2*a*d*(e*cos(c + d*x))**(sympy.S(7)/2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_671():
    f = (e*cos(c + d*x))**m*(I*a*tan(c + d*x) + a)**n
    F = -2**(-m/2 + n)*I*(e*cos(c + d*x))**m*(I*tan(c + d*x) + 1)**(m/2 - n)*(I*a*tan(c + d*x) + a)**n*hyper((-m/2, m/2 - n + 1), (1 - m/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_672():
    f = (e*cos(c + d*x))**m*(I*a*tan(c + d*x) + a)**2
    F = -2**(2 - m/2)*I*a**2*(e*cos(c + d*x))**m*(I*tan(c + d*x) + 1)**(m/2)*hyper((-m/2, m/2 - 1), (1 - m/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_673():
    f = (e*cos(c + d*x))**m*(I*a*tan(c + d*x) + a)
    F = -2**(1 - m/2)*I*a*(e*cos(c + d*x))**m*(I*tan(c + d*x) + 1)**(m/2)*hyper((-m/2, m/2), (1 - m/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_674():
    f = (e*cos(c + d*x))**m/(I*a*tan(c + d*x) + a)
    F = -2**(-m/2 - 1)*I*(e*cos(c + d*x))**m*(I*tan(c + d*x) + 1)**(m/2)*hyper((-m/2, m/2 + 2), (1 - m/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(a*d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_675():
    f = (e*cos(c + d*x))**m/(I*a*tan(c + d*x) + a)**2
    F = -2**(-m/2 - 2)*I*(e*cos(c + d*x))**m*(I*tan(c + d*x) + 1)**(m/2)*hyper((-m/2, m/2 + 3), (1 - m/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(a**2*d*m)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_676():
    f = (e*cos(c + d*x))**m*sqrt(I*a*tan(c + d*x) + a)
    F = -2**(sympy.S.Half - m/2)*I*a*(e*cos(c + d*x))**m*(I*tan(c + d*x) + 1)**(m/2 + sympy.S.Half)*hyper((-m/2, m/2 + sympy.S.Half), (1 - m/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_677():
    f = (e*cos(c + d*x))**m/sqrt(I*a*tan(c + d*x) + a)
    F = -2**(-m/2 + sympy.S(-1)/2)*I*(e*cos(c + d*x))**m*(I*tan(c + d*x) + 1)**(m/2 + sympy.S.Half)*hyper((-m/2, m/2 + sympy.S(3)/2), (1 - m/2,), -I*tan(c + d*x)/2 + sympy.S.Half)/(d*m*sqrt(I*a*tan(c + d*x) + a))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_678():
    f = (d*cos(e + f*x))**m*(a + b*tan(e + f*x))**3
    F = -a*(d*cos(e + f*x))**m*(-a**2*(1 - m) + 3*b**2)*(sec(e + f*x)**2)**(m/2)*tan(e + f*x)*hyper((sympy.S.Half, m/2 + 1), (sympy.S(3)/2,), -tan(e + f*x)**2)/(f*(1 - m)) + b*(d*cos(e + f*x))**m*(a + b*tan(e + f*x))**2/(f*(2 - m)) + b*(d*cos(e + f*x))**m*(a*b*m*(4 - m)*tan(e + f*x) + (1 - m)*(-2*a**2*(3 - m) + 2*b**2))/(f*m*(m**2 - 3*m + 2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_679():
    f = (d*cos(e + f*x))**m*(a + b*tan(e + f*x))**2
    F = -a*b*(d*cos(e + f*x))**m*(2 - m)/(f*m*(1 - m)) + b*(d*cos(e + f*x))**m*(a + b*tan(e + f*x))/(f*(1 - m)) + (d*cos(e + f*x))**m*(-a**2*(1 - m) + b**2)*sin(e + f*x)*cos(e + f*x)*hyper((sympy.S.Half, m/2 + sympy.S.Half), (m/2 + sympy.S(3)/2,), cos(e + f*x)**2)/(f*(1 - m)*(m + 1)*sqrt(sin(e + f*x)**2))
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_680():
    f = (d*cos(e + f*x))**m/(a + b*tan(e + f*x))
    F = b*(d*cos(e + f*x))**m*hyper((1, -m/2), (1 - m/2,), b**2*sec(e + f*x)**2/(a**2 + b**2))/(f*m*(a**2 + b**2)) + (d*cos(e + f*x))**m*(sec(e + f*x)**2)**(m/2)*tan(e + f*x)*appellf1(sympy.S.Half, 1, m/2 + 1, sympy.S(3)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(a*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_681():
    f = (d*cos(e + f*x))**m/(a + b*tan(e + f*x))**2
    F = 2*a*b*(d*cos(e + f*x))**m*hyper((2, -m/2), (1 - m/2,), b**2*sec(e + f*x)**2/(a**2 + b**2))/(f*m*(a**2 + b**2)**2) + (d*cos(e + f*x))**m*(sec(e + f*x)**2)**(m/2)*tan(e + f*x)*appellf1(sympy.S.Half, 2, m/2 + 1, sympy.S(3)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(a**2*f) + b**2*(d*cos(e + f*x))**m*(sec(e + f*x)**2)**(m/2)*tan(e + f*x)**3*appellf1(sympy.S(3)/2, 2, m/2 + 1, sympy.S(5)/2, b**2*tan(e + f*x)**2/a**2, -tan(e + f*x)**2)/(3*a**4*f)
    assert integrate(f, x) == F


def test_integrate_4_Trig_functions_4_3_Tangent_4_3_1_2_d_sec_pow_m_a_plus_b_tan_pow_n_682():
    f = (d*cos(e + f*x))**m*(a + b*tan(e + f*x))**n
    F = (d*cos(e + f*x))**m*(1 - (a + b*tan(e + f*x))/(a - sqrt(-b**2)))**(m/2 + 1)*(1 - (a + b*tan(e + f*x))/(a + sqrt(-b**2)))**(m/2 + 1)*(a + b*tan(e + f*x))**(n + 1)*cos(e + f*x)**2*appellf1(n + 1, m/2 + 1, m/2 + 1, n + 2, (a + b*tan(e + f*x))/(a - sqrt(-b**2)), (a + b*tan(e + f*x))/(a + sqrt(-b**2)))/(b*f*(n + 1))
    assert integrate(f, x) == F

